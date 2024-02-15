library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(harmony)

plan("multiprocess", workers = 10)
options(future.globals.maxSize = 8000 * 1024^3)

artery<-readRDS("/home/zeyy/dev/sdd/pyh/data/artery_QC.rds")

#/home/zeyy/dev/sdd/artery
rca.data <- Read10X(data.dir = "/home/zeyy/dev/sdd/artery/rca/")
na.data <- Read10X(data.dir = "/home/zeyy/dev/sdd/artery/na/")
ladb.data <- Read10X(data.dir = "/home/zeyy/dev/sdd/artery/ladb/")
cx.data <- Read10X(data.dir = "/home/zeyy/dev/sdd/artery/cx/")
rca<-CreateSeuratObject(counts = rca.data, project = "rca", min.cells = 3, min.features = 200)
na<-CreateSeuratObject(counts = na.data, project = "na", min.cells = 3, min.features = 200)
ladb<-CreateSeuratObject(counts = ladb.data, project = "ladb", min.cells = 3, min.features = 200)
cx<-CreateSeuratObject(counts = cx.data, project = "cx", min.cells = 3, min.features = 200)

#merge---------
artery <- merge(rca, y = c(na,ladb,cx), add.cell.ids = c("rca", "na","ladb","cx"), project = "ALL")
head(colnames(artery))

#harmony--------
artery<-artery%>%
  NormalizeData()%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)%>%
  ScaleData(features=rownames(artery))

artery <- RunPCA(artery, features = VariableFeatures(object = artery))
artery <- RunHarmony(artery,reduction = "pca",group.by.vars = "orig.ident")
pc.num=1:30
#artery <- subset(artery, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 25)
artery <- artery %>%  RunUMAP(dims=pc.num,reduction="harmony") %>%
  FindNeighbors(dims = pc.num,reduction="harmony") %>% FindClusters(resolution=0.3)

DimPlot(artery, reduction = "umap", label = T) + NoLegend()
#find markers
artery.markers <- FindAllMarkers(artery, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(artery.markers,'/home/kxm/pyh/result/ClusterMarker.csv')
FeaturePlot(artery,features="VWF")#Endo
FeaturePlot(artery,features="PECAM1")#Endo
FeaturePlot(artery,features="ECSCR")#Endo
FeaturePlot(artery,features="NKG7")#NK
FeaturePlot(artery,features="IL7R")#T
FeaturePlot(artery,features="CD163")#myeloid
FeaturePlot(artery,features="CD14")#myeloid
FeaturePlot(artery,features="DCN")#Fibro
FeaturePlot(artery,features=c("MPZ","PLP1"))#Neuro
FeaturePlot(artery,features="MYH11")#SMC
FeaturePlot(artery,features="MYH7")#Cardio


#annotate------
new.cluster.ids <- c("Myeloid",# 0
                     "Fibroblast",# 1
                     "Tcell",# 2
                     "Endothelial",# 3
                     "SMC",# 4
                     "Endothelial", # 5
                     "Myofibroblast", # 6
                     "Myeloid",# 7
                     "B/Plasma",# 8
                     "NK",# 9
                     "Endothelial",# 10
                     "Endothelial",# 11
                     "Fibroblast",# 12
                     "Neuronal",# 13
                     "Cardiomyocyte",# 14
                     "drop"# 15
                     )
names(new.cluster.ids) <- levels(artery)
artery <- RenameIdents(artery, new.cluster.ids)
artery$celltype <- Idents(artery)


#QC-------
artery[["percent.mt"]] <- PercentageFeatureSet(artery, pattern = "^MT-")
ribo.genes <- grep("^RPS|^RPL", rownames(artery), value = TRUE)
artery[["percent.rb"]] <- PercentageFeatureSet(artery, features = ribo.genes)
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- CaseMatch(HB.genes_total, rownames(artery))
artery[["percent.hb"]] <- PercentageFeatureSet(artery, features=HB_m)

##doubletfinder
library(DoubletFinder)
seu_list=c(rca,na,ladb,cx)
seu_list=lapply(seu_list,function(artery){
  artery <- NormalizeData(artery)
  artery <- FindVariableFeatures(artery, selection.method = "vst", nfeatures = 2000)
  artery <- ScaleData(artery, features = rownames(artery))
  artery <- RunPCA(artery, features = VariableFeatures(object = artery),reduction.name = "pca")
  artery <- RunUMAP(artery,reduction = "pca", dims = 1:30)
  artery=FindNeighbors(artery,reduction = "pca", dims = 1:30)
  artery=FindClusters(artery,resolution=0.1)
  sweep.res.list <- paramSweep_v3(artery, PCs = 1:30, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  DoubletRate = ncol(artery)*8*1e-6 
  homotypic.prop <- modelHomotypic(artery$seurat_clusters)  
  
  nExp_poi <- round(DoubletRate*length(artery$seurat_clusters)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  artery <- doubletFinder_v3(artery, PCs = 1:30, pN = 0.25, 
                                  pK = pK_bcmvn, nExp = nExp_poi.adj, 
                                  reuse.pANN = FALSE, sct = F)
})
###combine after finddoublet
seu_list <- lapply(seu_list, function(DF) {
  df_cols <- grep("DF.classifications_0.25", names(DF@meta.data), value = TRUE)
  DF@meta.data$DF.classifications <-DF@meta.data[, df_cols]
  DF@meta.data[, df_cols] <- NULL
  pann_cols<-grep("pANN_0.25", names(DF@meta.data), value = TRUE)
  DF@meta.data$pANN <-DF@meta.data[, pann_cols]
  DF@meta.data[, pann_cols] <- NULL
  return(DF)
})
artery_DF <- merge(seu_list[[1]], y = seu_list[2:4], add.cell.ids = c("rca", "na","ladb","cx"), project = "ALL")
artery$DF.classifications<-artery_DF$DF.classifications
DimPlot(artery, reduction = "umap",group.by="DF.classifications")

##scrublet
setwd("/home/zeyy/dev/sdd/artery")
files <- list.files(pattern = "cx|rca|ladb|na")
artery_list <- lapply(files, function(x) {
  data <- Read10X(x)  
  seurat_object <- CreateSeuratObject(data)
  return(seurat_object)
})
setwd("/home/kxm/pyh/result/scrublet_strict")
files <- list.files(pattern = "cx|rca|ladb|na.*doublet_strict\\.csv$")
doublet_list <- lapply(files, function(x) {
  data <- read.csv(x,row.names=1) 
  return(data)
})
for (i in 1:4){
  artery_list[[i]]<-AddMetaData(artery_list[[i]],doublet_list[[i]])
}
artery_df <- merge(artery_list[[1]], y = artery_list[2:4], add.cell.ids = c("cx", "ladb","na","rca"), project = "ALL")
artery$doublet_scores_strict<-artery_df$doublet_scores
artery$predicted_doublets_strict<-artery_df$predicted_doublets


##stress score
df <- read.csv("/home/kxm/pyh/stress.genes.csv")
lst_stress <- list()
lst_stress[[1]] <- df[, 1]
artery <- AddModuleScore(artery,
                       features = lst_stress,
                       ctrl = 100,
                       name = "stress.score")
FeaturePlot(artery,features = 'stress.score')

#drop low quality cells
Idents(artery)<-'RNA_snn_res.0.4'
keep_cells <- !(ident %in% c(15, 17))
artery_drop <- subset(artery, cells = which(keep_cells))
Idents(artery_drop)<-'celltype'
saveRDS(artery_drop,"/home/zeyy/dev/sdd/pyh/data/artery_dropDF.rds")

#DecontX,RNA contamination--------
library(celda)
counts <- artery@assays$RNA@counts
decontX_results <- decontX(counts) 
artery$Contamination =decontX_results$contamination
FeaturePlot(artery, features = 'Contamination', )
low_con_artery = artery[,artery$Contamination < 0.2]


#Dimplot-----
artery$celltype<-factor(artery$celltype,levels=c("Endothelial","Fibroblast","SMC","Myeloid",
                                         "NK",'Tcell',"Cardiomyocyte","B/Plasma","Neuronal","Myofibroblast"))
Idents(artery)<-'celltype'
lovely_color<-c("#00468BFF","#2673BF","#4DBBD5FF","#EBB365","#ED9A7C","#B15928","#00A087FF","#8F0021","#74C476","orange")
pdf("/home/zeyy/dev/sdd/pyh/result/artery/annotated.pdf",width=9)
DimPlot(artery, reduction = "umap", label = T, pt.size = 0.5,label.size=5.2,cols=lovely_color)+
  NoAxes()+theme(legend.text = element_text(size = 18))
dev.off()

artery<-RunUMAP(artery,dims=pc.num,reduction="harmony")
artery<-FindNeighbors(artery,dims = pc.num,reduction="harmony")
artery<-FindClusters(artery,resolution=0.5)
DimPlot(artery,label=T)

artery<-subset(artery,ident='9',invert=T)
  
FeaturePlot(artery,features='doublet_scores_strict')
FeaturePlot(artery,features='MYH11')
DimPlot(artery,group.by='orig.ident')

artery.markers <- FindAllMarkers(artery, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(artery.markers,'/home/kxm/pyh/result/ClusterMarker_artery.csv')

new.cluster.ids <- c("Myeloid",# 0
                     "Endothelial",# 1
                     "Fibroblast",# 2
                     "SMC",# 3
                     "T",# 4
                     "Myeloid", # 5
                     "NK", # 6
                     "Myofibroblast",# 7
                     "Fibroblast",# 8
                     "Endothelial",# 9
                     "Myeloid",# 10
                     "B/Plasma",# 11
                     "Cardiomyocyte",# 12
                     "Neuronal"# 13
)
names(new.cluster.ids) <- levels(artery)
artery <- RenameIdents(artery, new.cluster.ids)
artery$celltype <- Idents(artery)

saveRDS(artery,'/home/zeyy/dev/sdd/pyh/data/artery_reanno.rds')


#Tcell drop doublets-------
artery<-readRDS('/home/zeyy/dev/sdd/pyh/data/artery_dropCont.rds')

Tcell<-subset(artery,ident='Tcell')
Tcell<-Tcell%>%
  NormalizeData()%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)%>%
  ScaleData(features=rownames(Tcell))

Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
Tcell <- RunHarmony(Tcell,reduction = "pca",group.by.vars = "orig.ident")
pc.num=1:30
#Tcell <- subset(Tcell, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 25)
Tcell <- Tcell %>%  RunUMAP(dims=pc.num,reduction="harmony") %>%
  FindNeighbors(dims = pc.num,reduction="harmony") %>% FindClusters(resolution=0.8)
Tcell<-FindNeighbors(Tcell,dims = pc.num,reduction="harmony")
Tcell<-FindClusters(Tcell,resolution=0.8)
DimPlot(Tcell,label=T)

T.markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(T.markers,'/home/kxm/pyh/result/TClusterMarker.csv')

DimPlot(Tcell,group.by='predicted_doublets_strict')
DimPlot(Tcell,group.by='orig.ident')
FeaturePlot(Tcell,features='doublet_scores_strict')
FeaturePlot(Tcell,features='CD74')
bad_cells <- WhichCells(object = Tcell, idents = c("5","6"))

artery_n <- subset(x = artery_n, cells = setdiff(Cells(artery_n), bad_cells))

Tcell <- subset(x = Tcell, cells = setdiff(Cells(Tcell), bad_cells))
saveRDS(Tcell,'/home/zeyy/dev/sdd/pyh/data/artery_T.rds')

#Myeloid drop doublets------
Myeloid<-subset(artery,ident='Myeloid')
Myeloid<-Myeloid%>%
  NormalizeData()%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)%>%
  ScaleData(features=rownames(Myeloid))

Myeloid <- RunPCA(Myeloid, features = VariableFeatures(object = Myeloid))
Myeloid <- RunHarmony(Myeloid,reduction = "pca",group.by.vars = "orig.ident")
pc.num=1:30
#Myeloid <- subset(Myeloid, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 25)
Myeloid <- Myeloid %>%  RunUMAP(dims=pc.num,reduction="harmony") %>%
  FindNeighbors(dims = pc.num,reduction="harmony") %>% FindClusters(resolution=0.8)
Myeloid<-FindNeighbors(Myeloid,dims = pc.num,reduction="harmony")
Myeloid<-FindClusters(Myeloid,resolution=1.5)
DimPlot(Myeloid,label=T)

T.markers <- FindAllMarkers(Myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(T.markers,'/home/kxm/pyh/result/TClusterMarker.csv')

DimPlot(Myeloid,group.by='predicted_doublets_strict')
DimPlot(Myeloid,group.by='orig.ident')
FeaturePlot(Myeloid,features='doublet_scores_strict')
FeaturePlot(Myeloid,features='DCN')

bad_cells <- WhichCells(object = Myeloid, idents = c("0","3","10","25","11","15"))

artery <- subset(x = artery, cells = setdiff(Cells(artery), bad_cells))

Myeloid <- subset(x = Myeloid, cells = setdiff(Cells(Myeloid), bad_cells))
saveRDS(Myeloid,'/home/zeyy/dev/sdd/pyh/data/artery_Myeloid.rds')


#SMC drop doublets--------
SMC<-subset(artery,ident='SMC')
SMC<-SMC%>%
  NormalizeData()%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)%>%
  ScaleData(features=rownames(SMC))

SMC <- RunPCA(SMC, features = VariableFeatures(object = SMC))
SMC <- RunHarmony(SMC,reduction = "pca",group.by.vars = "orig.ident")
pc.num=1:30
#SMC <- subset(SMC, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 25)
SMC <- SMC %>%  RunUMAP(dims=pc.num,reduction="harmony") %>%
  FindNeighbors(dims = pc.num,reduction="harmony") %>% FindClusters(resolution=0.8)
SMC<-FindNeighbors(SMC,dims = pc.num,reduction="harmony")
SMC<-FindClusters(SMC,resolution=1.5)
DimPlot(SMC,label=T)

T.markers <- FindAllMarkers(SMC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(T.markers,'/home/kxm/pyh/result/SMCClusterMarker.csv')

DimPlot(SMC,group.by='predicted_doublets_strict')
DimPlot(SMC,group.by='orig.ident')
FeaturePlot(SMC,features='doublet_scores_strict')
FeaturePlot(SMC,features='ACTA2')
bad_cells <- WhichCells(object = SMC, idents = c("1","12","10","7","9","2"))

artery_n <- subset(x = artery_n, cells = setdiff(Cells(artery_n), bad_cells))

SMC <- subset(x = SMC, cells = setdiff(Cells(SMC), bad_cells))
saveRDS(SMC,'/home/zeyy/dev/sdd/pyh/data/artery_SMC.rds')

#endo drop doublets------
bad_cells<-readRDS('/home/zeyy/dev/sdd/gly/artery/endo/bad_ec.rds')

bad_cells<-WhichCells(bad_cells)

#add dmx metadata-------
artery<-subset(artery_n,subset=orig.ident!='na')
artery_dmx<-readRDS('/home/zeyy/dev/sdd/pyh/data/artery_dxm.rds')
metadata_info <- artery_dmx@meta.data[, c('cell.quality','cell.orig')]
artery <- AddMetaData(object = artery, metadata = metadata_info)
artery$cell.orig[artery$cell.orig == 'Re'] <- 'Recipient'
artery$cell.orig[artery$cell.orig == 'AP'] <- 'Donor'

saveRDS(artery,'/home/zeyy/dev/sdd/pyh/data/artery_dmx.rds')

saveRDS(artery_n,'/home/zeyy/dev/sdd/pyh/data/artery_dropDF.rds')

#artery myeloid cluster-----
artery_myeloid<-subset(artery_n,idents='Myeloid')
artery_myeloid<-artery_myeloid%>%
  NormalizeData()%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)%>%
  ScaleData(features=rownames(artery_myeloid))

artery_myeloid <- RunPCA(artery_myeloid, features = VariableFeatures(object = artery_myeloid))
artery_myeloid <- RunHarmony(artery_myeloid,reduction = "pca",group.by.vars = "orig.ident")
pc.num=1:30
#artery_myeloid <- subset(artery_myeloid, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 25)
artery_myeloid <- artery_myeloid %>%  RunUMAP(dims=pc.num,reduction="harmony") %>%
  FindNeighbors(dims = pc.num,reduction="harmony") %>% FindClusters(resolution=1)
artery_myeloid_0.3<-FindNeighbors(artery_myeloid_0.3,dims = pc.num,reduction="harmony")
artery_myeloid0.3<-FindClusters(artery_myeloid_0.3,resolution=1.5)
DimPlot(artery_myeloid,label=T,group.by='RNA_snn_res.0.3')
DimPlot(artery_myeloid,label=T)

FeaturePlot(artery_myeloid,features='SPP1')

mye.markers <- FindAllMarkers(artery_myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(mye.markers,'/home/kxm/pyh/result/heart_myeClusterMarker.csv')

#res=1
bad_cells <- WhichCells(object = artery_myeloid, idents = '9')
artery_myeloid <- subset(x = artery_myeloid, cells = setdiff(Cells(artery_myeloid), bad_cells))

#res=0.3
bad_cells <- WhichCells(object = artery_myeloid, idents = '0')
artery_myeloid <- subset(x = artery_myeloid, cells = setdiff(Cells(artery_myeloid), bad_cells))

DimPlot(artery_myeloid_1.5)
#annotate
new.cluster.ids <- c("Macro_SPP1",# 0
                     "Macro_FOLR2",# 1
                     "Macro_LVYE1",# 2
                     "Macro_FOLR2",# 3
                     "Mono_CD16",# 4
                     "Macro_SPP1",# 5
                     "Macro_SPP1", # 6
                     "Macro_MKI67",# 7
                     "Macro_FOLR2",# 8
                     "DC_CD1C",# 9
                     "Macro_FOLR2",# 10
                     "Macro_FOLR2",# 11
                     "Mono_S100A8/A9",# 12
                     "Macro_FOLR2",# 13
                     "Macro_FOLR2",# 14
                     "Macro_MKI67",# 15
                     "Mono_CD14",# 16
                     "Mono_CD16",# 17
                     "Macro_SPP1",# 18
                     "DC_CD1C",# 19
                     "Macro_MKI67",# 20
                     "Macro_SPP1",# 21
                     "Macro_SPP1",# 22
                     "Mono_CD16"# 23
)
names(new.cluster.ids) <- levels(artery_myeloid)
artery_myeloid <- RenameIdents(artery_myeloid, new.cluster.ids)
artery_myeloid$celltype_sub <- Idents(artery_myeloid)



