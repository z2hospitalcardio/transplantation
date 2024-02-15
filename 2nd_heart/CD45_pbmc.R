#second pbmc-----

library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(harmony)
library(scRepertoire)

pbmc.data <- Read10X(data.dir = "/home/zeyy/dev/sdd/second_human/PBMC")
oldw <- getOption("warn")
options(warn = -1)
pbmc<-CreateSeuratObject(counts = pbmc.data, project = 'pbmc',min.cells = 3, min.features = 200 )
options(warn = oldw)

#QC----
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
ribo.genes <- grep("^RPS|^RPL", rownames(pbmc), value = TRUE)
pbmc[["percent.rb"]] <- PercentageFeatureSet(pbmc, features = ribo.genes)
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- CaseMatch(HB.genes_total, rownames(pbmc))
pbmc[["percent.hb"]] <- PercentageFeatureSet(pbmc, features=HB_m)
stress.genes_total <- read.csv('/home/kxm/pyh/data/stress.genes.csv')
stress.genes_total <- stress.genes_total$gene
stress_m <- CaseMatch(stress.genes_total, rownames(pbmc))
pbmc[["percent.hb"]] <- PercentageFeatureSet(pbmc, features=stress_m)

VlnPlot(pbmc,features='nCount_RNA')
VlnPlot(pbmc,features='nFeature_RNA')
VlnPlot(pbmc,features='percent.mt')
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pbmc <- subset(pbmc, subset = nFeature_RNA>100)


#reduction-----
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.1)
pbmc <- RunUMAP(pbmc, dims = 1:30)
DimPlot(pbmc,label=T)

#annotate-----
FeaturePlot(pbmc,features='LILRA4')

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers,'/home/zeyy/dev/sdd/pyh/result/second_pbmc/ClusterMarker.csv')

pbmc<-subset(pbmc,subset=seurat_clusters!='12')
new.cluster.ids <- c("Neutrophil",# 0
                     "T",# 1
                     "NK",# 2
                     "T",# 3
                     "T",# 4
                     "Macro/Mono",# 5
                     "B", # 6
                     "T",# 7
                     "T",# 8
                     "T",# 9
                     "T",# 10
                     "DC"# 11
                     )
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$celltype<-Idents(pbmc)
saveRDS(pbmc,'/home/zeyy/dev/sdd/pyh/data/second_pbmc.rds')


#add dmx
demuxlet <- fread("/home/zeyy/dev/sdd/demuxlet/second/pbmc.best")
demuxlet<-na.omit(demuxlet)
dmx_new <- data.frame(barcode=demuxlet$BARCODE,
                      cell.quality=sapply(demuxlet$BEST,function(x){strsplit(x,"-")[[1]][[1]]}),
                      cell.orig=demuxlet$SNG.1ST,
                      rd.total=demuxlet$RD.TOTL,
                      rd.pass=demuxlet$RD.PASS,
                      rd.uniq=demuxlet$RD.UNIQ,
                      n.snp=demuxlet$N.SNP)
row.names(dmx_new) <- dmx_new$barcode
dmx_new$barcode <- NULL
second_pbmc<-AddMetaData(second_pbmc,dmx_new)

second_pbmc$cell.orig = unlist(lapply(second_pbmc$cell.orig, function(x) {
  if (is.na(x)) {
    return(NA)
  } else if (x == 'blood_2') {
    return('Recipient')
  } else if (x == 'donor_2' | x == 'host_2') {
    return('Donor')
  } else {
    return(x)
  }
}))

DimPlot(second_pbmc,group.by='cell.orig')

saveRDS(second_pbmc,'/home/zeyy/dev/sdd/pyh/data/second_pbmc.rds')


#sub cluster-------
second_pbmc<-FindNeighbors(second_pbmc,dims = 1:30)
second_pbmc<-FindClusters(second_pbmc,resolution=0.8)
pbmc.markers <- FindAllMarkers(second_pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers,'/home/zeyy/dev/sdd/pyh/result/second_pbmc/pbmcClusterMarker0.8.csv')

DimPlot(second_pbmc,label=T)+NoAxes()
Idents(second_pbmc)<-'RNA_snn_res.1.5'


FeaturePlot(second_pbmc,features='LEF1')+NoAxes()

second_pbmc<-subset(second_pbmc,subset=seurat_clusters!='9')

new.cluster.ids <- c("CD8_Temra",# 0
                     "Neutrophil",# 1
                     "CD4_Tn",# 2
                     "NK",# 3
                     "NK",# 4
                     "Myeloid",# 5
                     "Neutrophil", # 6
                     "CD8_Temra",# 7
                     "CD4_Temra",# 8
                     "CD4_Temra",# 10
                     "CD4_Temra",# 11
                     "CD8_Temra",# 12
                     "CD4_Tem/emra",# 13
                     "B",# 14
                     "Plasma",# 15
                     "CD8_Temra",# 16
                     "CD8_Tn",# 17
                     "CD4_Temra",# 18
                     "Myeloid",# 19
                     "Myeloid",# 20
                     "Myeloid"# 21
                     )
names(new.cluster.ids) <- levels(second_pbmc)
second_pbmc <- RenameIdents(second_pbmc, new.cluster.ids)
second_pbmc$subtype<-Idents(second_pbmc)

pbmc.markers <- FindAllMarkers(second_pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers,'/home/zeyy/dev/sdd/pyh/result/second_pbmc/pbmcClusterMarker3.csv')

#celltype
new.cluster.ids <- c("CD8_T",#CD8_Temra
                     "Neutrophil",# Neutrophil
                     "CD4_T",# CD4_Tn
                     "NK",# NK
                     "Myeloid",# Myeloid
                     "CD4_T",# CD4_Temra
                     "CD4_T",# CD4_Tem/emra
                     "B",# B
                     "Plasma",# Plasma
                     "CD8_T"# CD8_Tn
)
names(new.cluster.ids) <- levels(second_pbmc)
second_pbmc <- RenameIdents(second_pbmc, new.cluster.ids)
second_pbmc$celltype<-Idents(second_pbmc)


#umap-----
lovely_color<-c("#00468BFF","#2673BF","#4DBBD5FF","#00A087FF","#74C476","#EBB365","#ED9A7C","#B15928","#8F0021","orange","red","purple","yellow")

levels(second_pbmc$subtype)<-c("CD4_Tn","CD4_Temra","CD4_Tem/emra","CD8_Tn","CD8_Temra",
                               "NK","Myeloid","Neutrophil","B","Plasma")
Idents(second_pbmc)<-'subtype'

pdf("/home/zeyy/dev/sdd/pyh/result/second_pbmc/annotation.pdf",width=8)
DimPlot(second_pbmc, reduction = "umap", label = T, pt.size = 0.5,label.size=5.2,cols=lovely_color)+
  NoAxes()+theme(legend.text = element_text(size = 18))
dev.off()

pdf("/home/zeyy/dev/sdd/pyh/result/second_pbmc/hyper.pdf",width=10)
DimPlot(second_pbmc, reduction = "umap",group.by='CTaa', label = F, pt.size = 0.5)+NoAxes()
dev.off()


#T--------
Tcell<-subset(pbmc,ident='T')
Tcell <- NormalizeData(Tcell)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell, features = rownames(Tcell))
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
ElbowPlot(Tcell)
Tcell <- FindNeighbors(Tcell, dims = 1:30)
Tcell <- FindClusters(Tcell, resolution = 0.5)
Tcell <- RunUMAP(Tcell, dims = 1:30)
DimPlot(Tcell,label=T)
T.markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(T.markers,'/home/zeyy/dev/sdd/pyh/result/second_pbmc/TClusterMarker.csv')

Tcell<-subset(Tcell,subset=seurat_clusters!='4')
new.cluster.ids <- c("CD8_Temra/eff",# 0
                     "CD4_Tn",# 1
                     "CD4_Temra/eff",# 2
                     "CD8_T",# 3
                     "CD4_T",# 5
                     "CD8_Tem", # 6
                     "CD4_T",# 7
                     "γδT",# 8
                     "CD8_Tn",# 9
                     "CD8_T",# 10
                     "CD4_T",# 11
                     "CD4_T",# 12
                     "CD4_T"# 13
                     )
names(new.cluster.ids) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, new.cluster.ids)
Tcell$celltype<-Idents(Tcell)
DimPlot(Tcell,label=T)

saveRDS(Tcell,'/home/zeyy/dev/sdd/pyh/data/second_pbmcT.rds')

