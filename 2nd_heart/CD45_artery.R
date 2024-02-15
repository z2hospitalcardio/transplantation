#second_artery
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(harmony)

CX.data <- Read10X(data.dir = "/home/zeyy/dev/sdd/second_human/CX")
LAD.data <- Read10X(data.dir = "/home/zeyy/dev/sdd/second_human/LAD")
options(warn = 0)
CX<-CreateSeuratObject(counts = CX.data, project = 'cx',min.cells = 3, min.features = 200 )
LAD<-CreateSeuratObject(counts = LAD.data, project = 'lad',min.cells = 3, min.features = 200 )
artery <- merge(CX, y = c(LAD), add.cell.ids = c("cx", "lad"), project = "artery")
head(colnames(artery))

#harmony-----
artery <- NormalizeData(artery)
artery <- FindVariableFeatures(artery, selection.method = "vst", nfeatures = 2000)
artery <- ScaleData(artery, features = rownames(artery))
artery <- RunPCA(artery, features = VariableFeatures(object = artery))
artery <- RunHarmony(artery,reduction = "pca",group.by.vars = "orig.ident")
artery <- artery %>%  RunUMAP(dims=1:30,reduction="harmony") %>%
  FindNeighbors(dims = 1:30,reduction="harmony") %>% FindClusters(resolution=0.1)
DimPlot(artery,label=T)

#QC----
artery[["percent.mt"]] <- PercentageFeatureSet(artery, pattern = "^MT-")
ribo.genes <- grep("^RPS|^RPL", rownames(artery), value = TRUE)
artery[["percent.rb"]] <- PercentageFeatureSet(artery, features = ribo.genes)
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- CaseMatch(HB.genes_total, rownames(artery))
artery[["percent.hb"]] <- PercentageFeatureSet(artery, features=HB_m)
stress.genes_total <- read.csv('/home/kxm/pyh/data/stress.genes.csv')
stress.genes_total <- stress.genes_total$gene
stress_m <- CaseMatch(stress.genes_total, rownames(artery))
artery[["percent.hb"]] <- PercentageFeatureSet(artery, features=stress_m)

VlnPlot(artery,features='nCount_RNA')
VlnPlot(artery,features='nFeature_RNA')
VlnPlot(artery,features='percent.mt')
FeatureScatter(artery, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(artery, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#annotate-------
FeaturePlot(artery,features = 'CD4')
artery.markers <- FindAllMarkers(artery, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(artery.markers,'/home/zeyy/dev/sdd/pyh/result/second_artery/arteryClusterMarker.csv')

new.cluster.ids <- c("T",# 0
                     "Myeloid",# 1
                     "T",# 2
                     "NK",# 3
                     "B",# 4
                     "T",# 5
                     "T", # 6
                     "Plasma"# 7
)
names(new.cluster.ids) <- levels(artery)
artery <- RenameIdents(artery, new.cluster.ids)
artery$celltype<-Idents(artery)
saveRDS(artery,'/home/zeyy/dev/sdd/pyh/data/second_artery.rds')

#demuxlet--------
setwd("/home/zeyy/dev/sdd/demuxlet/second/")
files<-c('scx.best','lad.best')
second_artery<-readRDS('/home/zeyy/dev/sdd/pyh/data/second_artery.rds')
scelist<-SplitObject(second_artery,split.by='orig.ident')
pos<-c('cx','lad')
for (i in 1:2) {
  demuxlet <- fread(files[i])
  demuxlet$BARCODE<-paste0(pos,'_',demuxlet$BARCODE)
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
  scelist[[i]]<-AddMetaData(scelist[[i]],dmx_new)
}

second_artery_dmx <- merge(scelist[[1]], y = scelist[[2]])
metadata_info <- second_artery_dmx@meta.data[, c('cell.quality','cell.orig','rd.total','rd.pass','rd.uniq','n.snp')]
second_artery <- AddMetaData(object = second_artery, metadata = metadata_info)
DimPlot(second_artery,group.by='cell.orig')
second_artery$cell.orig = unlist(lapply(second_artery$cell.orig, function(x) {
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

second_artery_SNG<-subset(second_artery,subset=cell.quality=='SNG')


#sub clusters-----
second_artery<-FindNeighbors(second_artery,dims = 1:30,reduction="harmony")
second_artery<-FindClusters(second_artery,resolution=2.5)
DimPlot(second_artery,label=T)


FeaturePlot(second_artery,features='percent.mt')
second_artery<-subset(second_artery,subset=seurat_clusters!='22')

new.cluster.ids <- c("CD8_Tem",# 0
                     "CD4_Tn",# 1
                     "CD8_Tem",# 2
                     "CD4_Tn",# 3
                     "CD8_Tem",# 4
                     "CD4_Tem",# 5
                     "Tcm", # 6
                     "CD4_Tn",# 7
                     "CD8_Tn",# 8
                     "CD8_Tem",# 9
                     "Myeloid",# 10
                     "CD4_Tn",# 11
                     "Myeloid",# 12
                     "B",# 13
                     "Myeloid",# 14
                     "NK",# 15
                     "CD4_Treg",# 16
                     "Myeloid",# 17
                     "NK",# 18
                     "CD8_T_MKI67",# 19
                     "Myeloid",# 20
                     "Plasma"# 21
)
names(new.cluster.ids) <- levels(second_artery)
second_artery <- RenameIdents(second_artery, new.cluster.ids)
second_artery$celltype<-Idents(second_artery)

artery.markers <- FindAllMarkers(second_artery, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(artery.markers,'/home/zeyy/dev/sdd/pyh/result/second_artery/arteryClusterMarker3.csv')

new.cluster.ids <- c("CD4_Tn",# 0
                     "CD8_Tem",# 1
                     "CD4_Tem",# 2
                     "CD8_Tem",# 3
                     "CD8_Tem",# 4
                     "CD4_Tn",# 5
                     "CD4_Tn", # 6
                     "Tcm",# 7
                     "Myeloid",# 8
                     "CD4_Tem",# 9
                     "CD4_Tn",# 10
                     "Myeloid",# 11
                     "B",# 12
                     "NKT",# 13
                     "CD8_Tn",# 14
                     "CD8_Tem",# 15
                     "CD4_Tn",# 16
                     "NK",# 17
                     "Myeloid",# 18
                     "CD4_Treg",# 19
                     "CD8_Tem",# 20
                     "NK",# 21
                     "CD8_Tem",# 22
                     "CD8_Tem",# 23
                     "CD8_Tn",# 24
                     "CD8_Tem",# 25
                     "Myeloid",# 26
                     "CD8_T_MKI67",# 27
                     "Myeloid",# 28
                     "Myeloid",# 29
                     "Plasma",# 30
                     "Myeloid"# 31
)
names(new.cluster.ids) <- levels(second_artery)
second_artery <- RenameIdents(second_artery, new.cluster.ids)
second_artery$celltype2<-Idents(second_artery)

#umap-----
lovely_color<-c("#00468BFF","#2673BF","#4DBBD5FF","#00A087FF","#74C476","#EBB365","#ED9A7C","#B15928","#8F0021","orange","red")

pdf("/home/zeyy/dev/sdd/pyh/result/second_artery/annotation.pdf",width=9)
DimPlot(second_artery, reduction = "umap", label = T, pt.size = 0.5,label.size=5.2,cols=lovely_color)+
  NoAxes()+theme(legend.text = element_text(size = 18))
dev.off()

pdf("/home/zeyy/dev/sdd/pyh/result/second_artery/hyper.pdf",width=10)
DimPlot(second_artery, reduction = "umap",group.by='CTaa', label = F, pt.size = 0.5)+NoAxes()
dev.off()

