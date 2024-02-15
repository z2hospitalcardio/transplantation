#second heart-----

library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(harmony)
library(scRepertoire)

LVB.data <- Read10X(data.dir = "/home/zeyy/dev/sdd/second_human/LVB")
LVPW.data <- Read10X(data.dir = "/home/zeyy/dev/sdd/second_human/LVPW")
LVB<-CreateSeuratObject(counts = LVB.data, project = 'lvb',min.cells = 3, min.features = 200 )
LVPW<-CreateSeuratObject(counts = LVPW.data, project = 'lvpw',min.cells = 3, min.features = 200 )
heart <- merge(LVB, y = c(LVPW), add.cell.ids = c("lvb", "lvpw"), project = "heart")
head(colnames(heart))

#harmony-----
heart <- NormalizeData(heart)
heart <- FindVariableFeatures(heart, selection.method = "vst", nfeatures = 2000)
heart <- ScaleData(heart, features = rownames(heart))
heart <- RunPCA(heart, features = VariableFeatures(object = heart))
heart <- RunHarmony(heart,reduction = "pca",group.by.vars = "orig.ident")
heart <- heart %>%  RunUMAP(dims=1:30,reduction="harmony") %>%
  FindNeighbors(dims = 1:30,reduction="harmony") %>% FindClusters(resolution=0.1)
heart<-FindNeighbors(heart,dims = 1:30,reduction="harmony")
heart<-FindClusters(heart,resolution=0.3)
DimPlot(heart,label=T)

#QC----
heart[["percent.mt"]] <- PercentageFeatureSet(heart, pattern = "^MT-")
ribo.genes <- grep("^RPS|^RPL", rownames(heart), value = TRUE)
heart[["percent.rb"]] <- PercentageFeatureSet(heart, features = ribo.genes)
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- CaseMatch(HB.genes_total, rownames(heart))
heart[["percent.hb"]] <- PercentageFeatureSet(heart, features=HB_m)
stress.genes_total <- read.csv('/home/kxm/pyh/data/stress.genes.csv')
stress.genes_total <- stress.genes_total$gene
stress_m <- CaseMatch(stress.genes_total, rownames(heart))
heart[["percent.stress"]] <- PercentageFeatureSet(heart, features=stress_m)

VlnPlot(heart,features='nCount_RNA')
VlnPlot(heart,features='nFeature_RNA')
VlnPlot(heart,features='percent.mt')
FeatureScatter(heart, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(heart, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#annotate-------
FeaturePlot(heart,features = 'CD8A')
heart.markers <- FindAllMarkers(heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(heart.markers,'/home/zeyy/dev/sdd/gly/heartClusterMarker.csv')

heart<-subset(heart,subset=seurat_clusters!='7')
new.cluster.ids <- c("T",# 0
                     "Myeloid",# 1
                     "T",# 2
                     "T",# 3
                     "T",# 4
                     "B",# 5
                     "T", # 6
                     "Myeloid",# 8
                     "Mast",# 9
                     "Myeloid"# 10
                     )
names(new.cluster.ids) <- levels(heart)
heart <- RenameIdents(heart, new.cluster.ids)
heart$celltype<-Idents(heart)

DimPlot(heart,label=T)
saveRDS(heart,'/home/zeyy/dev/sdd/pyh/data/second_heart.rds')

#add dmx---------
files<-c('lvb.best','lvpw.best')
second_heart<-readRDS('/home/zeyy/dev/sdd/pyh/data/second_heart.rds')
scelist<-SplitObject(second_heart,split.by='orig.ident')
pos<-c('lvb','lvpw')
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

second_heart_dmx <- merge(scelist[[1]], y = scelist[[2]])
metadata_info <- second_heart_dmx@meta.data[, c('cell.quality','cell.orig','rd.total','rd.pass','rd.uniq','n.snp')]
second_heart <- AddMetaData(object = second_heart, metadata = metadata_info)
DimPlot(second_heart,group.by='cell.orig')
second_heart$cell.orig = unlist(lapply(second_heart$cell.orig, function(x) {
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
DimPlot(second_heart,group.by='cell.orig')
saveRDS(second_heart,'/home/zeyy/dev/sdd/pyh/data/second_heart.rds')

second_heart_SNG<-subset(second_heart,subset=cell.quality=='SNG')
pdf('/home/zeyy/dev/sdd/pyh/result/second_heart/cell_orig.pdf',width=8)
DimPlot(second_heart_SNG,group.by='cell.orig',cols = c("#ED855F","#548F9E"))+NoAxes()
dev.off()

#sub cluster-------
second_heart<-FindNeighbors(second_heart,dims = 1:30,reduction="harmony")
second_heart<-FindClusters(second_heart,resolution=1)
heart.markers <- FindAllMarkers(second_heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(heart.markers,'/home/zeyy/dev/sdd/pyh/result/second_heart/heartClusterMarker1.csv')

DimPlot(second_heart,label=T)
Idents(second_heart)<-'RNA_snn_res.1.5'

FeaturePlot(second_heart,features='GNLY')

new.cluster.ids <- c("CD4/CD8_Temra",# 0
                     "CD8_Tem/emra",# 1
                     "CD8_Tem",# 2
                     "CD4_Tn",# 3
                     "CD8_Tem",# 4
                     "CD4/CD8_Tem",# 5
                     "Myeloid", # 6
                     "CD4_Temra",# 7
                     "CD8_Tem",# 8
                     "CD8_Tn",# 9
                     "NK",# 10
                     "CD8_Tem",# 11
                     "NK",# 12
                     "Myeloid",# 13
                     "CD4_Tn",# 14
                     "CD4_Tn",# 15
                     "B",# 16
                     "Myeloid",# 17
                     "CD8_T_MKI67",# 18
                     "NK",# 19
                     "Myeloid",# 20
                     "Myeloid",# 21
                     "Plasma",# 22
                     "Mast",# 23
                     "Myeloid"# 24
)
names(new.cluster.ids) <- levels(second_heart)
second_heart <- RenameIdents(second_heart, new.cluster.ids)
second_heart$celltype<-Idents(second_heart)

heart.markers <- FindAllMarkers(second_heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(heart.markers,'/home/zeyy/dev/sdd/pyh/result/second_heart/heartClusterMarker3.csv')

#umap-----
lovely_color<-c("#00468BFF","#2673BF","#4DBBD5FF","#00A087FF","#74C476","#EBB365","#ED9A7C","#B15928","#8F0021","orange","red","purple","yellow")

pdf("/home/zeyy/dev/sdd/pyh/result/second_heart/annotation.pdf",width=9)
DimPlot(second_heart, reduction = "umap", label = T, pt.size = 0.5,label.size=5.2,cols=lovely_color)+
  NoAxes()+theme(legend.text = element_text(size = 18))
dev.off()

pdf("/home/zeyy/dev/sdd/pyh/result/second_heart/hyper.pdf",width=10)
DimPlot(second_heart, reduction = "umap",group.by='CTaa', label = F, pt.size = 0.5)+NoAxes()
dev.off()



