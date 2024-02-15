library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(harmony)

plan("multiprocess", workers = 10)
options(future.globals.maxSize = 8000 * 1024^3)

load("~/pyh/Heart/trans_heart.Rdata")
#cells_to_extract<-c("ap","la","lvib","lvna","ra","raf","rv","sp")
#Heart_new<-subset(Heart_filt,subset=orig.ident %in% cells_to_extract)
#normalization
trans_heart <- NormalizeData(trans_heart)
trans_heart <- FindVariableFeatures(trans_heart, selection.method = "vst", nfeatures = 2000)
trans_heart <- ScaleData(trans_heart, features = rownames(trans_heart))
#reduction
trans_heart <- RunPCA(trans_heart, features = VariableFeatures(object = trans_heart))
trans_heart <- RunHarmony(trans_heart,reduction = "pca",group.by.vars = "orig.ident")
pc.num=1:30
trans_heart <- trans_heart %>%  RunUMAP(dims=pc.num,reduction="harmony") %>%
  FindNeighbors(dims = pc.num,reduction="harmony") %>% FindClusters(resolution=0.1)
#umap
p1<-DimPlot(trans_heart, reduction = "umap", label = T) + NoLegend()
#find markers
trans_heart.markers <- FindAllMarkers(trans_heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(trans_heart.markers,'pyh/result/ClusterMarker.csv')

top10.markers <- trans_heart.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(top10.markers,'pyh/result/ClusterMarker_top.csv')

#annotation------
new.cluster.ids <- c("Endothelial","Fibroblast","SMC","Myeloid","Endothelial", "NK", 
                     "T", "Cardiomyocytes","Neutrophil", "Endothelial","Neuronal","Neuronal")
names(new.cluster.ids) <- levels(trans_heart)
trans_heart <- RenameIdents(trans_heart, new.cluster.ids)
trans_heart$celltype <- Idents(trans_heart)
lovely_color<-c("#00468BFF","#2673BF","#4DBBD5FF","#EBB365","#ED9A7C","#B15928","#00A087FF","#8F0021","#74C476")

pdf("pyh/result/annotated.pdf",width=10)
DimPlot(trans_heart, reduction = "umap", label = FALSE, pt.size = 0.5,label.size=5.2,cols=lovely_color)+
  NoAxes()+theme(legend.text = element_text(size = 18))
dev.off()

dir.create("data")
saveRDS(trans_heart,'~/pyh/data/trans_heart_annotationed.rds')


#drop doublet cells------
library(DoubletFinder)
seu_list=SplitObject(trans_heart,split.by = 'orig.ident')
seu_list=lapply(seu_list,function(trans_heart){
  trans_heart <- NormalizeData(trans_heart)
  trans_heart <- FindVariableFeatures(trans_heart, selection.method = "vst", nfeatures = 2000)
  trans_heart <- ScaleData(trans_heart, features = rownames(trans_heart))
  trans_heart <- RunPCA(trans_heart, features = VariableFeatures(object = trans_heart),reduction.name = "pca")
  trans_heart <- RunUMAP(trans_heart,reduction = "pca", dims = 1:30)
  trans_heart=FindNeighbors(trans_heart,reduction = "pca", dims = 1:30)
  trans_heart=FindClusters(trans_heart,resolution=0.1)
  sweep.res.list <- paramSweep_v3(trans_heart, PCs = 1:30, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  DoubletRate = ncol(trans_heart)*8*1e-6 
  homotypic.prop <- modelHomotypic(trans_heart$seurat_clusters)  
  nExp_poi <- round(DoubletRate*length(trans_heart$seurat_clusters)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  trans_heart <- doubletFinder_v3(trans_heart, PCs = 1:30, pN = 0.25, 
                                  pK = pK_bcmvn, nExp = nExp_poi.adj, 
                                  reuse.pANN = FALSE, sct = F)
})
saveRDS(seu_list,'~/pyh/data/DF.Rdata')


#combine after finddoublet------
seu_list <- lapply(seu_list, function(DF) {
  df_cols <- grep("DF.classifications_0.25", names(DF@meta.data), value = TRUE)
  DF@meta.data$DF.classifications <-DF@meta.data[, df_cols]
  DF@meta.data[, df_cols] <- NULL
  pann_cols<-grep("pANN_0.25", names(DF@meta.data), value = TRUE)
  DF@meta.data$pANN <-DF@meta.data[, pann_cols]
  DF@meta.data[, pann_cols] <- NULL
  return(DF)
})

trans_heart_DF <- Reduce(function(x, y) {merge(x, y)}, seu_list)
#trans_heart_DF <- merge(DF_list[[1]], y = DF_list[2:8])
Idents(trans_heart_DF)<-"celltype"


#change orig.ident name--------
#orig.ident:lvib to lvpw, lvna to lvaw
trans_heart$orig.ident=unlist(lapply(trans_heart$orig.ident,function(x){
  if(x=='lvib'){x='lvpw'}
  else if(x=='lvna'){x='lvaw'}
  else{x=as.character(x)}
}))
trans_heart$orig.ident<-as.factor(trans_heart$orig.ident)
trans_heart_DF$orig.ident=unlist(lapply(trans_heart_DF$orig.ident,function(x){
  if(x=='lvib'){x='lvpw'}
  else if(x=='lvna'){x='lvaw'}
  else{x=as.character(x)}
}))
trans_heart_DF$orig.ident<-as.factor(trans_heart_DF$orig.ident)


#ra_raf differential analysis---------
ra_raf<-subset(trans_heart,subset=orig.ident==c("ra","raf"))
ra_raf$celltype.group <- paste(ra_raf$celltype, ra_raf$orig.ident, sep = "_")
ra_raf$celltype <- Idents(ra_raf)
Idents(ra_raf) <- "celltype.group"
cellfordeg<-levels(ra_raf$celltype)
source("/home/kxm/code/function.R")
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(ra_raf, ident.1 = paste0(cellfordeg[i],"_ra"), ident.2 = paste0(cellfordeg[i],"_raf"), verbose = FALSE,min.pct = 0.1)
  CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
  write.csv(CELLDEG, file = paste0("~/pyh/result/diff/", cellfordeg[i], ".CSV"))
}


#sub celltypes-----

##immune cells-------
sce<-subset(trans_heart,idents=c('Myeloid','NK','T cell','Neutrophil'))
sce=FindNeighbors(sce,reduction = 'harmony')
sce=FindClusters(sce,resolution = 1.2)
#remove low quality & Dbt
sce1=subset(sce1,seurat_clusters %notin% c('11','14','15','16','18','19'))

#Myeloid
Myeloid=subset(sce1,idents = c('1','3','7','9','11','12','13')) 
Myeloid=RunUMAP(Myeloid,reduction = 'harmony',dims = 1:30)
Myeloid=RunUMAP(Myeloid,reduction = 'pca',dims = 1:30,reduction.name = 'umap_naive')
Myeloid=FindNeighbors(Myeloid,reduction = 'harmony',dims = 1:30)
Myeloid=FindClusters(Myeloid,resolution = 1.06)
DimPlot(Myeloid,group.by = 'DF.classifications')
mye2=subset(Myeloid,orig.ident %in% c('ra','raf'))
DimPlot(Myeloid,split.by = 'orig.ident',ncol = 3)
DimPlot(Myeloid,reduction = 'umap_naive',group.by = 'orig.ident')
FeaturePlot(myeloid,'TREM2')

DimPlot(Myeloid,group.by = 'RNA_snn_res.2',label = T)
raf=subset(Myeloid,orig.ident=='raf')
#raf=RunUMAP(raf,reduction = 'pca',dims = 1:30)
DimPlot(Myeloid,label = T)
FeaturePlot(raf,'ABL2')
sce.markers<-FindAllMarkers(Myeloid)
write.csv(sce.markers,'/home/yc/gly/data/Myeloid/Myeloid_markers_res1.8.csv')
Idents(Myeloid)<-'RNA_snn_res.1'
sce.markers<-FindAllMarkers(Myeloid)
write.csv(sce.markers,'/home/yc/gly/data/Myeloid/Myeloid_markers_res1.csv')
DimPlot(Myeloid,label = T)
saveRDS(Myeloid,'/home/zeyy/dev/sdd/gly/Myeloid/Myeloid.rds')
DimPlot(myeloid)
new.cluster.ids <- c("Macro_stressed",#0
                     "Macro_GPR183",#1
                     "Mono_CD16",#2
                     "Macro_FOLR2",#3
                     "Mono_CD16_CD14",#4 
                     "Macro_MARCO", #5
                     "Macro_MARCO", #6
                     "Low_quality",#7
                     "Macro_GPR183", #8
                     "Mono_CD14",#19
                     "DC_CD1C",#10
                     "Macro_APOE")#11
names(new.cluster.ids) <- levels(myeloid)
myeloid <- RenameIdents(myeloid, new.cluster.ids)
myeloid$celltype_subset=unlist(lapply(as.character(myeloid@active.ident),function(x){
  if(x=='Low_quality'){x='Macro_low_quality'}
  else{x=x}
}))
myeloid$celltype_subset=as.factor(myeloid$celltype_subset)
Idents(myeloid)<-'celltype_subset'
sc.markers=FindAllMarkers(myeloid)
write.csv(sc.markers,'/home/yc/gly/data/Myeloid/myeloid_markers.csv')

#T/NK cell
sce1=subset(sce1,seurat_clusters %notin% c('11','14','15','16','18','19'))
sce1=RunUMAP(sce1,reduction = 'harmony',dims = 1:30)
sce1=FindNeighbors(sce1,reduction = 'harmony',dims = 1:30)
sce1=FindClusters(sce1,resolution = 1)
DimPlot(sce1,label = T)
FeaturePlot(Tcell,'KLRB1')
Tcell=subset(sce1,idents = c('4','8','6','10'))
Tcell=RunUMAP(Tcell,reduction = 'harmony',dims = 1:30)
Tcell=FindNeighbors(Tcell,reduction = 'harmony',dims = 1:30)
Tcell=FindClusters(Tcell,resolution = 1.2)
DimPlot(Tcell,label = T,group.by = 'RNA_snn_res.1.2')
Idents(Tcell)<-'RNA_snn_res.1.2'
new.cluster.ids<-c('CD8+_Tem_GZMK',#0 CD8A CD8B GZMK 
                   'CD8+_Temra_CX3CR1',#1 CX3CR1 GZMB GZMH KLRD1 NKG7
                   'CD8+_Tem_GZMK',#2
                   'CD4+_Tem_GZMK',#3 CD4 CD40LG GZMK CD69 CCL5
                   'CD8+_Temra_CX3CR1',#4
                   'CD4+_Tn_LEF1',#5 LEF1 TCF7 SELL IL7R TRAC
                   'T_stressed',#6
                   'T_stressed',#7
                   'CD4+_Tem_GZMK',#8.
                   'T_low_quality',#9
                   'Mast_cell_HDC'#10
) 
names(new.cluster.ids) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, new.cluster.ids)
DimPlot(Tcell)
sce.markers<-FindAllMarkers(Tcell)
write.csv(sce.markers,'/home/zeyy/dev/sdd/gly/Tcell/Tcell_markers.csv')

##endothelial cells------
endothelial=subset(trans_heart,idents = 'Endothelial')
endothelial=RunUMAP(endothelial,dims = 1:30,reduction = 'harmony') %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution=1)
#rerun
{
  endothelial <- NormalizeData(endothelial)
  endothelial <- FindVariableFeatures(endothelial, selection.method = "vst", nfeatures = 2000)
  endothelial <- ScaleData(endothelial, features = rownames(endothelial))
  endothelial <- RunPCA(endothelial, features = VariableFeatures(object = endothelial),reduction.name = "pca")
  
  #endothelial <- RunUMAP(endothelial,reduction = "pca", dims = 1:30, reduction.name = "umap_naive")
  endothelial <- RunHarmony(endothelial,reduction = "pca",group.by.vars = "orig.ident")
  
  endothelial <- RunUMAP(endothelial,reduction = "harmony", dims = 1:30)
  
  endothelial=FindNeighbors(endothelial,reduction = "harmony", dims = 1:30)
  endothelial=FindClusters(endothelial,resolution=1)
  }
endothelial=AddModuleScore(endothelial,features= list(stress_genes$gene),name = 'Stress.score')
FeaturePlot(endothelial,'percent.mt')
DimPlot(endothelial,label = T)
DimPlot(endothelial,group.by = 'cell.orig',split.by = 'orig.ident',ncol = 4)
endothelial=FindClusters(endothelial,resolution = 2)
sc.markers=FindAllMarkers(endothelial)
write.csv(sc.markers,'/home/zeyy/dev/sdd/gly/Endo/Endo_scmarkers_res0.6.csv')
FeaturePlot(endothelial,'ACKR1')
Idents(endothelial)<-'RNA_snn_res.0.6'
new.cluster.ids<-c('EC_vein_NR2F2',#0 ACKR1 NR2F2 TGFBR2 PLVAP CD74
                   'EC_cap_RGCC',#1 RGCC CD36 CA4 LPL CD300LG
                   'EC_art_SEMA3G',#2 SEMA3G VEGFA DEPP1 HEY1 CLIC3
                   'EC_cap_Stressed',#3 DNAJB1  EGR1 APOLD1 ATF3 HSPA1B
                   'EC_low_quality',#4
                   'EC_cap_SPARC',#5 SPARC COL4A1 COL4A2 CTHRC1 TNFRSF4
                   'EC_art_IGFBP3',#6 IGFBP3 SERPINE2 FBLN5 GJA5 PCSK5
                   'Pericyte_RGS5',#7 RGS5 ACTA2 AGT ABCC9 CPE
                   'EC_vein_POSTN',#8 POSTN VCAN COL3A1 CDH11 C7
                   'EC_cap_IFIT1',#9 IFIT1 HMOX1 SPHK1 KLHL21 IFIH1
                   'EC_Fib_Dbt',#10
                   'Dbt_&_lowquality',#11
                   'Dbt_&_lowquality'#12
                   
)
names(new.cluster.ids) <- levels(endothelial)
endothelial <- RenameIdents(endothelial, new.cluster.ids)
endothelial$celltype_subset=endothelial@active.ident
table(endothelial$celltype_subset)
endothelial$celltype_subset=factor(endothelial$celltype_subset,levels=c('EC_art_SEMA3G',
                                                                        'EC_art_IGFBP3',
                                                                        'EC_cap_RGCC',
                                                                        'EC_cap_SPARC',
                                                                        'EC_cap_IFIT1',
                                                                        'EC_cap_Stressed',
                                                                        'EC_vein_NR2F2',
                                                                        'EC_vein_POSTN',
                                                                        'Pericyte_RGS5',
                                                                        'EC_low_quality',
                                                                        'EC_Fib_Dbt',
                                                                        'Dbt_&_lowquality'
))
Idents(endothelial)<-'celltype_subset'
DimPlot(endothelial)
endothelial=subset(endothelial,celltype_subset %notin% c('EC_low_quality',
                                                         'EC_Fib_Dbt',
                                                         'Dbt_&_lowquality'))
endothelial=RunUMAP(endothelial,dims = 1:30,reduction = 'harmony')
saveRDS(endothelial,'/home/zeyy/dev/sdd/gly/Endo/Endo.rds')
DimPlot(endothelial,split.by='orig.ident',group.by='cell.orig',ncol=4)

