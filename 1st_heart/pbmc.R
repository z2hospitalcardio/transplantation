library(scRepertoire)
library(Seurat)
library(SeuratObject)
library(tidyverse)


##########step1.add TCR information & CreateSeuratObject###########
S1 <- read.csv("/home/zeyy/dev/sdd/blood/tcr/filtered_contig_annotations.csv")
contig_list <- list(S1)

combined <- combineTCR(contig_list, 
                       samples ='S1', 
                       ID = 'pbmc')
sce=CreateSeuratObject(counts = Read10X('/home/zeyy/dev/sdd/blood/matrix'), 
                       project = 'pbmc' )


##########step2.pbmc QC reduction annotation##############
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
sce[["percent.rb"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
VlnPlot(sce, features = c("nFeature_RNA", "percent.rb", "percent.mt"), ncol = 3)
#filter
sce <- subset(sce, subset = nFeature_RNA > 500 & 
                nFeature_RNA < 5000 & 
                percent.mt < 20)
#normalize
sce <- NormalizeData(sce)
#high variable genes
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sce)
#scale
sce <- ScaleData(sce, features = all.genes)
#reduction
sce <- RunPCA(sce, npcs = 50)
sce <- FindNeighbors(sce, dims = 1:30)
sce <- FindClusters(sce, resolution = 0.2)
sce <- RunUMAP(sce, dims = 1:30)
head(sce@meta.data)

combined$S1_pbmc$barcode=gsub('S1_pbmc_','',combined$S1_pbmc$barcode)
sce <- combineExpression(combined, sce, 
                         cloneCall="gene", 
                         # group.by = "Sample", 
                         proportion = FALSE, 
                         cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
saveRDS(sce,'/home/zeyy/dev/sdd/gly/pbmc.rds')
library(SeuratDisk)
SaveH5Seurat(sce, filename = "/home/zeyy/dev/sdd/gly/pbmc.h5Seurat")
Convert("/home/zeyy/dev/sdd/gly/pbmc.h5Seurat", dest = "h5ad")

sce=readRDS('/home/zeyy/dev/sdd/gly/pbmc_scrublet.rds')
table(sce$doublet_scores)
DimPlot(sce,group.by = 'RNA_snn_res.0.8')
DimPlot(sce,label = T)
sce=FindClusters(sce,resolution = 0.4)
FeaturePlot(sce,'ZNF683')
sc.markers=FindAllMarkers(sce)
write.csv(sc.markers,'/home/zeyy/dev/sdd/gly/blood/pbmc_res0.4_markers.csv')
sce$doublet_scores
FeaturePlot(sce,'CD1C')
#res0.8
new.cluster.ids<-c('CD4+_T',#0
                   'CD4+_T',#1
                   'CD8+_T',#2
                   'Mono/Macro',#3
                   'B',#4
                   'NKT',#5 ZNF683 NKG7 CD2 CD3D 
                   'NK',#6
                   'B',#7
                   'Erythrocyte',#8
                   'Th17',#9
                   'CD8+_T',#10
                   'CD8+_T',#11
                   'CD8+_T',#12
                   'CD8+_T',#13 
                   'Plt',#14
                   'Endo_Dbt',#15
                   'CD8+_T',#16
                   'Low_quality',#17
                   'CD8+_T',#18
                   'DC',#19
                   'Plt'#20
)
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
DimPlot(sce)
sce$celltype=sce@active.ident
sce=subset(sce,celltype_subset %notin% c('Low_quality','Endo_Dbt','Plt')) %>% RunUMAP(dims=1:30)
saveRDS(sce,'/home/zeyy/dev/sdd/gly/blood/pbmc_anno.rds')
