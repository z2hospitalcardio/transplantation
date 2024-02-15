#apex---------
ap.data <- Read10X(data.dir = "/home/zeyy/dev/sdf/human/second_human/scRNA_cellranger/ap_nul/outs/filtered_feature_bc_matrix")
ap<-CreateSeuratObject(counts = ap.data, project = 'ap',min.cells = 3, min.features = 200 )

ap <- NormalizeData(ap)
ap <- FindVariableFeatures(ap, selection.method = "vst", nfeatures = 2000)
ap <- ScaleData(ap, features = rownames(ap))
ap <- RunPCA(ap, features = VariableFeatures(object = ap))
ElbowPlot(ap,ndims=30)

ap <- FindNeighbors(ap, dims = 1:20)
ap <- FindClusters(ap, resolution = 0.4)
ap <- RunUMAP(ap, dims = 1:20)
DimPlot(ap, reduction = "umap",label=T)

FeaturePlot(ap,features = 'DCN')

ap.markers <- FindAllMarkers(ap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ap.markers<-subset(ap.markers,subset=p_val_adj<0.05)
write.csv(ap.markers,'/home/zeyy/dev/sdd/pyh/result/second_heart/12.29/ap_nul_marker_0.4.csv')

ap[["percent.mt"]] <- PercentageFeatureSet(ap, pattern = "^MT-")
ribo.genes <- grep("^RPS|^RPL", rownames(ap), value = TRUE)
ap[["percent.rb"]] <- PercentageFeatureSet(ap, features = ribo.genes)
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- CaseMatch(HB.genes_total, rownames(ap))
ap[["percent.hb"]] <- PercentageFeatureSet(ap, features=HB_m)
stress.genes_total <- read.csv('/home/kxm/pyh/data/stress.genes.csv')
stress.genes_total <- stress.genes_total$gene
stress_m <- CaseMatch(stress.genes_total, rownames(ap))
ap[["percent.stress"]] <- PercentageFeatureSet(ap, features=stress_m)

new.cluster.ids <- c('Cardiomyocyte',#0
                     'Cardiomyocyte',#1
                     'Myeloid',#2
                     'Cardiomyocyte',#3
                     'Endothelial',#4
                     'T',#5
                     'Cardiomyocyte',#6
                     'Pericyte',#7
                     'Cardiomyocyte',#8
                     'B',#9
                     'Fibroblast'#10
)
names(new.cluster.ids) <- levels(ap)
ap <- RenameIdents(ap, new.cluster.ids)
ap$celltype<-Idents(ap)

DimPlot(ap,label=T)

##add dmx----
demuxlet <- fread('/home/zeyy/dev/sdd/demuxlet/second/ap_nul.best')
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
ap<-AddMetaData(ap,dmx_new)

ap$cell.orig = unlist(lapply(ap$cell.orig, function(x) {
  if (is.na(x)) {return(NA)}
  else if (x == 'blood_2') {return('Recipient')} 
  else if (x == 'donor_2' | x == 'host_2') {return('Donor')} else {
    return(x)
  }
}))

DimPlot(ap,group.by='cell.orig')

saveRDS(ap,'/home/zeyy/dev/sdd/pyh/data/second_ap.rds')

##cell.orig,umap-----
DimPlot(ap,group.by='cell.orig')+NoAxes()
ap_SNG<-subset(ap,subset=cell.quality=='SNG')
pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/second/heart_orig_umap.pdf',height = 6)
DimPlot(ap_SNG,group.by='cell.orig',cols = c("#ED855F","#548F9E"))+NoAxes()
dev.off()

##cell.orig,histogram------
tmp <- select(ap_SNG@meta.data, c("celltype","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
celltype<-unique(tmp$celltype)
celltype<-celltype[2:7]
celltype<-droplevels(celltype)
for(i in celltype){
  df_i <- subset(tmp, tmp$celltype==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("celltype", "cell.orig", "value")
  df <- rbind(df, df_i)
}
df$celltype <- factor(df$celltype, levels = c('Endothelial',"Fibroblast",
                                              "Pericyte",'Myeloid',"T","B"))
p <- ggplot(df, aes(x=celltype, y=value, fill=cell.orig)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/second/heart_orig.pdf',width = 6.5)
p
dev.off() 

##differential gene number,histogram--------
df<-data.frame()
for(i in celltype){
  seu_obj<-subset(ap_SNG,ident=i)
  Idents(seu_obj)<-'cell.orig'
  CELLDEG <- FindMarkers(seu_obj, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
  CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
  CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
  RE <- sum(CELLDEG$avg_log2FC > 0)
  DO <- sum(CELLDEG$avg_log2FC < 0)
  df_RE<-data.frame(celltype=i,cell.orig='Recipient',value=RE)
  df_DO<-data.frame(celltype=i,cell.orig='Donor',value=DO)
  df_i<-rbind(df_RE,df_DO)
  df<- rbind(df, df_i)
}
celltype_order<-df%>%group_by(celltype)%>%
  summarise(value = sum(value))%>%arrange(desc(value))
df$celltype <- factor(df$celltype, levels = celltype_order$celltype)

p <- ggplot(df, aes(x=celltype, y=value, fill=cell.orig)) + 
  geom_col() + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  #scale_fill_npg()+
  labs(x = '', y = '') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/second/heart_redo_DEG.pdf',width = 6.5)
p
dev.off()

##cellchat------
ap_drop<-subset(ap_SNG,celltype!='Cardiomyocyte')
metadata_info <- ap_drop@meta.data[, c('celltype','cell.orig')]
#combine cell.orig&celltype
metadata_info$label <- paste(metadata_info$celltype, metadata_info$cell.orig, sep = "_")
data.input <- GetAssayData(ap_drop, assay = "RNA", slot = "data")
data.input <- data.input[,rownames(metadata_info) ]

cellchat_heart <- createCellChat(object = data.input, meta = metadata_info, group.by = "label")

cellchat_heart@meta[["label"]] <- sub("_Recipient$", "_Re", cellchat_heart@meta[["label"]])
cellchat_heart@meta[["label"]] <- sub("_Donor$", "_Do", cellchat_heart@meta[["label"]])
cellchat_heart@idents <- as.factor(cellchat_heart@meta[['label']])

#preprocess
CellChatDB <- CellChatDB.human  
showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
# set the used database in the object
cellchat_heart@DB <- CellChatDB.use

#Ligand-Receptor interaction
cellchat_heart <- subsetData(cellchat_heart) 
future::plan("multisession", workers = 4) # do parallel
cellchat_heart <- identifyOverExpressedGenes(cellchat_heart)
cellchat_heart <- identifyOverExpressedInteractions(cellchat_heart)

#calculate probabilities
cellchat_heart <- computeCommunProb(cellchat_heart,population.size = TRUE)
cellchat_heart <- filterCommunication(cellchat_heart, min.cells = 10)
cellchat_heart <- computeCommunProbPathway(cellchat_heart)

cellchat_heart <- aggregateNet(cellchat_heart)

#visualization
groupSize <- as.numeric(table(cellchat_heart@idents))
pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/second/cellchat/heart_number.pdf')
netVisual_circle(cellchat_heart@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/second/cellchat/heart_weight.pdf')
netVisual_circle(cellchat_heart@net$weight, vertex.weight = groupSize, idents.use = c(),
                 weight.scale = T,  title.name = "Interaction weights/strength")
dev.off()

###LR interaction----
weight_matrix <- cellchat_heart@net[["weight"]]
df <- as.data.frame(as.table(weight_matrix))
names(df) <- c("Cell1", 'Cell2',"Weight")
top3 <- df %>% arrange(desc(Weight)) %>% head(3)
print(top3)
top5 <- df %>% arrange(desc(Weight)) %>% head(5)
print(top5)
prob_matrix <- cellchat_heart@net[["prob"]]
df <- as.data.frame(as.table(prob_matrix))
names(df) <- c("Cell1", "Cell2",'LR','Prob')
pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/second/cellchat/heart_LR.pdf')
netVisual_chord_gene(cellchat_heart, sources.use = c('Endothelial_Do','Pericyte_Do','Myeloid_Do'), 
                     targets.use = c('Myeloid_Re'), lab.cex = 0.5,legend.pos.y = 20)
dev.off()

saveRDS(cellchat_heart,'/home/zeyy/dev/sdd/pyh/data/cellchat_ap.rds')


#artery tissue--------
tissue.data <- Read10X(data.dir = "/home/zeyy/dev/sdf/human/second_human/scRNA_cellranger/tissue/outs/filtered_feature_bc_matrix")
tissue<-CreateSeuratObject(counts = tissue.data, project = 'tissue',min.cells = 3, min.features = 200 )

tissue <- NormalizeData(tissue)
tissue <- FindVariableFeatures(tissue, selection.method = "vst", nfeatures = 2000)
tissue <- ScaleData(tissue, features = rownames(tissue))
tissue <- RunPCA(tissue, features = VariableFeatures(object = tissue))
ElbowPlot(tissue,ndims=30)

tissue <- FindNeighbors(tissue, dims = 1:20)
tissue <- FindClusters(tissue, resolution = 0.4)
tissue <- RunUMAP(tissue, dims = 1:20)
DimPlot(tissue, reduction = "umap",label=T)+NoAxes()

FeaturePlot(tissue,features = 'CCL3')+NoAxes()

tissue.markers <- FindAllMarkers(tissue, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tissue.markers<-tissue.markers[order(-tissue.markers$avg_log2FC),]
tissue.markers<-subset(tissue.markers,subset=p_val_adj<0.05)
write.csv(tissue.markers,'/home/zeyy/dev/sdd/pyh/result/second_heart/12.29/tissue_nul_marker_0.4.csv')

tissue[["percent.mt"]] <- PercentageFeatureSet(tissue, pattern = "^MT-")
ribo.genes <- grep("^RPS|^RPL", rownames(tissue), value = TRUE)
tissue[["percent.rb"]] <- PercentageFeatureSet(tissue, features = ribo.genes)
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- CaseMatch(HB.genes_total, rownames(tissue))
tissue[["percent.hb"]] <- PercentageFeatureSet(tissue, features=HB_m)
stress.genes_total <- read.csv('/home/kxm/pyh/data/stress.genes.csv')
stress.genes_total <- stress.genes_total$gene
stress_m <- CaseMatch(stress.genes_total, rownames(tissue))
tissue[["percent.stress"]] <- PercentageFeatureSet(tissue, features=stress_m)

new.cluster.ids <- c('Endothelial',#0
                     'Fibroblast',#1
                     'Cardiomyocyte',#2
                     'T',#3
                     'Myeloid',#4
                     'SMC',#5
                     'Endothelial',#6
                     'Adipocyte',#7
                     'Endothelial',#8
                     'Pericyte',#9
                     'B',#10
                     'Neuronal',#11
                     'Fibroblast',#12
                     'Mesothelial'#13
)
names(new.cluster.ids) <- levels(tissue)
tissue <- RenameIdents(tissue, new.cluster.ids)
tissue$celltype<-Idents(tissue)

##celltype,umap--------
lovely_color<-c("#00468BFF","#2673BF","#4DBBD5FF","#00A087FF","#74C476","#EBB365","#ED9A7C","#B15928","#8F0021",'#BCB9DA','yellow')
DimPlot(tissue,label=T)+NoAxes()
pdf("/home/zeyy/dev/sdd/pyh/result/trans_heart/second/artery_annotated.pdf",width=8.5)
DimPlot(tissue, reduction = "umap", label = T, pt.size = 0.5,label.size=5.2,cols=lovely_color)+
  NoAxes()+theme(legend.text = element_text(size = 18))
dev.off()

##add dmx----
demuxlet <- fread('/home/zeyy/dev/sdd/demuxlet/second/tissue.best')
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
tissue<-AddMetaData(tissue,dmx_new)

tissue$cell.orig = unlist(lapply(tissue$cell.orig, function(x) {
  if (is.na(x)) {return(NA)}
  else if (x == 'blood_2') {return('Recipient')} 
  else if (x == 'donor_2' | x == 'host_2') {return('Donor')} else {
    return(x)
  }
}))
##cell.orig,umap-----
DimPlot(tissue,group.by='cell.orig')+NoAxes()
tissue_SNG<-subset(tissue,subset=cell.quality=='SNG')
pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/second/artery_orig_umap.pdf',width=7.5)
DimPlot(tissue_SNG,group.by='cell.orig',cols = c("#ED855F","#548F9E"))+NoAxes()
dev.off()

saveRDS(tissue,'/home/zeyy/dev/sdd/pyh/data/second_tissue.rds')

##cell.orig,histogram------
tmp <- select(tissue_SNG@meta.data, c("celltype","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
celltype<-unique(tmp$celltype)
celltype<-celltype[c(1,4:6,8:10)]
celltype<-droplevels(celltype)
for(i in celltype){
  df_i <- subset(tmp, tmp$celltype==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("celltype", "cell.orig", "value")
  df <- rbind(df, df_i)
}
df$celltype <- factor(df$celltype, levels = c('SMC','Endothelial',"Fibroblast",
                                              "Pericyte",'Myeloid',"T","B"))
p <- ggplot(df, aes(x=celltype, y=value, fill=cell.orig)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/second/artery_orig.pdf',width = 6.5)
p
dev.off() 

##differential gene number,histogram--------
df<-data.frame()
for(i in celltype){
  seu_obj<-subset(tissue_SNG,ident=i)
  Idents(seu_obj)<-'cell.orig'
  CELLDEG <- FindMarkers(seu_obj, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
  CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
  CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
  RE <- sum(CELLDEG$avg_log2FC > 0)
  DO <- sum(CELLDEG$avg_log2FC < 0)
  df_RE<-data.frame(celltype=i,cell.orig='Recipient',value=RE)
  df_DO<-data.frame(celltype=i,cell.orig='Donor',value=DO)
  df_i<-rbind(df_RE,df_DO)
  df<- rbind(df, df_i)
}
celltype_order<-df%>%group_by(celltype)%>%
  summarise(value = sum(value))%>%arrange(desc(value))
df$celltype <- factor(df$celltype, levels = celltype_order$celltype)

p <- ggplot(df, aes(x=celltype, y=value, fill=cell.orig)) + 
  geom_col() + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  #scale_fill_npg()+
  labs(x = '', y = '') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/second/artery_redo_DEG.pdf',width = 6.5)
p
dev.off()

##endothelial differential analysis-------
#tissue_SNG<-subset(tissue,cell.quality=='SNG')
tissue_endo<-subset(tissue_SNG,celltype=='Endothelial')
Idents(tissue_endo)<-'cell.orig'
CELLDEG <- FindMarkers(tissue_endo, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
write.csv(CELLDEG,'/home/zeyy/dev/sdd/pyh/result/trans_heart/second/tissue_endo_DEG.csv')

##cellchat------
tissue_drop<-subset(tissue_SNG,celltype%in%c('Fibroblast','Endothelial','SMC','T',
                                             'Myeloid','B','Pericyte'))
metadata_info <- tissue_drop@meta.data[, c('celltype','cell.orig')]
#combine cell.orig&celltype
metadata_info$label <- paste(metadata_info$celltype, metadata_info$cell.orig, sep = "_")
data.input <- GetAssayData(tissue_drop, assay = "RNA", slot = "data")
data.input <- data.input[,rownames(metadata_info) ]

cellchat_artery <- createCellChat(object = data.input, meta = metadata_info, group.by = "label")

cellchat_artery@meta[["label"]] <- sub("_Recipient$", "_Re", cellchat_artery@meta[["label"]])
cellchat_artery@meta[["label"]] <- sub("_Donor$", "_Do", cellchat_artery@meta[["label"]])
cellchat_artery@idents <- as.factor(cellchat_artery@meta[['label']])

#preprocess
CellChatDB <- CellChatDB.human  
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat_artery@DB <- CellChatDB.use

cellchat_artery <- subsetData(cellchat_artery) 
future::plan("multisession", workers = 4) # do parallel
cellchat_artery <- identifyOverExpressedGenes(cellchat_artery)
cellchat_artery <- identifyOverExpressedInteractions(cellchat_artery)

#calculate probabilities
cellchat_artery <- computeCommunProb(cellchat_artery,population.size = TRUE)
cellchat_artery <- filterCommunication(cellchat_artery, min.cells = 10)
cellchat_artery <- computeCommunProbPathway(cellchat_artery)
cellchat_artery <- aggregateNet(cellchat_artery)

#visualization
groupSize <- as.numeric(table(cellchat_artery@idents))
pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/second/cellchat/artery_number.pdf')
netVisual_circle(cellchat_artery@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/second/cellchat/artery_weight.pdf')
netVisual_circle(cellchat_artery@net$weight, vertex.weight = groupSize, idents.use = c(),
                 weight.scale = T,  title.name = "Interaction weights/strength")
dev.off()
###LR interaction----
weight_matrix <- cellchat_artery@net[["weight"]]
df <- as.data.frame(as.table(weight_matrix))
names(df) <- c("Cell1", 'Cell2',"Weight")
top3 <- df %>% arrange(desc(Weight)) %>% head(3)
print(top3)
top5 <- df %>% arrange(desc(Weight)) %>% head(5)
print(top5)
prob_matrix <- cellchat_artery@net[["prob"]]
df <- as.data.frame(as.table(prob_matrix))
names(df) <- c("Cell1", "Cell2",'LR','Prob')
top10_LR_1 <- df %>%
  filter(Cell1 == "Endothelial_Re", Cell2 == "Endothelial_Do") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
top10_LR_2 <- df %>%
  filter(Cell1 == "T_Re", Cell2 == "Endothelial_Re") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
top10_LR_3 <- df %>%
  filter(Cell1 == "Endothelial_Re", Cell2 == "Endothelial_Re") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
combined_LR <- c(top10_LR_1, top10_LR_2, top10_LR_3)
combined_LR <- as.character(combined_LR)
unique_LR <- unique(combined_LR)#去重
unique_LR <- data.frame(interaction_name = unique_LR)
pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/second/cellchat/artery_LR.pdf')
netVisual_chord_gene(cellchat_artery, sources.use = c('Endothelial_Re','T_Re'), 
                     targets.use = c('Endothelial_Re','Endothelial_Do'), lab.cex = 0.5,legend.pos.y = 30,pairLR.use = unique_LR)
dev.off()

saveRDS(cellchat_artery,'/home/zeyy/dev/sdd/pyh/data/cellchat_tissue.rds')


