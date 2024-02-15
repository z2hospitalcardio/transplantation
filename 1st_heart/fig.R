#cell type composition analysis-----
#bar_color<-c("#3B4992FF","#4DBBD5FF","#00A087FF","#91D1C2FF","#74C476","#F59899","#FCBF6E","#B15928","#AD002AFF")
bar_color<-c("#00468BFF","#2673BF","#4DBBD5FF","#00A087FF","#74C476","#EBB365","#ED9A7C","#B15928","#8F0021")
tmp <- select(artery@meta.data, c("orig.ident", "celltype"))
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("sample", "celltype", "value")
  df <- rbind(df, df_i)
}
df_artery<-df
df_artery$celltype <- as.character(df_artery$celltype)
df_artery$celltype[df_artery$celltype == 'Tcell'] <- 'T'
df_artery$celltype <- as.factor(df_artery$celltype)
df_heart$sample<-factor(df_heart$sample,levels=c('raf','ra',"la",'rv',"lvaw","lvpw","sp",'ap'))
df_artery$sample<-factor(df_artery$sample,levels=c('cx','ladb','rca'))
df<-rbind(df_heart,df_artery)
p <- ggplot(df, aes(x=sample, y=value, fill=celltype)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = bar_color) +
  #scale_fill_npg()+
  labs(x = '', y = '') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf("/home/zeyy/dev/sdd/pyh/result/proportion_all.pdf",width=8)
p
dev.off()

#immune and stroma cell composition--------
new.celltype.ids <- c("stroma","stroma","stroma","immune","immune","immune",
                      "immune","stroma")
names(new.celltype.ids) <- levels(artery)
artery <- RenameIdents(artery, new.celltype.ids)
artery$celltype_diplo <- Idents(artery)

tmp <- select(artery@meta.data, c("orig.ident", "celltype_diplo"))
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype_diplo) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("sample", "celltype_diplo", "value")
  df <- rbind(df, df_i)
}
df_artery_diplo<-df
df_heart_diplo$sample<-factor(df_heart_diplo$sample,levels=c('raf','ra',"la",'rv',"lvaw","lvpw","sp",'ap'))
df_artery_diplo$sample<-factor(df_artery_diplo$sample,levels=c('cx','ladb','rca'))
df<-rbind(df_heart_diplo,df_artery_diplo)
new_colors<-c("#00468BFF","#8F0021")
p <- ggplot(df, aes(x=sample, y=value, fill=celltype_diplo)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = new_colors) +
  #scale_fill_npg()+
  labs(x = '', y = '') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf("/home/zeyy/dev/sdd/pyh/result/proportion_diplo.pdf",width = 7.5)
p
dev.off()


#recipient/donor histogram in all orig.ident--------
#heart
heart<-readRDS('/home/zeyy/dev/sdd/pyh/data/trans_heart_dmx.rds')
heart_drop<-subset(heart,subset=cell.quality=='SNG')
heart_drop$cell.orig[heart_drop$cell.orig == 'Re'] <- 'Recipient'
heart_drop$cell.orig[heart_drop$cell.orig == 'AP'] <- 'Donor'
heart_drop$orig.ident<-factor(heart_drop$orig.ident,levels=c('raf','ra','la','rv','lvaw','lvpw','ap','sp'))
tmp <- select(heart_drop@meta.data, c("orig.ident","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("orig.ident", "cell.orig", "value")
  df <- rbind(df, df_i)
}
heart_df<-df

#artery
artery_drop<-subset(artery,subset=cell.quality=='SNG')
tmp <- select(artery_drop@meta.data, c("orig.ident","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("orig.ident", "cell.orig", "value")
  df <- rbind(df, df_i)
}
artery_df<-df

#pbmc
pbmc<-readRDS('/home/zeyy/dev/sdd/gly/blood/pbmc_anno.rds')
pbmc_dmx<-readRDS('/home/zeyy/dev/sdd/pyh/data/pbmc_dxm.rds')
metadata_info <- pbmc_dmx@meta.data[, c('cell.quality','cell.orig')]
pbmc <- AddMetaData(object = pbmc, metadata = metadata_info)
saveRDS(pbmc,'/home/zeyy/dev/sdd/pyh/data/pbmc_dmx.rds')
pbmc$cell.orig[pbmc$cell.orig == 'Re'] <- 'Recipient'
pbmc$cell.orig[pbmc$cell.orig == 'AP'] <- 'Donor'
pbmc_drop<-subset(pbmc,subset=cell.quality=='SNG')

tmp <- select(pbmc_drop@meta.data, c("cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
df <- tmp %>% pull(cell.orig) %>% table()
df <- data.frame(sample=rep('blood', length(df)),value=df)
names(df) <- c("orig.ident","cell.orig", "value")
pbmc_df<-df

df<-rbind(heart_df,artery_df,pbmc_df)
df$orig.ident <- factor(df$orig.ident, levels = c('raf','ra','la','rv','lvaw','lvpw','sp','ap',
                                                'cx','ladb','rca','blood'))

p <- ggplot(df, aes(x=orig.ident, y=value, fill=cell.orig)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/figE.pdf')
p
dev.off()  


#raf，heart，artery，pbmc,Recipient&Donor_umap-------
raf<-subset(heart_drop,subset=orig.ident=='raf')
pdf('/home/zeyy/dev/sdd/pyh/result/figD/raf.pdf',width=8)
DimPlot(raf,group.by='cell.orig',cols = c("#ED855F","#548F9E"))+NoAxes()
dev.off()

pdf('/home/zeyy/dev/sdd/pyh/result/figD/heart.pdf',width=8)
DimPlot(heart_drop,group.by='cell.orig',cols = c("#ED855F","#548F9E"))+NoAxes()
dev.off()

pdf('/home/zeyy/dev/sdd/pyh/result/figD/artery.pdf',width=8)
DimPlot(artery_drop,group.by='cell.orig',cols = c("#ED855F","#548F9E"))+NoAxes()
dev.off()

pdf('/home/zeyy/dev/sdd/pyh/result/figD/pbmc.pdf',width=8)
DimPlot(pbmc_drop,group.by='cell.orig',cols = c("#ED855F","#548F9E"))+NoAxes()
dev.off()


#atrium,ventricle,artery,Recipient&Donor_umap-------
#atrium
pos<-c('ra','la')
atrium<-subset(heart_drop,subset=orig.ident %in% pos)
pdf('/home/zeyy/dev/sdd/pyh/result/figF/atrium_umap.pdf',width=8)
p<-DimPlot(atrium,group.by='cell.orig',cols = c("#ED855F","#548F9E"))+NoAxes()
p
dev.off()

#ventricle
pos<-c('rv','lvaw','lvpw')
ventricle<-subset(heart_drop,subset=orig.ident %in% pos)
pdf('/home/zeyy/dev/sdd/pyh/result/figF/ventricle_umap.pdf',width=8)
p<-DimPlot(ventricle,group.by='cell.orig',cols = c("#ED855F","#548F9E"))+NoAxes()
p
dev.off()

#artery
pos<-artery_drop$orig.ident
for(i in unique(pos)){
  object<-subset(artery_drop,subset=orig.ident==i)
  pdf(paste0('/home/zeyy/dev/sdd/pyh/result/figF/',i,'_umap.pdf'),width=8)
  p<-DimPlot(object,group.by='cell.orig',cols = c("#ED855F","#548F9E"))+NoAxes()
  p
  print(p)
  dev.off()
}

#atrium,ventricle,artery,Recipient&Donor_celltype_histogram-------
#ventricle
pos<-c('rv','lvaw','lvpw')
heart_v<-subset(heart_drop,subset=orig.ident %in% pos)
tmp <- select(ventricle@meta.data, c("celltype","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
for(i in unique(tmp$celltype)){
  df_i <- subset(tmp, tmp$celltype==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("celltype", "cell.orig", "value")
  df <- rbind(df, df_i)
}
df$celltype <- factor(df$celltype, levels = c('SMC','Endothelial',"Fibroblast",
                                              "Neutrophil",'Myeloid',"NK","T"))
p <- ggplot(df, aes(x=celltype, y=value, fill=cell.orig)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/figG/ventricle.pdf',width = 6.5)
p
dev.off()  

#atrium
pos<-c('ra','la')
heart_a<-subset(heart_drop,subset=orig.ident %in% pos)
tmp <- select(atrium@meta.data, c("celltype","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
for(i in unique(tmp$celltype)){
  df_i <- subset(tmp, tmp$celltype==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("celltype", "cell.orig", "value")
  df <- rbind(df, df_i)
}
df$celltype <- factor(df$celltype, levels = c('SMC','Endothelial',"Fibroblast",
                                              "Neutrophil",'Myeloid',"NK","T"))
p <- ggplot(df, aes(x=celltype, y=value, fill=cell.orig)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/figG/atrium.pdf',width=6.5)
p
dev.off()  

#artery
tmp <- select(artery_SNG@meta.data, c("celltype","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
for(i in unique(tmp$celltype)){
  df_i <- subset(tmp, tmp$celltype==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("celltype", "cell.orig", "value")
  df <- rbind(df, df_i)
}
df$celltype <- factor(df$celltype, levels = c('SMC','Endothelial',"Fibroblast",'Myofibroblast',
                                              'Myeloid',"NK","Tcell","B/Plasma"))
p <- ggplot(df, aes(x=celltype, y=value, fill=cell.orig)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/figG/artery.pdf',width = 6.5)
p
dev.off() 

#pbmc
tmp <- select(pbmc_SNG@meta.data, c("celltype","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
for(i in unique(tmp$celltype)){
  df_i <- subset(tmp, tmp$celltype==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("celltype", "cell.orig", "value")
  df <- rbind(df, df_i)
}
#df$celltype <- factor(df$celltype, levels = c('Cardiomyocyte','SMC','Endothelial',"Fibroblast",'Myofibroblast','Neuronal',
#                                              'Myeloid',"NK","Tcell","B/Plasma"))
p <- ggplot(df, aes(x=celltype, y=value, fill=cell.orig)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/figG/pbmc.pdf')
p
dev.off() 

#endothelial in atrium,ventricle,artery,Recipient&Donor_subtype_histogram-------
heart_endo<-readRDS('/home/zeyy/dev/sdd/gly/Endo/heart_Endo.rds')
#ventricle
pos<-c('rv','lvaw','lvpw')
heart_endo_v<-subset(heart_endo,subset=orig.ident %in% pos)
heart_endo_v<-subset(heart_endo_v,subset=cell.quality=='SNG')
heart_endo_v$cell.orig[heart_endo_v$cell.orig == 'Re'] <- 'Recipient'
heart_endo_v$cell.orig[heart_endo_v$cell.orig == 'AP'] <- 'Donor'
tmp <- select(heart_endo_v@meta.data, c("celltype_subset","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
for(i in unique(tmp$celltype_subset)){
  df_i <- subset(tmp, tmp$celltype_subset==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("celltype_subset", "cell.orig", "value")
  df <- rbind(df, df_i)
}
#df$celltype_subset <- factor(df$celltype_subset, levels = c('Cardiomyocytes','SMC','Endothelial',"Fibroblast",'Neuronal',
#                                              "Neutrophil",'Myeloid',"NK","T"))
p <- ggplot(df, aes(x=celltype_subset, y=value, fill=cell.orig)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/figH/ventricle.pdf')
p
dev.off()  

#atrium
pos<-c('ra','la')
heart_endo_a<-subset(heart_endo,subset=orig.ident %in% pos)
heart_endo_a<-subset(heart_endo_a,subset=cell.quality=='SNG')
heart_endo_a$cell.orig[heart_endo_a$cell.orig == 'Re'] <- 'Recipient'
heart_endo_a$cell.orig[heart_endo_a$cell.orig == 'AP'] <- 'Donor'
tmp <- select(heart_endo_a@meta.data, c("celltype_subset","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
for(i in unique(tmp$celltype_subset)){
  df_i <- subset(tmp, tmp$celltype_subset==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("celltype_subset", "cell.orig", "value")
  df <- rbind(df, df_i)
}
#df$celltype_subset <- factor(df$celltype_subset, levels = c('Cardiomyocytes','SMC','Endothelial',"Fibroblast",'Neuronal',
#                                              "Neutrophil",'Myeloid',"NK","T"))
p <- ggplot(df, aes(x=celltype_subset, y=value, fill=cell.orig)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/figH/atrium.pdf')
p
dev.off()  

#artery
artery_endo<-readRDS('/home/zeyy/dev/sdd/gly/artery/endo/rpca/art_endo.rds')
tmp <- select(artery_endo@meta.data, c("celltype_subset","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
for(i in unique(tmp$celltype_subset)){
  df_i <- subset(tmp, tmp$celltype_subset==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("celltype_subset", "cell.orig", "value")
  df <- rbind(df, df_i)
}
p <- ggplot(df, aes(x=celltype_subset, y=value, fill=cell.orig)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/figH/artery.pdf')
p
dev.off()  


#myeloid in atrium,ventricle,artery,Recipient&Donor_subtype_histogram---------
#artery
artery_dmx<-readRDS('/home/zeyy/dev/sdd/pyh/data/artery_SNG.rds')
artery_myeloid$cell.orig<-artery_dmx$cell.orig
artery_myeloid$cell.quality<-artery_dmx$cell.quality
artery_myeloid_drop<-subset(artery_myeloid,subset=cell.quality=='SNG')
tmp <- select(artery_myeloid_drop@meta.data, c("celltype_sub","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
for(i in unique(tmp$celltype_sub)){
  df_i <- subset(tmp, tmp$celltype_sub==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("celltype_sub", "cell.orig", "value")
  df <- rbind(df, df_i)
}
p <- ggplot(df, aes(x=celltype_sub, y=value, fill=cell.orig)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/myeloid/artery.pdf')
p
dev.off()  

heart_myeloid_test<-heart_myeloid
heart_myeloid_test[["umap_naive"]] <- NULL

#ventricle
tmp <- select(ventricle_myeloid@meta.data, c("celltype_subset","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
for(i in unique(tmp$celltype_subset)){
  df_i <- subset(tmp, tmp$celltype_subset==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("celltype_subset", "cell.orig", "value")
  df <- rbind(df, df_i)
}
p <- ggplot(df, aes(x=celltype_subset, y=value, fill=cell.orig)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/myeloid/ventricle.pdf')
p
dev.off()  

#atrium
tmp <- select(atrium_myeloid@meta.data, c("celltype_subset","cell.orig"))
tmp<-na.omit(tmp)
df <- data.frame()
for(i in unique(tmp$celltype_subset)){
  df_i <- subset(tmp, tmp$celltype_subset==i) %>% pull(cell.orig) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("celltype_subset", "cell.orig", "value")
  df <- rbind(df, df_i)
}
df$cell.orig<- factor(df$cell.orig, levels=c('Donor','Recipient'))
p <- ggplot(df, aes(x=celltype_subset, y=value, fill=cell.orig)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/myeloid/atrium.pdf')
p
dev.off()  


#endothelial differential analysis--------
#recipient/donor
heart_endo$cell.orig<-heart_drop$cell.orig
RE_DO<-subset(heart_endo,subset=cell.orig %in% c("Recipient","Donor"))
Idents(RE_DO) <- "cell.orig"
CELLDEG <- FindMarkers(RE_DO, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
write.csv(CELLDEG, file = "/home/zeyy/dev/sdd/pyh/result/figI/RE_DO.csv")
#atrium
heart_endo<-subset(heart_SNG,ident='Endothelial')
pos<-c('ra','la')
atrium_endo<-subset(heart_endo,subset=orig.ident%in%pos)
Idents(atrium_endo)<-'cell.orig'
CELLDEG <- FindMarkers(atrium_endo, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
saveRDS(CELLDEG,'/home/zeyy/dev/sdd/pyh/result/figI/atrium_RE_DO.rds')
write.csv(CELLDEG, file = "/home/zeyy/dev/sdd/pyh/result/figI/atrium_RE_DO.csv")
#entricle
pos<-c('rv','lvaw','lvpw')
ventricle_endo<-subset(heart_endo,subset=orig.ident%in%pos)
Idents(ventricle_endo)<-'cell.orig'
CELLDEG <- FindMarkers(ventricle_endo, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
saveRDS(CELLDEG,'/home/zeyy/dev/sdd/pyh/result/figI/ventricle_RE_DO.rds')
write.csv(CELLDEG, file = "/home/zeyy/dev/sdd/pyh/result/figI/ventricle_RE_DO.csv")
#artery
artery_endo<-subset(artery_SNG,ident='Endothelial')
Idents(artery_endo)<-'cell.orig'
CELLDEG <- FindMarkers(artery_endo, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
saveRDS(CELLDEG,'/home/zeyy/dev/sdd/pyh/result/figI/artery_RE_DO.rds')
write.csv(CELLDEG, file = "/home/zeyy/dev/sdd/pyh/result/figI/artery_RE_DO.csv")

#raf/orig.ident
Idents(atrium)<-'orig.ident'
cellfordeg<-levels(atrium$orig.ident)
cellfordeg <- setdiff(cellfordeg, "raf")
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(atrium, ident.1 = cellfordeg[i], ident.2 ='raf', verbose = FALSE,min.pct = 0.1)
  CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
  CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
  write.csv(CELLDEG, file = paste0("/home/zeyy/dev/sdd/pyh/result/figI/", cellfordeg[i], "_raf.csv"))
}


#myeloid differential analysis-------
#recipient/donor
heart_Myeloid<-subset(heart_drop,ident='Myeloid')
Idents(heart_Myeloid) <- "cell.orig"
CELLDEG <- FindMarkers(heart_Myeloid, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
write.csv(CELLDEG, file = "/home/zeyy/dev/sdd/pyh/result/figK/RE_DO.csv")
#atrium
CELLDEG <- FindMarkers(atrium_myeloid, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
saveRDS(CELLDEG,'/home/zeyy/dev/sdd/pyh/result/figK/atrium_RE_DO.rds')
#ventricle
CELLDEG <- FindMarkers(ventricle_myeloid, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
saveRDS(CELLDEG,'/home/zeyy/dev/sdd/pyh/result/figK/ventricle_RE_DO.rds')
#artery
CELLDEG <- FindMarkers(artery_myeloid_drop, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
saveRDS(CELLDEG,'/home/zeyy/dev/sdd/pyh/result/figK/artery_RE_DO.rds')
write.csv(CELLDEG, file = "/home/zeyy/dev/sdd/pyh/result/figK/artery_RE_DO.csv")

#raf/orig.ident
Idents(heart_Myeloid)<-'orig.ident'
cellfordeg<-levels(heart_Myeloid$orig.ident)
cellfordeg <- setdiff(cellfordeg, "raf")
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(heart_Myeloid, ident.1 = cellfordeg[i], ident.2 ='raf', verbose = FALSE,min.pct = 0.1)
  CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]#对数据进行排序
  CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
  write.csv(CELLDEG, file = paste0("/home/zeyy/dev/sdd/pyh/result/figK/", cellfordeg[i], "_raf.csv"))
}


#recipient vs. donor,differential gene number histogram-------
#atrium
pos<-c('la','ra')
atrium<-subset(heart_SNG,subset=orig.ident %in% pos)

cellfordeg<-unique(atrium$celltype)
df_atrium<-data.frame()
for(i in cellfordeg){
  seu_obj<-subset(atrium,ident=i)
  Idents(seu_obj)<-'cell.orig'
  CELLDEG <- FindMarkers(seu_obj, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
  CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
  CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
  RE <- sum(CELLDEG$avg_log2FC > 0)
  DO <- sum(CELLDEG$avg_log2FC < 0)
  df_RE<-data.frame(celltype=i,cell.orig='Recipient',value=RE)
  df_DO<-data.frame(celltype=i,cell.orig='Donor',value=DO)
  df_i<-rbind(df_RE,df_DO)
  df_atrium <- rbind(df_atrium, df_i)
}
df_atrium$celltype <- factor(df_atrium$celltype, levels = c('NK','Endothelial','Fibroblast',
                                                            'SMC','Myeloid','T',
                                                            'Neutrophil'))

p <- ggplot(df_atrium, aes(x=celltype, y=value, fill=cell.orig)) + 
  geom_col() + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  #scale_fill_npg()+
  labs(x = '', y = '') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/number_deg_REDO/atrium.pdf',width = 6.5)
p
dev.off()

#ventricle
pos<-c('rv','lvaw','lvpw')
ventricle<-subset(heart_SNG,subset=orig.ident %in% pos)

cellfordeg<-unique(ventricle$celltype)
df_ventricle<-data.frame()
for(i in cellfordeg){
  seu_obj<-subset(ventricle,ident=i)
  Idents(seu_obj)<-'cell.orig'
  CELLDEG <- FindMarkers(seu_obj, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
  CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
  CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
  RE <- sum(CELLDEG$avg_log2FC > 0)
  DO <- sum(CELLDEG$avg_log2FC < 0)
  df_RE<-data.frame(celltype=i,cell.orig='Recipient',value=RE)
  df_DO<-data.frame(celltype=i,cell.orig='Donor',value=DO)
  df_i<-rbind(df_RE,df_DO)
  df_ventricle <- rbind(df_ventricle, df_i)
}

df_ventricle$celltype <- factor(df_ventricle$celltype, levels = c('Endothelial','SMC',
                                                                  'NK','Fibroblast','Myeloid',
                                                                  'Neutrophil','T'))

p <- ggplot(df_ventricle, aes(x=celltype, y=value, fill=cell.orig)) + 
  geom_col() + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  #scale_fill_npg()+
  labs(x = '', y = '') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/number_deg_REDO/ventricle.pdf',width = 6.5)
p
dev.off()

#artery
cellfordeg<-unique(artery_SNG$celltype)
df_artery_SNG<-data.frame()
for(i in cellfordeg){
  seu_obj<-subset(artery_SNG,ident=i)
  Idents(seu_obj)<-'cell.orig'
  CELLDEG <- FindMarkers(seu_obj, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
  CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
  CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
  RE<- sum(CELLDEG$avg_log2FC > 0)
  DO <- sum(CELLDEG$avg_log2FC < 0)
  df_RE<-data.frame(celltype=i,cell.orig='Recipient',value=RE)
  df_DO<-data.frame(celltype=i,cell.orig='Donor',value=DO)
  df_i<-rbind(df_RE,df_DO)
  df_artery_SNG <- rbind(df_artery_SNG, df_i)
}
df_artery_SNG$celltype <- factor(df_artery_SNG$celltype, levels = c('Fibroblast','Endothelial','Myeloid',
                                                                    'SMC','Myofibroblast','Tcell',
                                                                    'NK','B/Plasma'))

p <- ggplot(df_artery_SNG, aes(x=celltype, y=value, fill=cell.orig)) + 
  geom_col() + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  #scale_fill_npg()+
  labs(x = '', y = '') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/number_deg_REDO/artery.pdf')
p
dev.off()

#pbmc
cellfordeg<-unique(pbmc_SNG$celltype)
df_pbmc_SNG<-data.frame()
for(i in cellfordeg){
  seu_obj<-subset(pbmc_SNG,ident=i)
  Idents(seu_obj)<-'cell.orig'
  CELLDEG <- FindMarkers(seu_obj, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
  CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
  CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
  RE<- sum(CELLDEG$avg_log2FC > 0)
  DO <- sum(CELLDEG$avg_log2FC < 0)
  df_RE<-data.frame(celltype=i,cell.orig='Recipient',value=RE)
  df_DO<-data.frame(celltype=i,cell.orig='Donor',value=DO)
  df_i<-rbind(df_RE,df_DO)
  df_pbmc_SNG <- rbind(df_pbmc_SNG, df_i)
}
df_pbmc_SNG$celltype <- factor(df_pbmc_SNG$celltype, levels = c('DC','CD4+_T','CD8+_T',
                                                                'NK','NKT','Mono/Macro',
                                                                'B','Th17'))

p <- ggplot(df_pbmc_SNG, aes(x=celltype, y=value, fill=cell.orig)) + 
  geom_col() + 
  scale_fill_manual(values = c("#ED855F","#548F9E")) +
  #scale_fill_npg()+
  labs(x = '', y = '') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45,hjust = 1,size=20),axis.text.y = element_text(size=20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size=22),
        axis.ticks.x = element_blank())
pdf('/home/zeyy/dev/sdd/pyh/result/number_deg_REDO/pbmc.pdf')
p
dev.off()

cellfordeg<-unique(pbmc_SNG$celltype)
df_pbmc_SNG<-data.frame()
for(i in 5:length(cellfordeg)){
  seu_obj<-subset(pbmc_SNG,ident=cellfordeg[i])
  Idents(seu_obj)<-'cell.orig'
  CELLDEG <- FindMarkers(seu_obj, ident.1 = 'Recipient', ident.2 = 'Donor', verbose = FALSE,min.pct = 0.1)
  CELLDEG <- CELLDEG[order(CELLDEG$avg_log2FC, decreasing = TRUE),]
  CELLDEG <- subset(CELLDEG,CELLDEG$p_val_adj<0.05)
  write.csv(CELLDEG, file = paste0("/home/zeyy/dev/sdd/pyh/result/pbmc_DEG/", cellfordeg[i], ".CSV"))
}


#each orig.ident,recipient&donor umap------
artery_SNG<-subset(artery,subset=cell.quality=='SNG')
cx<-subset(artery_SNG,subset=orig.ident=='cx')
rca<-subset(artery_SNG,subset=orig.ident=='rca')
ladb<-subset(artery_SNG,subset=orig.ident=='ladb')
p1<-DimPlot(cx,group.by='cell.orig',cols = c('#E11E24','#92D0E6'))+NoAxes()+ggtitle('cx')
p2<-DimPlot(rca,group.by='cell.orig',cols = c('#E11E24','#92D0E6'))+NoAxes()+ggtitle('rca')
p3<-DimPlot(ladb,group.by='cell.orig',cols = c('#E11E24','#92D0E6'))+NoAxes()+ggtitle('ladb')
p<-p1+p2+p3
pdf('/home/zeyy/dev/sdd/pyh/result/supp/artery_REDO.pdf',width=15,height=5)
p
dev.off()

#cell marker------
artery.markers <- FindAllMarkers(artery, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(artery.markers,'/home/kxm/pyh/result/ClusterMarker_artery.csv')

pbmc<-subset(pbmc,subset=celltype!='Erythrocyte')
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers,'/home/kxm/pyh/result/ClusterMarker_pbmc.csv')


library(RColorBrewer)
my_palette <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100))
p1<-DotPlot(artery,features = c('VWF','PECAM1','ACKR1','EGFL7','CLDN5',#Endo
                                'DCN','CFD','GSN','LUM','FBLN1',#Fibro
                                'MYH11','TAGLN','ACTA2','MYL9','MUSTN1',#SMC
                                'CD163','C1QA','LYZ','FCER1G','CD68',#Myeloid
                                'GNLY','NKG7','GZMB','FGFBP2','KLRD1',#NK
                                'IL7R','CD3D','CD69','RUNX3','CD3E',#T
                                'MYL3','MYH7','MYL7','CSRP3','ACTN2',#Cardio
                                'IGHG1','CD79A','MZB1','MS4A1','TNFRSF13C',#B
                                'NRXN1','S100B','PLP1','GPM6B','LGI4',#Neuronal
                                'COL1A1','COL1A2','POSTN'#Myofib
))+
  #coord_flip()+
  RotatedAxis()+
  scale_color_gradientn(colors = my_palette)+ #theme_bw()+
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black"), panel.grid = element_line(colour = "grey92"), 
        strip.background = element_rect(fill = "grey85", 
                                        colour = "grey20"),panel.background = element_rect(colour = 'black'), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        legend.key.size = unit(0.5, "cm"), legend.position = "right"
  )+
  RotatedAxis()+
  ggtitle('artery')

p2<-DotPlot(pbmc,features = c('CD3D','IL7R','CD40LG','LEF1','MAL',#CD4T
                              'CD8A','CD8B','CCL5','GZMH',#CD8T
                              'LYZ','CD68','CD14','VCAN','S100A8',#Mono
                              'CD74','CD79A','MS4A1','HLA-DRA','IGHM',#B
                              'ZNF683','TRBV6-2','FGFBP2',#NKT
                              'FCGR3A','FCGR3B','GNLY','NKG7','NCAM1',#NK
                              'KLRB1','RORC','IL23R',#Th17
                              'FCER1A','CLEC10A','PLD4'#DC
))+
  #coord_flip()+
  RotatedAxis()+
  scale_color_gradientn(colors = my_palette)+ #theme_bw()+
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black"), panel.grid = element_line(colour = "grey92"), 
        strip.background = element_rect(fill = "grey85", 
                                        colour = "grey20"),panel.background = element_rect(colour = 'black'), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        legend.key.size = unit(0.5, "cm"), legend.position = "right"
  )+
  RotatedAxis()+
  ggtitle('pbmc')

p<-p1+p2+plot_layout(ncol = 1, byrow = FALSE)

pdf('/home/zeyy/dev/sdd/pyh/result/supp/marker.pdf',width=15,height=12)
p
dev.off()




