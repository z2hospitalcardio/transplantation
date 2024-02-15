#T/immune_endo correlation------
heart_1<-readRDS('/home/zeyy/dev/sdd/pyh/data/trans_heart_SNG.rds')
unique(Idents(heart_1))
heart_1<-subset(heart_1,celltype%in%c('Endothelial','Fibroblast','SMC','Myeloid','NK','T','Neutrophil'))
tmp <- select(heart_1@meta.data, c("orig.ident", "celltype","celltype_diplo","cell.orig"))
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  #T proportion
  # df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype) %>% table()
  # df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  # names(df_i) <- c("sample", "celltype", "value")
  # sum_i<-sum(df_i$value)
  # df_i<-subset(df_i,celltype=='T')
  # df_i$value<-df_i$value/sum_i
  
  #immune_infiltration
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype_diplo) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("sample", "celltype_diplo", "value")
  sum_i<-sum(df_i$value)
  df_i<-subset(df_i,celltype_diplo=='immune')
  df_i$value<-df_i$value/sum_i
  
  #endothelial
  df_j<-subset(tmp, tmp$orig.ident==i) %>% subset(celltype=='Endothelial') %>% pull(cell.orig) %>% table()
  df_j <- as.data.frame(df_j)
  names(df_j) <- c("cell.orig", "value")
  sum_j<-sum(df_j$value)
  df_j<-subset(df_j,cell.orig=='Recipient')
  df_i$endo<-df_j$value/sum_j
  df <- rbind(df, df_i)
}
df_heart1<-df


artery_1<-readRDS('/home/zeyy/dev/sdd/pyh/data/artery_SNG.rds')
unique(Idents(artery_1))
artery_1<-subset(artery_1,celltype%in%c('Endothelial','Fibroblast','SMC','Myeloid','NK','Tcell','B/Plasma','Myofibroblast'))
unique(Idents(artery_1))
new.celltype.ids <- c("stroma","stroma","stroma","immune","immune","immune","immune","stroma")
names(new.celltype.ids) <- levels(artery_1)
artery_1 <- RenameIdents(artery_1, new.celltype.ids)
artery_1$celltype_diplo <- Idents(artery_1)
tmp <- select(artery_1@meta.data, c("orig.ident", "celltype","celltype_diplo","cell.orig"))
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  #T比例
  # df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype) %>% table()
  # df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  # names(df_i) <- c("sample", "celltype", "value")
  # sum_i<-sum(df_i$value)
  # df_i<-subset(df_i,celltype=='Tcell')
  # df_i$value<-df_i$value/sum_i
  
  #immune infiltration
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype_diplo) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("sample", "celltype_diplo", "value")
  sum_i<-sum(df_i$value)
  df_i<-subset(df_i,celltype_diplo=='immune')
  df_i$value<-df_i$value/sum_i
  
  #endothelial
  df_j<-subset(tmp, tmp$orig.ident==i) %>% subset(celltype=='Endothelial') %>% pull(cell.orig) %>% table()
  df_j <- as.data.frame(df_j)
  names(df_j) <- c("cell.orig", "value")
  sum_j<-sum(df_j$value)
  df_j<-subset(df_j,cell.orig=='Recipient')
  df_i$endo<-df_j$value/sum_j
  df <- rbind(df, df_i)
}
df_artery1<-df

df1<-rbind(df_heart1,df_artery1)
df1<-df1[,c('sample','value','endo')]
names(df1)<-c('sample','Tcell','endo')

df2<-rbind(df_heart1,df_artery1)
df2<-df2[,c('sample','value','endo')]
names(df2)<-c('sample','immune','endo')

library(ggpubr)
pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/T_endo_corr.pdf',height = 4,width = 6)
ggplot(data=df1, aes(x=endo, y=Tcell))+geom_point(color="red")+stat_smooth(method="lm",se=T)+stat_cor(data=df1, method = "pearson")
dev.off()

pdf('/home/zeyy/dev/sdd/pyh/result/trans_heart/immune_endo_corr.pdf',height = 4,width = 6)
ggplot(data=df2, aes(x=endo, y=immune))+geom_point(color="red")+stat_smooth(method="lm",se=T)+stat_cor(data=df2, method = "pearson")
dev.off()


