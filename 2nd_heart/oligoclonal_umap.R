#oligoclonal-------
#second_heart
second_heart<-readRDS('/home/zeyy/dev/sdd/pyh/data/second_heart.rds')
second_artery<-readRDS('/home/zeyy/dev/sdd/pyh/data/second_artery.rds')
second_pbmc<-readRDS('/home/zeyy/dev/sdd/pyh/data/second_pbmc.rds')

#heart add vdj----
heartT<-subset(second_heart,ident='T')
S1<-read.csv("/home/zeyy/dev/sdd/second_human/TCR_BCR/TCR/LVB/filtered_contig_annotations.csv",stringsAsFactors = F)
S2<-read.csv("/home/zeyy/dev/sdd/second_human/TCR_BCR/TCR/LVPW/filtered_contig_annotations.csv",stringsAsFactors = F)
contig_list <- list(S1,S2)
combined <- combineTCR(contig_list, 
                       samples =c('lvb','lvpw'))
heartT <- combineExpression(combined, heartT, 
                               cloneCall="gene", 
                               proportion = FALSE)

#artery add vdj----
arteryT<-subset(second_artery,ident='T')
S1<-read.csv("/home/zeyy/dev/sdd/second_human/TCR_BCR/TCR/CX/filtered_contig_annotations.csv",stringsAsFactors = F)
S2<-read.csv("/home/zeyy/dev/sdd/second_human/TCR_BCR/TCR/LAD/filtered_contig_annotations.csv",stringsAsFactors = F)
contig_list <- list(S1,S2)
combined <- combineTCR(contig_list, 
                       samples =c('cx','lad'), 
                       ID = c('artery','artery'))
combined$cx_artery$barcode=gsub('cx_artery','cx',combined$cx_artery$barcode)
combined$lad_artery$barcode=gsub('lad_artery','lad',combined$lad_artery$barcode)
arteryT <- combineExpression(combined, arteryT, 
                            cloneCall="gene", 
                            proportion = FALSE)

#pbmc add VDJ------
pbmcT<-subset(second_pbmc,ident='T')
contig<-read.csv("/home/zeyy/dev/sdd/second_human/TCR_BCR/TCR/PBMC/filtered_contig_annotations.csv",stringsAsFactors = F)
contig_list <- list(contig)
combined <- combineTCR(contig_list, 
                       samples ='contig')
pbmcT <- combineExpression(combined, pbmcT, 
                            cloneCall="gene", 
                            proportion = FALSE)


Hartery<-subset(arteryT,subset=CTaa%in%c('CAVESPPDDKIIF_CASSLTIAGGTGELFF',
                                         'NA_CASSLTIAGGTGELFF'))
Hheart<-subset(heartT,subset=CTaa%in%c('CAVESPPDDKIIF_CASSLTIAGGTGELFF',
                                       'NA_CASSLTIAGGTGELFF',
                                       'NA_CASSQRRGATEAFF',
                                       'CALVYYNTDKLIF_CASSQRRGATEAFF'))
Hpbmc<-subset(pbmcT,subset=CTaa%in%c('CAVESPPDDKIIF_CASSLTIAGGTGELFF',
                                     'NA_CASSLTIAGGTGELFF'))


# add CTaa metadata to second_heart------
metadatainfo <- heartT$CTaa
new_metadata <- rep(NA, ncol(second_heart))
names(new_metadata) <- colnames(second_heart)
new_metadata[names(metadatainfo)] <- metadatainfo
second_heart <- AddMetaData(second_heart, metadata = new_metadata, col.name = 'CTaa')
pdf('/home/zeyy/dev/sdd/pyh/result/second_heart/hyper.pdf',width = 10)
DimPlot(second_heart,group.by='CTaa')
dev.off()
saveRDS(second_heart,'/home/zeyy/dev/sdd/pyh/data/second_heart.rds')


# add CTaa metadata to second_pbmc-----
second_pbmc<-readRDS('/home/zeyy/dev/sdd/pyh/data/second_pbmc.rds')
metadatainfo <- Hpbmc$CTaa
new_metadata <- rep(NA, ncol(second_pbmc))
names(new_metadata) <- colnames(second_pbmc)
new_metadata[names(metadatainfo)] <- metadatainfo
second_pbmc <- AddMetaData(second_pbmc, metadata = new_metadata, col.name = 'CTaa')
pdf('/home/zeyy/dev/sdd/pyh/result/second_pbmc/hyper.pdf',width = 10)
DimPlot(second_pbmc,group.by='CTaa')
dev.off()
saveRDS(second_pbmc,'/home/zeyy/dev/sdd/pyh/data/second_pbmc.rds')


# add CTaa metadata to second_artery------
second_artery<-readRDS('/home/zeyy/dev/sdd/pyh/data/second_artery.rds')
metadatainfo <- Hartery$CTaa
new_metadata <- rep(NA, ncol(second_artery))
names(new_metadata) <- colnames(second_artery)
new_metadata[names(metadatainfo)] <- metadatainfo
second_artery <- AddMetaData(second_artery, metadata = new_metadata, col.name = 'CTaa')
pdf('/home/zeyy/dev/sdd/pyh/result/second_artery/hyper.pdf',width = 10)
DimPlot(second_artery,group.by='CTaa')
dev.off()
saveRDS(second_artery,'/home/zeyy/dev/sdd/pyh/data/second_artery.rds')

