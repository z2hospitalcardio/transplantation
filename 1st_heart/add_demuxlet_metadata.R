library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)

setwd("/home/zeyy/dev/sdd/demuxlet")

#heart demuxlet---------
#Get all files ending in '.best' in the current folder
files <- list.files(pattern = "\\.best$")
files<-files[1:8]
trans_heart<-readRDS('/home/zeyy/dev/sdd/pyh/data/trans_heart_annotationed.rds')
scelist<-SplitObject(trans_heart,split.by='orig.ident')
#change
names(scelist)[names(scelist) == "lvpw"] <- "lvib"
names(scelist)[names(scelist) == "lvaw"] <- "lvna"

scelist_n<-lapply(files, function(file) {
  #数Get the part of the file name before '.best'
  pos <- sub("\\.best$", "", file)
  demuxlet <- fread(file)
  #change barcode
  demuxlet$BARCODE<-paste0(pos,'_',demuxlet$BARCODE)
  #drop empty cell
  demuxlet<-na.omit(demuxlet)
  dmx_new <- data.frame(barcode=demuxlet$BARCODE,cell.quality=sapply(demuxlet$BEST,function(x){strsplit(x,"-")[[1]][[1]]}),
                    cell.orig=sapply(demuxlet$BEST,function(x){strsplit(x,"-")[[1]][[2]]}))
  row.names(dmx_new) <- dmx_new$barcode
  dmx_new$barcode <- NULL
  
  i <- which(names(scelist) == pos)
  object<-AddMetaData(scelist[[i]],dmx_new)
  return(object)
})
names(scelist_n)<-c('ap','la','lvpw','lvaw','ra','raf','rv','sp')

#change back
# names(scelist)[names(scelist) == "lvib"] <- "lvpw"
# names(scelist)[names(scelist) == "lvna"] <- "lvaw"
source("/home/kxm/code/function.R")
orig.ident<-names(scelist_n)

trans_heart_dmx <- merge(scelist_n[[1]], y = scelist_n[2:8])

##extract metadata from dmx and add to trans_heart---------
metadata_info <- trans_heart_dmx@meta.data[, c('cell.quality','cell.orig')]
trans_heart <- AddMetaData(object = trans_heart, metadata = metadata_info)


setwd("/home/zeyy/dev/sdd/demuxlet/artery/")

#artery_demuxlet---------
files <- list.files(pattern = "\\.best$")
scelist<-SplitObject(artery,split.by='orig.ident')

scelist_n<-lapply(files, function(file) {
  pos <- sub("\\.best$", "", file)
  demuxlet <- fread(file)
  demuxlet$BARCODE<-paste0(pos,'_',demuxlet$BARCODE)
  demuxlet<-na.omit(demuxlet)
  dmx_new <- data.frame(barcode=demuxlet$BARCODE,cell.quality=sapply(demuxlet$BEST,function(x){strsplit(x,"-")[[1]][[1]]}),
                        cell.orig=sapply(demuxlet$BEST,function(x){strsplit(x,"-")[[1]][[2]]}))
  row.names(dmx_new) <- dmx_new$barcode
  dmx_new$barcode <- NULL
  i <- which(names(scelist) == pos)
  object<-AddMetaData(scelist[[i]],dmx_new)
  return(object)
})
names(scelist_n)<-c('cx','ladb','rca')

##addmetadata-----
artery_dmx <- merge(scelist_n[[1]], y = scelist_n[2:3])
metadata_info <- artery_dmx@meta.data[, c('cell.quality','cell.orig')]
artery <- AddMetaData(object = artery, metadata = metadata_info)


#pbmc_demuxlet------
demuxlet <- fread("/home/zeyy/dev/sdd/demuxlet/blood.best")
demuxlet<-na.omit(demuxlet)
dmx_new <- data.frame(barcode=demuxlet$BARCODE,
                      cell.quality=sapply(demuxlet$BEST,function(x){strsplit(x,"-")[[1]][[1]]}),
                      cell.orig=sapply(demuxlet$BEST,function(x){strsplit(x,"-")[[1]][[2]]}),
                      rd.total=demuxlet$RD.TOTL,rd.pass=demuxlet$RD.PASS,
                      rd.uniq=demuxlet$RD.UNIQ,n.snp=demuxlet$N.SNP)
row.names(dmx_new) <- dmx_new$barcode
dmx_new$barcode <- NULL
pbmc<-AddMetaData(pbmc,dmx_new)



