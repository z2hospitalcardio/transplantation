library(magrittr)
library(dplyr)
library(ggalluvial)
col21 = c("#92D0E6","#F5949B","#E11E24","#FBB96F","#007AB8","#A2D184","#00A04E",
          "#BC5627","#0080BD","#EBC379","#A74D9D","#F57E21","#00998F","#F4E708",
          "#F184B3","#BCB9DA","#6CC7C0","#F37C75","#5FB0D4","#F9AD62",rep('grey',nrow(df)-20))


#first heart:pbmc_sc,artery_bulk,heart_bulk-------
##TRB------
#before_pbmc_sc<-read.csv("/home/zeyy/dev/sdf/human/vdj/tcr/outs/filtered_contig_annotations.csv")
artery_bulk<-read.csv("/home/zeyy/dev/sdf/human/vdj/TCR测序/2023-12-26/浙江大学附属第二医院TRA测序/result/01.Excel/SP_A.csv")
heart_bulk<-read.csv("/home/zeyy/dev/sdf/human/vdj/TCR测序/2023-11-14/浙江大学附属第二医院TRA测序/result/01.Excel/arteryAP_A.csv")

#pbmc
pbmc_sc<-subset(before_pbmc_sc,subset=chain=='TRB')
pbmc_sc<-pbmc_sc[,c('v_gene','d_gene','j_gene','c_gene','cdr3','cdr3_nt')]
df_pbmc_sc<-pbmc_sc%>%pull(cdr3)%>%table()
df_pbmc_sc<-df_pbmc_sc[order(df_pbmc_sc,decreasing = T)]
df_pbmc_sc<-data.frame(sample=rep('pbmc_sc',length(df_pbmc_sc)),value=df_pbmc_sc)
names(df_pbmc_sc) <- c("sample", "CDR3.aa", "Clone.Count")
df_pbmc_sc$Clone.Proportion<-df_pbmc_sc$Clone.Count/sum(df_pbmc_sc$Clone.Count)*100

#artery
df_artery_bulk<-artery_bulk[,c('CDR3.AA.Sequence','Clone.Count')]
df_artery_bulk <- df_artery_bulk %>%
  group_by(`CDR3.AA.Sequence`) %>%
  summarise(`Clone.Count` = sum(`Clone.Count`)) %>%
  arrange(desc(`Clone.Count`))
df_artery_bulk$sample<-rep('artery_bulk',nrow(df_artery_bulk))
df_artery_bulk<-cbind(df_artery_bulk[,3],df_artery_bulk[1:2])
names(df_artery_bulk) <- c("sample", "CDR3.aa", "Clone.Count")
df_artery_bulk$Clone.Proportion<-df_artery_bulk$Clone.Count/sum(df_artery_bulk$Clone.Count)*100
df_artery_bulk$CDR3.aa<-factor(df_artery_bulk$CDR3.aa,levels = df_artery_bulk$CDR3.aa)

#heart
df_heart_bulk<-heart_bulk[,c('CDR3.AA.Sequence','Clone.Count')]

df_heart_bulk <- df_heart_bulk %>%
  group_by(`CDR3.AA.Sequence`) %>%
  summarise(`Clone.Count` = sum(`Clone.Count`)) %>%
  arrange(desc(`Clone.Count`))
df_heart_bulk$sample<-rep('heart_bulk',nrow(df_heart_bulk))
df_heart_bulk<-cbind(df_heart_bulk[,3],df_heart_bulk[1:2])
names(df_heart_bulk) <- c("sample", "CDR3.aa", "Clone.Count")
df_heart_bulk$Clone.Proportion<-df_heart_bulk$Clone.Count/sum(df_heart_bulk$Clone.Count)*100
df_heart_bulk$CDR3.aa<-factor(df_heart_bulk$CDR3.aa,levels = df_heart_bulk$CDR3.aa)

df<-rbind(df_pbmc_sc,df_artery_bulk,df_heart_bulk)
df$sample<-factor(df$sample,levels=c('pbmc_sc','artery_bulk','heart_bulk'))

hyper<-df$CDR3.aa[1:20]

pdf("/home/zeyy/dev/sdd/pyh/result/TCR/first_heart_all3_TRB_aa.pdf")
p<-ggplot(df, aes(x = sample, fill = CDR3.aa, 
                  group = CDR3.aa, stratum = CDR3.aa, alluvium = CDR3.aa,
                  y = Clone.Proportion, label = CDR3.aa)) + theme_classic() +
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = col21,breaks = hyper)+
  geom_stratum() + geom_flow(stat = "alluvium")
p
dev.off()







