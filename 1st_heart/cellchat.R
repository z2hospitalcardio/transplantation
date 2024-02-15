#cellchat first_heart------
library(dplyr)
library(CellChat)

#artery-------
metadata_info <- artery_SNG@meta.data[, c('celltype','cell.orig')]
#combine cell.orig&celltype
metadata_info$label <- paste(metadata_info$celltype, metadata_info$cell.orig, sep = "_")

data.input <- GetAssayData(artery_SNG, assay = "RNA", slot = "data")
data.input <- data.input[,rownames(metadata_info) ]

cellchat_artery <- createCellChat(object = data.input, meta = metadata_info, group.by = "label")

cellchat_artery@meta[["label"]] <- sub("_Recipient$", "_Re", cellchat_artery@meta[["label"]])
cellchat_artery@meta[["label"]] <- sub("_Donor$", "_Do", cellchat_artery@meta[["label"]])
cellchat_artery@idents <- as.factor(cellchat_artery@meta[['label']])

##preprocess------
CellChatDB <- CellChatDB.human  
showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# set the used database in the object
cellchat_artery@DB <- CellChatDB.use

#recognition ligand receptor interaction
cellchat_artery <- subsetData(cellchat_artery) 
future::plan("multisession", workers = 4) # do parallel
cellchat_artery <- identifyOverExpressedGenes(cellchat_artery)
cellchat_artery <- identifyOverExpressedInteractions(cellchat_artery)

#calculating probabilities
cellchat_artery <- computeCommunProb(cellchat_artery,population.size = TRUE)#Compute the communication probability/strength,in net slot
cellchat_artery <- filterCommunication(cellchat_artery, min.cells = 10)
cellchat_artery <- computeCommunProbPathway(cellchat_artery)#in netP slot
#The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively

df.net <- subsetCommunication(cellchat_artery)
write.csv(df.net, "/home/zeyy/dev/sdd/pyh/result/cellchat/artery/artery_cc_communications.all.csv")

df.net <- subsetCommunication(cellchat_artery,slot.name='netP')
write.csv(df.net, "/home/zeyy/dev/sdd/pyh/result/cellchat/artery/artery_sig_pathway.all.csv")


##network-------
#calculate the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability
cellchat_artery <- aggregateNet(cellchat_artery)

#visualization
groupSize <- as.numeric(table(cellchat_artery@idents))
pdf('/home/zeyy/dev/sdd/pyh/result/cellchat/artery/artery_number.pdf')
netVisual_circle(cellchat_artery@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf('/home/zeyy/dev/sdd/pyh/result/cellchat/artery/artery_weight.pdf')
netVisual_circle(cellchat_artery@net$weight, vertex.weight = groupSize, 
                 weight.scale = T,  title.name = "Interaction weights/strength")
dev.off()

saveRDS(cellchat_artery,'/home/zeyy/dev/sdd/pyh/data/cellchat_artery.rds')

##LR------
# Get Weight Matrix
weight_matrix <- cellchat_artery@net[["weight"]]
df <- as.data.frame(as.table(weight_matrix))
names(df) <- c("Cell1", 'Cell2',"Weight")
top3 <- df %>% arrange(desc(Weight)) %>% head(3)
print(top3)
top5 <- df %>% arrange(desc(Weight)) %>% head(5)
print(top5)
# Get probability Matrix
prob_matrix <- cellchat_artery@net[["prob"]]
df <- as.data.frame(as.table(prob_matrix))
names(df) <- c("Cell1", "Cell2",'LR','Prob')
top10_LR_Myeloid_Re_Myeloid_Re <- df %>%
  filter(Cell1 == "Myeloid_Re", Cell2 == "Myeloid_Re") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
top10_LR_Fibroblast_Do_Myeloid_Re <- df %>%
  filter(Cell1 == "Fibroblast_Do", Cell2 == "Myeloid_Re") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
top10_LR_Myeloid_Re_Endothelial_Do <- df %>%
  filter(Cell1 == "Myeloid_Re", Cell2 == "Endothelial_Do") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
combined_LR <- c(top10_LR_Myeloid_Re_Myeloid_Re, top10_LR_Fibroblast_Do_Myeloid_Re, top10_LR_Myeloid_Re_Endothelial_Do)
combined_LR <- as.character(combined_LR)
unique_LR <- unique(combined_LR)
unique_LR <- data.frame(interaction_name = unique_LR)
pdf('/home/zeyy/dev/sdd/pyh/result/cellchat/artery/LR.pdf')
netVisual_chord_gene(cellchat_artery, sources.use = c('Myeloid_Re','Fibroblast_Do'), 
                     targets.use = c('Myeloid_Re','Endothelial_Do'), lab.cex = 0.5,legend.pos.y = 30,pairLR.use = unique_LR)
dev.off()


#artium--------
metadata_info <- atrium@meta.data[, c('celltype','cell.orig')]
metadata_info$label <- paste(metadata_info$celltype, metadata_info$cell.orig, sep = "_")

data.input <- GetAssayData(atrium, assay = "RNA", slot = "data")
data.input <- data.input[,rownames(metadata_info) ]

celchat_atrium <- createCellChat(object = data.input, meta = metadata_info, group.by = "label")

celchat_atrium@meta[["label"]] <- sub("_Recipient$", "_Re", celchat_atrium@meta[["label"]])
celchat_atrium@meta[["label"]] <- sub("_Donor$", "_Do", celchat_atrium@meta[["label"]])
celchat_atrium@idents <- as.factor(celchat_atrium@meta[['label']])

##preprocess------
celchat_atrium@DB <- CellChatDB.use

celchat_atrium <- subsetData(celchat_atrium) 
celchat_atrium <- identifyOverExpressedGenes(celchat_atrium)
celchat_atrium <- identifyOverExpressedInteractions(celchat_atrium)

celchat_atrium <- computeCommunProb(celchat_atrium,population.size = TRUE)
celchat_atrium <- filterCommunication(celchat_atrium, min.cells = 10)
celchat_atrium <- computeCommunProbPathway(celchat_atrium)

df.net <- subsetCommunication(celchat_atrium)
write.csv(df.net, "/home/zeyy/dev/sdd/pyh/result/cellchat/atrium/atrium_cc_communications.all.csv")

df.net <- subsetCommunication(celchat_atrium,slot.name='netP')
write.csv(df.net, "/home/zeyy/dev/sdd/pyh/result/cellchat/atrium/atrium_sig_pathway.all.csv")

##network--------
celchat_atrium <- aggregateNet(celchat_atrium)

#visualization
groupSize <- as.numeric(table(celchat_atrium@idents))
pdf('/home/zeyy/dev/sdd/pyh/result/cellchat/atrium/atrium_number.pdf')
netVisual_circle(celchat_atrium@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf('/home/zeyy/dev/sdd/pyh/result/cellchat/atrium/atrium_weight.pdf')
netVisual_circle(celchat_atrium@net$weight, vertex.weight = groupSize, 
                 weight.scale = T,  title.name = "Interaction weights/strength")
dev.off()

saveRDS(celchat_atrium,'/home/zeyy/dev/sdd/pyh/data/cellchat_atrium.rds')

##LR----
weight_matrix <- cellchat_atrium@net[["weight"]]
df <- as.data.frame(as.table(weight_matrix))
names(df) <- c("Cell1", 'Cell2',"Weight")
top3 <- df %>% arrange(desc(Weight)) %>% head(3)
print(top3)
top5 <- df %>% arrange(desc(Weight)) %>% head(5)
print(top5)
prob_matrix <- cellchat_atrium@net[["prob"]]
df <- as.data.frame(as.table(prob_matrix))
names(df) <- c("Cell1", "Cell2",'LR','Prob')
top10_LR_Fibroblast_Do_Myeloid_Re <- df %>%
  filter(Cell1 == "Fibroblast_Do", Cell2 == "Myeloid_Re") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
top10_LR_Endothelial_Do_Myeloid_Re <- df %>%
  filter(Cell1 == "Endothelial_Do", Cell2 == "Myeloid_Re") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
top10_LR_Myeloid_Re_Myeloid_Re <- df %>%
  filter(Cell1 == "Myeloid_Re", Cell2 == "Myeloid_Re") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
combined_LR <- c(top10_LR_Fibroblast_Do_Myeloid_Re, top10_LR_Endothelial_Do_Myeloid_Re, top10_LR_Myeloid_Re_Myeloid_Re)
combined_LR <- as.character(combined_LR)
unique_LR <- unique(combined_LR)#去重
unique_LR <- data.frame(interaction_name = unique_LR)
pdf('/home/zeyy/dev/sdd/pyh/result/cellchat/atrium/LR.pdf')
netVisual_chord_gene(cellchat_atrium, sources.use = c('Fibroblast_Do','Endothelial_Do','Myeloid_Re'), 
                     targets.use = 'Myeloid_Re', lab.cex = 0.5,legend.pos.y = 30,pairLR.use = unique_LR)
dev.off()


#ventricle-------
metadata_info <- ventricle@meta.data[, c('celltype','cell.orig')]
metadata_info$label <- paste(metadata_info$celltype, metadata_info$cell.orig, sep = "_")

data.input <- GetAssayData(ventricle, assay = "RNA", slot = "data")
data.input <- data.input[,rownames(metadata_info) ]
#
celchat_ventricle <- createCellChat(object = data.input, meta = metadata_info, group.by = "label")

celchat_ventricle@meta[["label"]] <- sub("_Recipient$", "_Re", celchat_ventricle@meta[["label"]])
celchat_ventricle@meta[["label"]] <- sub("_Donor$", "_Do", celchat_ventricle@meta[["label"]])
celchat_ventricle@idents <- as.factor(celchat_ventricle@meta[['label']])

##preprocess------
celchat_ventricle@DB <- CellChatDB.use

#识别配受体互作
celchat_ventricle <- subsetData(celchat_ventricle) 
celchat_ventricle <- identifyOverExpressedGenes(celchat_ventricle)
celchat_ventricle <- identifyOverExpressedInteractions(celchat_ventricle)

#计算概率
celchat_ventricle <- computeCommunProb(celchat_ventricle,population.size = TRUE)
celchat_ventricle <- filterCommunication(celchat_ventricle, min.cells = 10)
celchat_ventricle <- computeCommunProbPathway(celchat_ventricle)

df.net <- subsetCommunication(celchat_ventricle)
write.csv(df.net, "/home/zeyy/dev/sdd/pyh/result/cellchat/ventricle/ventricle_cc_communications.all.csv")

df.net <- subsetCommunication(celchat_ventricle,slot.name='netP')
write.csv(df.net, "/home/zeyy/dev/sdd/pyh/result/cellchat/ventricle/ventricle_sig_pathway.all.csv")

##network--------
celchat_ventricle <- aggregateNet(celchat_ventricle)

#visualization
groupSize <- as.numeric(table(celchat_ventricle@idents))
pdf('/home/zeyy/dev/sdd/pyh/result/cellchat/ventricle/ventricle_number.pdf')
netVisual_circle(celchat_ventricle@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf('/home/zeyy/dev/sdd/pyh/result/cellchat/ventricle/ventricle_weight.pdf')
netVisual_circle(celchat_ventricle@net$weight, vertex.weight = groupSize, 
                 weight.scale = T,  title.name = "Interaction weights/strength")
dev.off()

saveRDS(celchat_ventricle,'/home/zeyy/dev/sdd/pyh/data/cellchat_ventricle.rds')

##LR------
weight_matrix <- cellchat_ventricle@net[["weight"]]
df <- as.data.frame(as.table(weight_matrix))
names(df) <- c("Cell1", 'Cell2',"Weight")
top3 <- df %>% arrange(desc(Weight)) %>% head(3)
print(top3)
top5 <- df %>% arrange(desc(Weight)) %>% head(5)
print(top5)
prob_matrix <- cellchat_ventricle@net[["prob"]]
df <- as.data.frame(as.table(prob_matrix))
names(df) <- c("Cell1", "Cell2",'LR','Prob')
top10_LR_Endothelial_Re_Endothelial_Re <- df %>%
  filter(Cell1 == "Endothelial_Re", Cell2 == "Endothelial_Re") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
top10_LR_Fibroblast_Do_Endothelial_Re <- df %>%
  filter(Cell1 == "Fibroblast_Do", Cell2 == "Endothelial_Re") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
top10_LR_SMC_Do_Endothelial_Re <- df %>%
  filter(Cell1 == "SMC_Do", Cell2 == "Endothelial_Re") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
combined_LR <- c(top10_LR_Endothelial_Re_Endothelial_Re, top10_LR_Fibroblast_Do_Endothelial_Re, top10_LR_SMC_Do_Endothelial_Re)
combined_LR <- as.character(combined_LR)
unique_LR <- unique(combined_LR)#去重
unique_LR <- data.frame(interaction_name = unique_LR)
pdf('/home/zeyy/dev/sdd/pyh/result/cellchat/ventricle/LR.pdf')
netVisual_chord_gene(cellchat_ventricle, sources.use = c('Endothelial_Re','Fibroblast_Do','SMC_Do'), 
                     targets.use = 'Endothelial_Re', lab.cex = 0.5,legend.pos.y = 30,pairLR.use = unique_LR)
dev.off()


#pbmc------
metadata_info <- pbmc_SNG@meta.data[, c('celltype','cell.orig')]
metadata_info$label <- paste(metadata_info$celltype, metadata_info$cell.orig, sep = "_")

data.input <- GetAssayData(pbmc_SNG, assay = "RNA", slot = "data")
data.input <- data.input[,rownames(metadata_info) ]

celchat_pbmc <- createCellChat(object = data.input, meta = metadata_info, group.by = "label")

celchat_pbmc@meta[["label"]] <- sub("_Recipient$", "_Re", celchat_pbmc@meta[["label"]])
celchat_pbmc@meta[["label"]] <- sub("_Donor$", "_Do", celchat_pbmc@meta[["label"]])
celchat_pbmc@idents <- as.factor(celchat_pbmc@meta[['label']])

##preprocess------
celchat_pbmc@DB <- CellChatDB.use

celchat_pbmc <- subsetData(celchat_pbmc) 
celchat_pbmc <- identifyOverExpressedGenes(celchat_pbmc)
celchat_pbmc <- identifyOverExpressedInteractions(celchat_pbmc)

#计算概率
celchat_pbmc <- computeCommunProb(celchat_pbmc,population.size = TRUE)
celchat_pbmc <- filterCommunication(celchat_pbmc, min.cells = 10)
celchat_pbmc <- computeCommunProbPathway(celchat_pbmc)
df.net <- subsetCommunication(celchat_pbmc)
write.csv(df.net, "/home/zeyy/dev/sdd/pyh/result/cellchat/pbmc/pbmc_cc_communications.all.csv")

df.net <- subsetCommunication(celchat_pbmc,slot.name='netP')
write.csv(df.net, "/home/zeyy/dev/sdd/pyh/result/cellchat/pbmc/pbmc_sig_pathway.all.csv")

##network--------
celchat_pbmc <- aggregateNet(celchat_pbmc)

#visualization
groupSize <- as.numeric(table(celchat_pbmc@idents))
pdf('/home/zeyy/dev/sdd/pyh/result/cellchat/pbmc/pbmc_number.pdf')
netVisual_circle(celchat_pbmc@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf('/home/zeyy/dev/sdd/pyh/result/cellchat/pbmc/pbmc_weight.pdf')
netVisual_circle(celchat_pbmc@net$weight, vertex.weight = groupSize, 
                 weight.scale = T,  title.name = "Interaction weights/strength")
dev.off()

saveRDS(pbmc_SNG,'/home/zeyy/dev/sdd/pyh/data/cellchat_pbmc.rds')

##LR------
weight_matrix <- cellchat_pbmc@net[["weight"]]
df <- as.data.frame(as.table(weight_matrix))
names(df) <- c("Cell1", 'Cell2',"Weight")
top3 <- df %>% arrange(desc(Weight)) %>% head(3)
print(top3)
top5 <- df %>% arrange(desc(Weight)) %>% head(5)
print(top5)
prob_matrix <- cellchat_pbmc@net[["prob"]]
df <- as.data.frame(as.table(prob_matrix))
names(df) <- c("Cell1", "Cell2",'LR','Prob')
top10_LR_CD4_T_Re_MonoMacro_Re <- df %>%
  filter(Cell1 == "CD4+_T_Re", Cell2 == "Mono/Macro_Re") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
top10_LR_CD4_T_Re_B_Re <- df %>%
  filter(Cell1 == "CD4+_T_Re", Cell2 == "B_Re") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
top10_LR_CD8_T_Re_MonoMacro_Re <- df %>%
  filter(Cell1 == "CD8+_T_Re", Cell2 == "Mono/Macro_Re") %>%
  arrange(desc(Prob)) %>%
  head(10) %>%
  pull(LR)
combined_LR <- c(top10_LR_CD4_T_Re_MonoMacro_Re, top10_LR_CD4_T_Re_B_Re, top10_LR_CD8_T_Re_MonoMacro_Re)
combined_LR <- as.character(combined_LR)
unique_LR <- unique(combined_LR)#去重
unique_LR <- data.frame(interaction_name = unique_LR)
pdf('/home/zeyy/dev/sdd/pyh/result/cellchat/pbmc/LR.pdf')
netVisual_chord_gene(cellchat_pbmc, sources.use = c('CD4+_T_Re','CD8+_T_Re'), 
                     targets.use = c('Mono/Macro_Re','B_Re'), lab.cex = 0.5,legend.pos.y = 30,pairLR.use = unique_LR)
dev.off()




