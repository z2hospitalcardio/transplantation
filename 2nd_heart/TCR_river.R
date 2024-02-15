#tcr-----
library(scRepertoire)

#second heart------
pbmc=read.csv("/home/zeyy/dev/sdd/second_human/TCR_BCR/TCR/PBMC/filtered_contig_annotations.csv")
cx=read.csv("/home/zeyy/dev/sdd/second_human/TCR_BCR/TCR/CX/filtered_contig_annotations.csv")
lad=read.csv("/home/zeyy/dev/sdd/second_human/TCR_BCR/TCR/LAD/filtered_contig_annotations.csv")
lvb=read.csv("/home/zeyy/dev/sdd/second_human/TCR_BCR/TCR/LVB/filtered_contig_annotations.csv")
lvpw=read.csv("/home/zeyy/dev/sdd/second_human/TCR_BCR/TCR/LVPW/filtered_contig_annotations.csv")

contig_list <- list(pbmc, cx, lad, lvb,lvpw)
combined1 <- combineTCR(contig_list, 
                       samples = c("blood", "cx", "lad", "lvb", "lvpw"), 
                       ID = c("re", "do-a", "do-a", "do-h", "do-h")
)

compareClonotypes(combined1, numbers = 5, 
                  samples = c("blood_re","cx_do-a","lad_do-a","lvb_do-h","lvpw_do-h"), 
                  cloneCall="aa", graph = "alluvial")

pdf('/home/zeyy/dev/sdd/pyh/result/TCR/second_heart.pdf',width=11.5)
compareClonotypes_mine(combined1,
                  samples = c("blood_re","cx_do-a","lad_do-a","lvb_do-h","lvpw_do-h"), 
                  cloneCall="aa", graph = "alluvial",n=5)+
  theme(axis.text.x = element_text(size=20),axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size = 14),
        legend.title = element_text(size=18))
dev.off()


#Modify the function source code
trace(compareClonotypes,edit=T)

#compareClonotypes_mine-----
compareClonotypes_mine<-function (df, cloneCall = "strict", chain = "both", samples = NULL, 
                                  clonotypes = NULL, numbers = NULL, split.by = NULL, graph = "alluvial", 
                                  exportTable = FALSE, n = 10) 
{
  library(magrittr)
  library(ggplot2)
  library(ggalluvial)
  #list.input.return
  list.input.return<-function (df, split.by) 
  {
    if (inherits(x = df, what = "Seurat") | inherits(x = df, 
                                                     what = "SummarizedExperiment")) {
      if (is.null(split.by)) {
        split.by <- "cluster"
      }
      df <- expression2List(df, split.by)
    }
    return(df)
  }
  #theCall
  theCall<-function (x) 
  {
    if (x %in% c("CTnt", "CTgene", "CTaa", "CTstrict")) {
      x <- x
    }
    else if (x == "gene" | x == "genes") {
      x <- "CTgene"
    }
    else if (x == "nt" | x == "nucleotide") {
      x <- "CTnt"
    }
    else if (x == "aa" | x == "amino") {
      x <- "CTaa"
    }
    else if (x == "gene+nt" | x == "strict") {
      x <- "CTstrict"
    }
    return(x)
  }
  #checkBlanks
  checkBlanks<-function (df, cloneCall) 
  {
    count <- NULL
    for (i in seq_along(df)) {
      if (length(df[[i]][, cloneCall]) == length(which(is.na(df[[i]][, 
                                                                     cloneCall]))) | length(which(!is.na(df[[i]][, cloneCall]))) == 
          0 | nrow(df[[i]]) == 0) {
        count <- c(i, count)
      }
      else {
        (next)()
      }
    }
    if (!is.null(count)) {
      df <- df[-count]
    }
    return(df)
  }
  
  df <- list.input.return(df, split.by)
  cloneCall <- theCall(cloneCall)
  df <- checkBlanks(df, cloneCall)
  if (!is.null(numbers) & !is.null(clonotypes)) {
    stop("Make sure your inputs are either numbers or clonotype sequences.")
  }
  Con.df <- NULL
  for (i in seq_along(df)) {
    if (chain != "both") {
      df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
    }
    tbl <- as.data.frame(table(df[[i]][, cloneCall]))
    tbl[, 2] <- tbl[, 2]/sum(tbl[, 2])
    colnames(tbl) <- c("Clonotypes", "Proportion")
    tbl$Sample <- names(df[i])
    Con.df <- rbind.data.frame(Con.df, tbl)
  }
  if (!is.null(samples)) {
    Con.df <- Con.df[Con.df$Sample %in% samples, ]
  }
  if (!is.null(clonotypes)) {
    Con.df <- Con.df[Con.df$Clonotypes %in% clonotypes, 
    ]
  }
  if (!is.null(numbers)) {
    top <- Con.df %>% dplyr::group_by(Con.df[, 3]) %>% dplyr::slice_max(n = numbers, 
                                                          order_by = Proportion, with_ties = FALSE)
    Con.df <- Con.df[Con.df$Clonotypes %in% top$Clonotypes, 
    ]
  }
  if (nrow(Con.df) < length(unique(Con.df$Sample))) {
    stop("Reasses the filtering strategies here, there is not \n            enough clonotypes to examine.")
  }
  if (exportTable == TRUE) {
    return(Con.df)
  }
  
  
  #Con.df <- Con.df[order(Con.df$Proportion), ]
  Con.df <- dplyr::arrange(Con.df, desc(Proportion))

  Con.df$Clonotypes <- factor(Con.df$Clonotypes, levels = unique(Con.df$Clonotypes))
  #Con.df$Clonotypes <- factor(Con.df$Clonotypes, levels = c(hyper_levels, other_levels))
  
  #find top n oligoclonal，and name other cells as 'others'
  hyper <- unique(Con.df$Clonotypes)[1:n]
  
  Con.df$Clonotype_group <- ifelse(Con.df$Clonotypes %in% hyper, Con.df$Clonotypes, "Other")
  
  color_palette <- c("#ED9A7C","#EBB365","#74C476","#4DBBD5FF","#00468BFF", rep("grey", nrow(Con.df) - n))
  
  plot <- ggplot(Con.df, aes(x = Sample, fill = Clonotypes, 
                             group = Clonotypes, stratum = Clonotypes, alluvium = Clonotypes, 
                             y = Proportion, label = Clonotypes)) + theme_classic() + 
    theme(axis.title.x = element_blank())+
    scale_fill_manual(values = color_palette,breaks = hyper)
  
  
  if (graph == "alluvial") {
    plot <- plot + geom_stratum() + geom_flow(stat = "alluvium")
  }
  else if (graph == "area") {
    plot <- plot + geom_area(aes(group = Clonotypes), color = "black")
  }
  return(plot)
  
}

  
