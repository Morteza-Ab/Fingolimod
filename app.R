library(DESeq2)
library(readxl)
library(dplyr)
library(pheatmap)
library(data.table)
library(ggplot2)
library(openxlsx) # not needed here but kept for other analysis
library(plotly)
library(htmlwidgets)
library(RColorBrewer) 
library(ggrepel)
library(tidyverse)
library(reshape2)
library(shiny)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
#library(pathfindR)
library(enrichR)
library(biomaRt)
library(RWeka)
library(mice)
library(Hmisc)
library(VIM)
library(readxl)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(UpSetR)
library(writexl)
library(pheatmap)
library(DESeq2)
library(readxl)
library(dplyr)
library(pheatmap)
library(data.table)
library(ggplot2)
library(openxlsx) # not needed here but kept for other analysis
library(plotly)
library(htmlwidgets)
library(RColorBrewer) 
library(ggrepel)
library(tidyverse)
library(reshape2)
library(shiny)
library(DOSE)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
#library(pathfindR)
library(enrichR)
library(biomaRt)

#setwd("/Users/mortezaabyadeh/Desktop/senescence omics data/used for app")


base_url <- "https://github.com/Morteza-Ab/Shiny-Senescence/raw/main/"

transcriptome_Eto_and_CTRL1 <- read.table(paste0(base_url, "clean_transcriptome_Eto_vs_CTRL.txt"), header = FALSE, fill = TRUE)
transcriptome_FTY_and_CTRL1 <- read.table(paste0(base_url, "clean_transcriptome_FTY_vs_CTRL.txt"), header = FALSE, fill = TRUE)
transcriptome_FTY_Eto_and_CTRL1 <- read.table(paste0(base_url, "clean_transcriptome_FTY_Eto_vs_CTRL.txt"), header = FALSE, fill = TRUE)
transcriptome_FTY_and_Eto1 <- read.table(paste0(base_url, "clean_transcriptome_FTY_vs_Eto.txt"), header = FALSE, fill = TRUE)
transcriptome_FTY_Eto_and_Eto1 <- read.table(paste0(base_url, "clean_transcriptome_FTY_Eto_vs_Eto.txt"), header = FALSE, fill = TRUE)
transcriptome_FTY_Eto_and_FTY1 <- read.table(paste0(base_url, "clean_transcriptome_FTY_Eto_vs_FTY.txt"), header = FALSE, fill = TRUE)

proteome_Eto_and_CTRL1 <- read.table(paste0(base_url, "set1.txt"), header = FALSE, fill = TRUE)
proteome_FTY_and_CTRL1 <- read.table(paste0(base_url, "set2.txt"), header = FALSE, fill = TRUE)
proteome_FTY_Eto_and_CTRL1 <- read.table(paste0(base_url, "set3.txt"), header = FALSE, fill = TRUE)
proteome_FTY_and_Eto1 <- read.table(paste0(base_url, "set4.txt"), header = FALSE, fill = TRUE)
proteome_FTY_Eto_and_Eto1 <- read.table(paste0(base_url, "set5.txt"), header = FALSE, fill = TRUE)
proteome_FTY_Eto_and_FTY1 <- read.table(paste0(base_url, "set6.txt"), header = FALSE, fill = TRUE)

transcriptome_Eto_vs_CTRL <- read.xlsx("https://github.com/Morteza-Ab/Shiny-Senescence/raw/main/clean_transcriptome_Eto_vs_CTRL.xlsx")
transcriptome_FTY_vs_CTRL <- read.xlsx("https://github.com/Morteza-Ab/Shiny-Senescence/raw/main/clean_transcriptome_FTY_vs_CTRL.xlsx")
transcriptome_FTY_Eto_vs_CTRL <- read.xlsx("https://github.com/Morteza-Ab/Shiny-Senescence/raw/main/clean_transcriptome_FTY_Eto_vs_CTRL.xlsx")
transcriptome_FTY_vs_Eto <- read.xlsx("https://github.com/Morteza-Ab/Shiny-Senescence/raw/main/clean_transcriptome_FTY_vs_Eto.xlsx")
transcriptome_FTY_Eto_vs_Eto <- read.xlsx("https://github.com/Morteza-Ab/Shiny-Senescence/raw/main/clean_transcriptome_FTY_Eto_vs_Eto.xlsx")
transcriptome_FTY_Eto_vs_FTY <- read.xlsx("https://github.com/Morteza-Ab/Shiny-Senescence/raw/main/clean_transcriptome_FTY_Eto_vs_FTY.xlsx")

proteome_Eto_vs_CTRL <- read.xlsx(paste0(base_url, "set1.xlsx"))
proteome_FTY_vs_CTRL <- read.xlsx(paste0(base_url, "set2.xlsx"))
proteome_FTY_Eto_vs_CTRL <- read.xlsx(paste0(base_url, "set3.xlsx"))
proteome_FTY_vs_Eto <- read.xlsx(paste0(base_url, "set4.xlsx"))
proteome_FTY_Eto_vs_Eto <- read.xlsx(paste0(base_url, "set5.xlsx"))
proteome_FTY_Eto_vs_FTY <- read.xlsx(paste0(base_url, "set6.xlsx"))

colnames(proteome_Eto_vs_CTRL)[which(names(proteome_Eto_vs_CTRL) == "p-adj")] <- "p_adj"
colnames(proteome_FTY_vs_CTRL)[which(names(proteome_FTY_vs_CTRL) == "p-adj")] <- "p_adj"
colnames(proteome_FTY_Eto_vs_CTRL)[which(names(proteome_FTY_Eto_vs_CTRL) == "p-adj")] <- "p_adj"
colnames(proteome_FTY_vs_Eto)[which(names(proteome_FTY_vs_Eto) == "p-adj")] <- "p_adj"
colnames(proteome_FTY_Eto_vs_Eto)[which(names(proteome_FTY_Eto_vs_Eto) == "p-adj")] <- "p_adj"
colnames(proteome_FTY_Eto_vs_FTY)[which(names(proteome_FTY_Eto_vs_FTY) == "p-adj")] <- "p_adj"


datasets1 <- list(
  transcriptome1 = list(
    Eto_and_CTRL = transcriptome_Eto_and_CTRL1,
    FTY_and_CTRL = transcriptome_FTY_and_CTRL1,
    FTY_Eto_and_CTRL = transcriptome_FTY_Eto_and_CTRL1,
    FTY_and_Eto = transcriptome_FTY_and_Eto1,
    FTY_Eto_and_Eto = transcriptome_FTY_Eto_and_Eto1,
    FTY_Eto_and_FTY = transcriptome_FTY_Eto_and_FTY1
  ),
  proteome1 = list(
    Eto_and_CTRL = proteome_Eto_and_CTRL1,
    FTY_and_CTRL = proteome_FTY_and_CTRL1,
    FTY_Eto_and_CTRL = proteome_FTY_Eto_and_CTRL1,
    FTY_and_Eto = proteome_FTY_and_Eto1,
    FTY_Eto_and_Eto = proteome_FTY_Eto_and_Eto1,
    FTY_Eto_and_FTY = proteome_FTY_Eto_and_FTY1
  )
)



########################### bssWssFast function ############################
bssWssFast <- function (X, givenClassArr, numClass=2)
  # between squares / within square feature selection
{
  classVec <- matrix(0, numClass, length(givenClassArr))
  for (k in 1:numClass) {
    temp <- rep(0, length(givenClassArr))
    temp[givenClassArr == (k - 1)] <- 1
    classVec[k, ] <- temp
  }
  classMeanArr <- rep(0, numClass)
  ratio <- rep(0, ncol(X))
  for (j in 1:ncol(X)) {
    overallMean <- sum(X[, j]) / length(X[, j])
    for (k in 1:numClass) {
      classMeanArr[k] <- 
        sum(classVec[k, ] * X[, j]) / sum(classVec[k, ])
    }
    classMeanVec <- classMeanArr[givenClassArr + 1]
    bss <- sum((classMeanVec - overallMean)^2)
    wss <- sum((X[, j] - classMeanVec)^2)
    ratio[j] <- bss/wss
  }
  sort(ratio, decreasing = TRUE, index = TRUE)
}


################################## feature selector function  ###################################

bwf_feature_selection <- function(df, n){
  # n will be selected throu a slider within the app title"number of features" ranging from 10 to 50 interval 5
  mrnaNorm <- df[-1,]
  mrnaIDs <- df[1,]
  mrnaIDs <- mrnaIDs[,-1]
  samp <- lapply(as.list(t(mrnaIDs)), function(t) {
    parts <- unlist(strsplit(t, "_"))
    if (length(parts) > 2) {
      paste(parts[1], parts[2], sep = "_")
    } else {
      parts[1]
    }
  })
  sampleType <- as.data.frame(samp)
  sampClass <- lapply(samp, function(t) (if (t == sampleType[[1]] ) return("1") else return("0")))
  mrnaClass <- as.data.frame(sampClass)
  table(unlist(sampClass))
  sampClassNum <- lapply(samp, function(t) (if (t == sampleType[[1]] ) return(1) else return(0)))
  mrnaClassNum <- as.data.frame(sampClassNum) 
  geneNames <- mrnaNorm[1] # extract the gene names from mrnaNorm as its first column
  mrnaData <- mrnaNorm[, -1]
  mrnaData <- as.data.frame(lapply(mrnaData, as.numeric))
  mrnaData <- na.omit(mrnaData)
  mrnaData = t(mrnaData)
  bss <- bssWssFast(mrnaData, t(mrnaClassNum), 2)
  whole_gene_name <- bss$ix
  gene_name_whole <- geneNames[,1][whole_gene_name]
  data_bss_list <- mrnaData[, whole_gene_name]
  transpose_data_bss_list <- t(data_bss_list)
  colnames(transpose_data_bss_list)
  whole_genes_expression <- cbind(Gene = gene_name_whole, transpose_data_bss_list)
  bss_data <- as.data.frame(bss)
  bss_data$"Gene" <- geneNames[,1][whole_gene_name]
  result <- merge(whole_genes_expression, bss_data, by="Gene")
  sorted <- result[order(-result$x), ]
  colnames(sorted) <- c("Gene", rep(sampleType[[1]], 3), rep(sampleType[[4]], 3), "Score", "index")
  #return(bss_data)
  #return(sorted)
  sorted_subset <- sorted[1:n,2:7]
  rowname <- sorted[1:n, 1]
  numeric_subset <- apply(sorted_subset, 2, function(x) as.numeric(as.character(x)))
  
  # Convert to a matrix
  i_matrix <- as.matrix(numeric_subset)
  
  # Handle any potential NA values (if needed)
  i_matrix[is.na(i_matrix)] <- 0  # or use another method to handle NAs
  
  # Apply the log2 transformation
  #i_log <- log2(i_matrix + 1)
  
  i_log <- log2(i_matrix + 1)
  i <- t(scale(t(i_log), scale = FALSE))
  head(i)
  i_annotation_col <- as.data.frame(t(sampleType))
  rownames(i_annotation_col) <- NULL
  colnames(i)
  dim(i_annotation_col)
  head(i_annotation_col)
  colnames(i)
  rownames(i_annotation_col) <- colnames(i)
  row.names(i) <- rowname
  i_head <- pheatmap(i,
                     cluster_cols = TRUE,  
                     scale = "row",  
                     color = colorRampPalette(c("green", "black", "red"))(100), 
                     main = "Heatmap of Normalized Counts",  
                     fontsize = 13,  
                     fontsize_row = 12,  
                     fontsize_col = 15,  
                     border_color = NA,
                     rownames_force = TRUE,
                     #angle_col = 45,
                     annotation_col = i_annotation_col,
                     width = 10,
                     height = 13
  )
}













datasets <- list(
  transcriptome = list(
    Eto_vs_CTRL = transcriptome_Eto_vs_CTRL,
    FTY_vs_CTRL = transcriptome_FTY_vs_CTRL,
    FTY_Eto_vs_CTRL = transcriptome_FTY_Eto_vs_CTRL,
    FTY_vs_Eto = transcriptome_FTY_vs_Eto,
    FTY_Eto_vs_Eto = transcriptome_FTY_Eto_vs_Eto,
    FTY_Eto_vs_FTY = transcriptome_FTY_Eto_vs_FTY
  ),
  proteome = list(
    Eto_vs_CTRL = proteome_Eto_vs_CTRL,
    FTY_vs_CTRL = proteome_FTY_vs_CTRL,
    FTY_Eto_vs_CTRL = proteome_FTY_Eto_vs_CTRL,
    FTY_vs_Eto = proteome_FTY_vs_Eto,
    FTY_Eto_vs_Eto = proteome_FTY_Eto_vs_Eto,
    FTY_Eto_vs_FTY = proteome_FTY_Eto_vs_FTY
  )
)

# Function for Volcano Plot
volcano_maker <- function(df){
  df <- df[!df$p_adj == 0, ]
  df$diffexpressed <- ifelse(df$`p_adj` < 0.05 & df$`log2FC` > 0, "Up",
                             ifelse(df$`p_adj` < 0.05 & df$`log2FC` < 0, "Down", "No"))
  df$delabel <- df$Gene
  df <- df[complete.cases(df), ]
  counts <- table(df$diffexpressed)
  max_y <- round(max(-log10(df$p_adj)))
  max_x <- round(max(df$log2FC))
  downregulated_count <- counts["Down"]
  not_significant_count <- counts["No"]
  upregulated_count <- counts["Up"]
  
  ggplot(data = df, aes(x = log2FC, y = -log10(p_adj), col = diffexpressed)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point(size = 3) + 
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), 
                       labels = c(paste("Downregulated (", downregulated_count, ")", sep = ""),
                                  paste("Not significant (", not_significant_count, ")", sep = ""),
                                  paste("Upregulated (", upregulated_count, ")", sep = ""))) +
    coord_cartesian(ylim = c(0, max_y), xlim = c(-max_x, max_x)) +
    scale_x_continuous(breaks = seq(-max_x, max_x, 1)) + 
    ggtitle('Volcano plot') +
    theme_bw() + 
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 15)) 
}

# Function for PCA Plot
pca_maker <- function(df){
  df_filtered <- df[df$`p_adj` < 0.05, ]
  numberic_column_dataset1 <- df_filtered[,2:7]
  dataset1_numeric <- data.frame(lapply(numberic_column_dataset1, function(x) if(is.numeric(x)) as.integer(x) else x))
  
  dataset1_numeric <- na.omit(dataset1_numeric)
  colnames(dataset1_numeric)
  #gr <- sub("\\..*", "", colnames(dataset1_numeric))
  gr <- sub("\\-.*", "", colnames(dataset1_numeric))
  gr <- as.factor(gr)
  gr
  levels(gr)
  dataset1_numeric <- na.omit(dataset1_numeric)
  colnames(dataset1_numeric) <- gr
  #rownames(colData) <- colnames(dataset1_numeric)
  colnames(dataset1_numeric) <- gr
  dim(dataset1_numeric)
  cds <- DESeqDataSetFromMatrix(countData = dataset1_numeric, colData <- data.frame(Group = gr), design = ~ Group)
  #cds <- DESeq(cds)
  data.norm <- log2(1+counts(cds, normalized=F))
  data.mean.center <- t(scale(t(data.norm), scale = F))
  pc <- prcomp(data.mean.center)
  pcr <- data.frame(pc$r)
  pcr$group <- gr
  pc_var <- pc$sdev^2 / sum(pc$sdev^2)
  
  ggplot(pcr, aes(PC1, PC2, color=group)) + geom_point(size=10, alpha=0.6, shape = 19) + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15),
          legend.position = "right") + 
    labs(
      x = paste0("PC1 (", round(pc_var[1] * 100, 2), "% variance)"),
      y = paste0("PC2 (", round(pc_var[2] * 100, 2), "% variance)"),
      title = "PCA Plot"
    )
}

# Function for Heatmap
heatmap_maker <- function(df){
  df <- df
  df_filtered <- df[df$`p_adj` < 0.05, ]
  df_filtered <- df_filtered[complete.cases(df_filtered), ]
  i_subset <- df_filtered[,2:7]
  i_matrix <- as.matrix(i_subset)
  i_log <- log2(i_matrix + 1)
  i <- t(scale(t(i_log), scale = FALSE))
  i_annotation_col <- data.frame(Group= c(sub("-.*", "", colnames(i))))
  rownames(i_annotation_col) <- colnames(i)
  i_head <- pheatmap(i,
                     cluster_cols = TRUE,  
                     scale = "row",  
                     color = colorRampPalette(c("green", "black", "red"))(100), 
                     main = "Heatmap of Normalized Counts",  
                     fontsize = 13,  
                     fontsize_row = 12,  
                     fontsize_col = 13,  
                     border_color = NA,  
                     #angle_col = 45,
                     annotation_col = i_annotation_col,
                     width = 10,
                     height = 7
  )
}


# Function for 3D PCA

pca3D_maker <- function(df){
  df_filtered <- df[df$`p_adj` < 0.05, ]
  numberic_column_dataset1 <- df_filtered[,2:7]
  dataset1_numeric <- data.frame(lapply(numberic_column_dataset1, function(x) if(is.numeric(x)) as.integer(x) else x))
  
  dataset1_numeric <- na.omit(dataset1_numeric)
  colnames(dataset1_numeric)
  #gr <- sub("\\..*", "", colnames(dataset1_numeric))
  gr <- sub("\\-.*", "", colnames(dataset1_numeric))
  gr <- as.factor(gr)
  gr
  levels(gr)
  dataset1_numeric <- na.omit(dataset1_numeric)
  colnames(dataset1_numeric) <- gr
  #rownames(colData) <- colnames(dataset1_numeric)
  colnames(dataset1_numeric) <- gr
  dim(dataset1_numeric)
  cds <- DESeqDataSetFromMatrix(countData = dataset1_numeric, colData <- data.frame(Group = gr), design = ~ Group)
  #cds <- DESeq(cds)
  data.norm <- log2(1+counts(cds, normalized=F))
  data.mean.center <- t(scale(t(data.norm), scale = F))
  pc <- prcomp(data.mean.center)
  pcr <- data.frame(pc$r)
  pcr$group <- gr
  pc_var <- pc$sdev^2 / sum(pc$sdev^2)
  plot_ly(data = pcr, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers", color = ~group) %>%
    layout(
      scene = list(
        xaxis = list(title = paste0("PC1 (", round(pc_var[1] * 100, 2), "% variance)")),
        yaxis = list(title = paste0("PC2 (", round(pc_var[2] * 100, 2), "% variance)")),
        zaxis = list(title = paste0("PC3 (", round(pc_var[3] * 100, 2), "% variance)"))
      ),
      title = "3D PCA Plot"
    )
}



# Function for correlation

correlation_maker <- function(df) {
  df1 <- df
  df_filtered1 <- df1[df1$`p_adj` < 0.05, ]
  df_filtered1 <- df_filtered1[complete.cases(df_filtered1), ]
  i_subset1 <- df_filtered1[,2:7]
  i_matrix1 <- as.matrix(i_subset1)
  i_log1 <- log2(i_matrix1 + 1)
  i1 <- t(scale(t(i_log1), scale = FALSE))
  
  corr_matrix <- cor(i1)
  corr_melted <- melt(corr_matrix)
  
  ggplot(data = corr_melted, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab", 
                         name = "Correlation") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    coord_fixed() +
    ggtitle("Correlation") +
    labs(x = "", y = "")
}

### KEGG Down Function 

KEGG_maker_down <- function(df){
  df$diffexpressed <- ifelse(df$`p_adj` < 0.05 & df$`log2FC` > 0, "Up",
                             ifelse(df$`p_adj` < 0.05 & df$`log2FC` < 0, "Down", "No"))
  df$delabel <- df$Gene
  df <- df[complete.cases(df), ]
  df <- df[df$diffexpressed=="Down", ]
  gene <- df$Gene
  gene <- as.data.frame(gene)
  gene <- as.list(gene)
  class(gene)
  listEnrichrDbs()
  dbs <- listEnrichrDbs()
  db_names <- dbs$libraryName
  selected_dbs <- c("KEGG_2021_Human")
  enriched <- enrichr(df$Gene, selected_dbs)
  kegg_result <- as.data.frame(enriched[["KEGG_2021_Human"]])
  kegg_filtered <- kegg_result[kegg_result$Adjusted.P.value < 0.05, ]
  kegg_sorted <- kegg_filtered[order(kegg_filtered$Adjusted.P.value), ]
  top_10_kegg <- head(kegg_sorted, 10)
  ggplot(top_10_kegg, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    xlab("KEGG Pathways") +
    ylab("-log10(Adjusted P-value)") +
    ggtitle("KEGG Pathway Enrichment Analysis") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
}


KEGG_maker_down(transcriptome_Eto_vs_CTRL)

### KEGG Up Function  


KEGG_maker_up <- function(df){
  df$diffexpressed <- ifelse(df$`p_adj` < 0.05 & df$`log2FC` > 0, "Up",
                             ifelse(df$`p_adj` < 0.05 & df$`log2FC` < 0, "Down", "No"))
  df$delabel <- df$Gene
  df <- df[complete.cases(df), ]
  df <- df[df$diffexpressed=="Up", ]
  gene <- df$Gene
  gene <- as.data.frame(gene)
  gene <- as.list(gene)
  class(gene)
  listEnrichrDbs()
  dbs <- listEnrichrDbs()
  db_names <- dbs$libraryName
  selected_dbs <- c("KEGG_2021_Human")
  enriched <- enrichr(df$Gene, selected_dbs)
  kegg_result <- as.data.frame(enriched[["KEGG_2021_Human"]])
  kegg_filtered <- kegg_result[kegg_result$Adjusted.P.value < 0.05, ]
  kegg_sorted <- kegg_filtered[order(kegg_filtered$Adjusted.P.value), ]
  top_10_kegg <- head(kegg_sorted, 10)
  ggplot(top_10_kegg, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
    geom_bar(stat = "identity", fill = "red") +
    coord_flip() +
    xlab("KEGG Pathways") +
    ylab("-log10(Adjusted P-value)") +
    ggtitle("KEGG Pathway Enrichment Analysis") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
}




### Molecular Function Down Function 

MF_maker_down <- function(df){
  df$diffexpressed <- ifelse(df$`p_adj` < 0.05 & df$`log2FC` > 0, "Up",
                             ifelse(df$`p_adj` < 0.05 & df$`log2FC` < 0, "Down", "No"))
  df$delabel <- df$Gene
  df <- df[complete.cases(df), ]
  df1 <- df[df$diffexpressed=="Down", ]
  gene <- df1$Gene
  gene <- as.data.frame(gene)
  gene <- as.list(gene)
  class(gene)
  listEnrichrDbs()
  dbs <- listEnrichrDbs()
  db_names <- dbs$libraryName
  selected_dbs <- c("GO_Molecular_Function_2023")
  enriched <- enrichr(df1$Gene, selected_dbs)
  mf_result <- as.data.frame(enriched[["GO_Molecular_Function_2023"]])
  mf_filtered <- mf_result[mf_result$Adjusted.P.value < 0.05, ]
  mf_sorted <- mf_filtered[order(mf_filtered$Adjusted.P.value), ]
  top_10_mf <- head(mf_sorted, 10)
  
  ggplot(top_10_mf, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    xlab("Molecular Function") +
    ylab("-log10(Adjusted P-value)") +
    ggtitle("Molecular Function Enrichment Analysis") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
}





### Molecular Function Up Function 

MF_maker_up <- function(df){
  df$diffexpressed <- ifelse(df$`p_adj` < 0.05 & df$`log2FC` > 0, "Up",
                             ifelse(df$`p_adj` < 0.05 & df$`log2FC` < 0, "Down", "No"))
  df$delabel <- df$Gene
  df <- df[complete.cases(df), ]
  df1 <- df[df$diffexpressed=="Up", ]
  gene <- df1$Gene
  gene <- as.data.frame(gene)
  gene <- as.list(gene)
  class(gene)
  listEnrichrDbs()
  dbs <- listEnrichrDbs()
  db_names <- dbs$libraryName
  selected_dbs <- c("GO_Molecular_Function_2023")
  enriched <- enrichr(df1$Gene, selected_dbs)
  mf_result <- as.data.frame(enriched[["GO_Molecular_Function_2023"]])
  mf_filtered <- mf_result[mf_result$Adjusted.P.value < 0.05, ]
  mf_sorted <- mf_filtered[order(mf_filtered$Adjusted.P.value), ]
  top_10_mf <- head(mf_sorted, 10)
  ggplot(top_10_mf, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
    geom_bar(stat = "identity", fill = "red") +
    coord_flip() +
    xlab("Molecular Function") +
    ylab("-log10(Adjusted P-value)") +
    ggtitle("Molecular Function Enrichment Analysis") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
}




### Biological Process Down Function 

BP_maker_down <- function(df){
  df$diffexpressed <- ifelse(df$`p_adj` < 0.05 & df$`log2FC` > 0, "Up",
                             ifelse(df$`p_adj` < 0.05 & df$`log2FC` < 0, "Down", "No"))
  df$delabel <- df$Gene
  df <- df[complete.cases(df), ]
  df2 <- df[df$diffexpressed=="Down", ]
  gene <- df2$Gene
  gene <- as.data.frame(gene)
  gene <- as.list(gene)
  class(gene)
  listEnrichrDbs()
  dbs <- listEnrichrDbs()
  db_names <- dbs$libraryName
  selected_dbs <- c("GO_Biological_Process_2023")
  enriched <- enrichr(df2$Gene, selected_dbs)
  bp_result <- as.data.frame(enriched[["GO_Biological_Process_2023"]])
  bp_filtered <- bp_result[bp_result$Adjusted.P.value < 0.05, ]
  bp_sorted <- bp_filtered[order(bp_filtered$Adjusted.P.value), ]
  top_10_bp <- head(bp_sorted, 10)
  ggplot(top_10_bp, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    xlab("Biological Process") +
    ylab("-log10(Adjusted P-value)") +
    ggtitle("Biological Process Analysis") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
}



### Biological Process Up Function 

BP_maker_up <- function(df){
  df$diffexpressed <- ifelse(df$`p_adj` < 0.05 & df$`log2FC` > 0, "Up",
                             ifelse(df$`p_adj` < 0.05 & df$`log2FC` < 0, "Down", "No"))
  df$delabel <- df$Gene
  df <- df[complete.cases(df), ]
  df2 <- df[df$diffexpressed=="Up", ]
  gene <- df2$Gene
  gene <- as.data.frame(gene)
  gene <- as.list(gene)
  class(gene)
  listEnrichrDbs()
  dbs <- listEnrichrDbs()
  db_names <- dbs$libraryName
  selected_dbs <- c("GO_Biological_Process_2023")
  enriched <- enrichr(df2$Gene, selected_dbs)
  bp_result <- as.data.frame(enriched[["GO_Biological_Process_2023"]])
  bp_filtered <- bp_result[bp_result$Adjusted.P.value < 0.05, ]
  bp_sorted <- bp_filtered[order(bp_filtered$Adjusted.P.value), ]
  top_10_bp <- head(bp_sorted, 10)
  ggplot(top_10_bp, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
    geom_bar(stat = "identity", fill = "red") +
    coord_flip() +
    xlab("Biological Process") +
    ylab("-log10(Adjusted P-value)") +
    ggtitle("Biological Process Analysis") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
} 


##################################### App  ####################################


library(shiny)

ui <- navbarPage("Omics Data App",
                 
                 # First Page - Expression Analysis and Enrichment
                 tabPanel("Expression & Enrichment Analysis",
                          fluidPage(
                            sidebarLayout(
                              sidebarPanel(
                                tags$h3("Expression Analysis"),
                                selectInput("data_type", "Select Data Type", choices = c("Transcriptome", "Proteome")),
                                selectInput("dataset", "Select Dataset:", choices = c("Eto_vs_CTRL", "FTY_vs_CTRL", "FTY_Eto_vs_CTRL", "FTY_vs_Eto", "FTY_Eto_vs_Eto", "FTY_Eto_vs_FTY")),
                                uiOutput("dataset_selector"),
                                selectInput("plot_type", "Select Plot Type", choices = c("Volcano Plot", "Heatmap", "PCA Plot", "Correlation Plot"))
                              ),
                              mainPanel(
                                plotOutput("plot", width = "70%", height = "50%")
                              )
                            ),
                            sidebarLayout(
                              sidebarPanel(
                                tags$h3("Enrichment Analysis"),
                                selectInput("dataset_type", "Dataset Type:", choices = c("transcriptome", "proteome")),
                                selectInput("comparison", "Comparison:", choices = c("Eto_vs_CTRL", "FTY_vs_CTRL", "FTY_Eto_vs_CTRL", "FTY_vs_Eto", "FTY_Eto_vs_Eto", "FTY_Eto_vs_FTY")),
                                selectInput("regulation", "Regulation:", choices = c("Down", "Up")),
                                selectInput("enrichment", "Functional Enrichment:", choices = c("KEGG", "Molecular Function", "Biological Process"))
                              ),
                              mainPanel(
                                plotOutput("additional_plot", width = "70%", height = "500px")
                              )
                            )
                          )
                 ),
                 
                 # Second Page - Feature Selection and Heatmap
                 tabPanel("Feature Selection & Heatmap",
                          fluidPage(
                            sidebarLayout(
                              sidebarPanel(
                                tags$h3("Feature Selection"),
                                selectInput("data_type", "Select Data Type", choices = c("Transcriptome", "Proteome")),
                                uiOutput("dataset_selector"),
                                sliderInput("n_features", "Number of features", min = 10, max = 50, step = 5, value = 20)
                              ),
                              mainPanel(
                                plotOutput("heatmap_plot", width = "70%", height = "500px")
                              )
                            )
                          )
                 )
)



server <- function(input, output, session) {
  output$dataset_selector <- renderUI({
    if (input$data_type == "Transcriptome") {
      selectInput("dataset", "Select Dataset", choices = names(datasets$transcriptome))
    } else {
      selectInput("dataset", "Select Dataset", choices = names(datasets$proteome))
    }
  })
  
  output$plot <- renderPlot({
    req(input$data_type, input$dataset, input$plot_type)
    
    df <- if (input$data_type == "Transcriptome") {
      datasets$transcriptome[[input$dataset]]
    } else {
      datasets$proteome[[input$dataset]]
    }
    if (input$plot_type == "Volcano Plot") {
      volcano_maker(df)
    } else if (input$plot_type == "Heatmap") {
      heatmap_maker(df)
    } else if (input$plot_type == "PCA Plot") {
      pca_maker(df)
    } else if (input$plot_type == "Correlation Plot") {
      correlation_maker(df)
    }
  }, height = 600, width = 700, res = 96)
  
  output$additional_plot <- renderPlot({
    req(input$dataset_type, input$comparison, input$regulation, input$enrichment)
    
    dataset <- datasets[[input$dataset_type]][[input$comparison]]
    
    if (input$regulation == "Down" && input$enrichment == "KEGG") {
      KEGG_maker_down(dataset)
    } else if (input$regulation == "Up" && input$enrichment == "KEGG") {
      KEGG_maker_up(dataset)
    } else if (input$regulation == "Down" && input$enrichment == "Molecular Function") {
      MF_maker_down(dataset)
    } else if (input$regulation == "Up" && input$enrichment == "Molecular Function") {
      MF_maker_up(dataset)
    } else if (input$regulation == "Down" && input$enrichment == "Biological Process") {
      BP_maker_down(dataset)
    } else if (input$regulation == "Up" && input$enrichment == "Biological Process") {
      BP_maker_up(dataset)
    }
  })
  
  # Second Page - Feature Selection and Heatmap
  output$dataset_selector <- renderUI({
    if (input$data_type == "Transcriptome") {
      selectInput("dataset", "Select Dataset", choices = names(datasets1$transcriptome1))
    } else {
      selectInput("dataset", "Select Dataset", choices = names(datasets1$proteome1))
    }
  })
  
  output$heatmap_plot <- renderPlot({
    req(input$data_type, input$dataset, input$n_features)
    df <- if (input$data_type == "Transcriptome") {
      datasets1$transcriptome1[[input$dataset]]
    } else {
      datasets1$proteome1[[input$dataset]]
    }
    selected_features <- bwf_feature_selection(df, input$n_features)
  }, height = 1100, res = 96)
}


# Run the application 
shinyApp(ui = ui, server = server)
