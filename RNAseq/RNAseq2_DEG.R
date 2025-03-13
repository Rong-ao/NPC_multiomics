library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(stringr)
library(magrittr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
############# Differential Expression Gene analysis #############

ori_counts = read_csv("/RNAseq2/raw_counts.csv") # path to RNAseq2 dataset
grouplist = c('D1','D1','D2','D2','D3_Ctrl','D3_Ctrl','D3_HU48h','D3_HU48h','D3_HU24h','D3_HU24h','D4_Ctrl','D4_Ctrl')
col_data <- data.frame(
  sampleName = colnames(ori_counts)[7:length(colnames(ori_counts))],
  group_list = grouplist  # Adjust to your experimental design
)
counts_data <- ori_counts[, c(1, 7:length(colnames(ori_counts)))]

exp_list <- c("D2",
              "D3_Ctrl",
              "D3_HU24h",
              "D3_HU48h",
              "D3_HU24h",
              "D3_HU48h",
              "D3_HU24h")

ctr_list <- c("D1",
              "D2",
              "D3_HU48h",
              "D3_Ctrl",
              "D3_Ctrl",
              "D4_Ctrl",
              "D4_Ctrl")

for (i in c(1:7)) {
  exp <- exp_list[i]
  ctr <- ctr_list[i]
  log2FC = 1.5
  padj = 0.05
  col_data2 <- subset(col_data, group_list %in% c(exp, ctr))
  col_data2$group_list <- factor(col_data2$group_list)
  dds <- DESeqDataSetFromMatrix(countData = counts_data[,c(1, which(colnames(counts_data) %in% col_data2$sampleName))], 
                                colData = col_data2,
                                design = ~ group_list, tidy = T)
  keep <- rowSums(counts(dds)) >= 1.5*ncol(counts_data)
  dds <- dds[keep,]
  dds <- DESeq(dds, quiet=F)
  result <- results(dds, contrast = c("group_list", exp, ctr))
  resOrdered <- result[order(result$padj),]
  tempDEG <- as.data.frame(resOrdered)
  DEG <-subset(tempDEG,padj < 0.05 & (log2FoldChange > log2FC | log2FoldChange < -log2FC))
  nona_DEseq2 <- na.omit(tempDEG)
  result1 <- nona_DEseq2[c(2,6)]
  result1[,c(1,2)] <- as.numeric(unlist(result1[,c(1,2)]))
  result1$change = ifelse(result1$padj < 0.05 & abs(result1$log2FoldChange) >= log2FC, ifelse(result1$log2FoldChange> log2FC ,'Up','Down'),'Stable')
  
  
  g <- ggplot(data=result1, aes(x=result1$log2FoldChange, y=-log10(result1$padj), color=change)) +
    geom_point(alpha=0.8, size=0.8) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
    geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
    xlab(bquote(Log[2]*FoldChange))+
    ylab(bquote(-Log[10]*italic(P.adj)))+
    labs(title=paste0(exp, ' vs. ', ctr))+
    theme_classic(base_size = 14) +
    scale_color_manual('',labels=c(paste0("Down(",table(result1$change)[[1]],')'),'Stable',
                                   paste0("Up(",table(result1$change)[[3]],')' )),
                       values=c("blue", "grey","red") )+
    guides(color=guide_legend(override.aes = list(size=3, alpha=1)))
  ggsave(paste0("/RNAseq2/DEG/", exp, ".vs.", ctr, "_volcano.pdf"),g, width=8,height=6)
  
  result1$Label = ""
  upordown <- subset(result1, change=="Up"|change=='Down')
  upordown.gene <- c(as.character(rownames(upordown)))
  result1$Label[match(upordown.gene, rownames(result1))] <- upordown.gene
  upordown$Label[match(upordown.gene, rownames(upordown))] <- upordown.gene
  write.csv(upordown,paste0('/RNAseq2/DEG/', exp, ".vs.", ctr, '.csv'),row.names = TRUE)
}


############ Separeated GO enrichment #############

for (i in c(1:7)) {
  exp <- exp_list[i]
  ctr <- ctr_list[i]
  upordown <- read.csv(paste0('/RNAseq2/DEG/', exp, ".vs.", ctr, '.csv'))
  down <- subset(upordown, change=="Down")
  ids=bitr(down$X, 'SYMBOL','ENTREZID','org.Hs.eg.db')
  xx2 = enrichGO(ids$ENTREZID,
                 OrgDb = "org.Hs.eg.db",   # Annotation database must correspond to species origin
                 keyType = 'ENTREZID',      
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = T)
  xx2_df <- xx2@result
  write.csv(xx2_df,paste0('/RNAseq2/DEG/down_', exp, ".vs.", ctr, '_GO(BP).csv'),row.names = TRUE)
  # significant_threshold <- 0.05
  significant_threshold <- 1  # for checking all possible enrichment pathways
  significant_down  <- xx2_df[xx2_df$p.adjust < significant_threshold, ]

  # Sort the significant data by p.adjust (or another metric of choice)
  significant_down  <- significant_down[order(significant_down$p.adjust), ]
  # Choose how many top terms to display
  top_n <- 10
  top_terms_down  <- head(significant_down , top_n)
  top_terms_down$neg_log_pvalue <- -log(top_terms_down$pvalue, base = 10)
  top_terms_down$Description <- reorder(top_terms_down$Description, top_terms_down$neg_log_pvalue)
  top_terms_down$Description <- sapply(top_terms_down$Description, function(x) str_wrap(x, width = 30)) 
  top_terms_down$Description <- reorder(top_terms_down$Description, top_terms_down$neg_log_pvalue)
  # Plotting with ggplot2
  plot <- ggplot(top_terms_down, aes(y = Description, x = neg_log_pvalue)) +
    geom_bar(stat = "identity", fill = '#000080') + # navy blue
    #  coord_flip() +  # Flip the coordinates to make it a horizontal bar plot
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),  # White background
      axis.line = element_line(color = "black"),  # Solid black axis lines
      axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),  # Adjust text angle if necessary
      axis.text.y = element_text(size = 16, colour = "black"),  # Increase size of y-axis labels
      axis.ticks = element_line(color = "black")  # Black axis ticks
    ) +
    labs(x = "-log10(p-value)", y = "", title = "Significant GO Terms for Down-regulated Genes \n (org.Hs.eg.db)") +
    scale_x_continuous(expand = c(0, 0))
  plot
  ggsave(paste0('/RNAseq2/DEG/down_', exp, ".vs.", ctr, '_GO(BP).pdf'), plot, width=8,height=10)
  
  print("Up")
  up <- subset(upordown, change=="Up")
  ids=bitr(up$X, 'SYMBOL','ENTREZID','org.Hs.eg.db')
  xx2 = enrichGO(ids$ENTREZID,
                 OrgDb = "org.Hs.eg.db",   # Annotation database must correspond to species origin
                 keyType = 'ENTREZID',      
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = T)
  xx2_df <- xx2@result
  write.csv(xx2_df,paste0('/RNAseq2/DEG/up_', exp, ".vs.", ctr, '_GO(BP).csv'),row.names = TRUE)
  # significant_threshold <- 0.05
  significant_threshold <- 1
  significant_up  <- xx2_df[xx2_df$p.adjust < significant_threshold, ]
  # Sort the significant data by p.adjust (or another metric of choice)
  significant_up  <- significant_up[order(significant_up$p.adjust), ]
  # Choose how many top terms to display
  top_n <- 10
  top_terms_up  <- head(significant_up , top_n)
  top_terms_up$neg_log_pvalue <- -log(top_terms_up$pvalue, base = 10)
  top_terms_up$Description <- reorder(top_terms_up$Description, top_terms_up$neg_log_pvalue)
  top_terms_up$Description <- sapply(top_terms_up$Description, function(x) str_wrap(x, width = 30)) 
  top_terms_up$Description <- reorder(top_terms_up$Description, top_terms_up$neg_log_pvalue)
  # Plotting with ggplot2
  plot <- ggplot(top_terms_up, aes(y = Description, x = neg_log_pvalue)) +
    geom_bar(stat = "identity", fill = '#CD2626') + # wine red
    #  coord_flip() +  # Flip the coordinates to make it a horizontal bar plot
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),  # White background
      axis.line = element_line(color = "black"),  # Solid black axis lines
      axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),  # Adjust text angle if necessary
      axis.text.y = element_text(size = 16, colour = "black"),  # Increase size of y-axis labels
      axis.ticks = element_line(color = "black")  # Black axis ticks
    ) +
    labs(x = "-log10(p-value)", y = "", title = "Significant GO Terms for Up-regulated Genes \n (org.Hs.eg.db)") +
    scale_x_continuous(expand = c(0, 0))
  plot
  ggsave(paste0('/RNAseq2/DEG/up_', exp, ".vs.", ctr, '_GO(BP).pdf'), plot, width=8,height=10)
  
}


################### PCA #####################

col_data$group_list <- factor(col_data$group_list)
dds <- DESeqDataSetFromMatrix(countData = counts_data[,c(1, which(colnames(counts_data) %in% col_data$sampleName))],
                              colData = col_data,
                              design = ~ group_list, tidy = T)

# Run DESeq to normalize the data

dds <- DESeq(dds)

# Perform PCA

vsd <- vst(dds, blind = FALSE)  # Variance stabilizing transformation
pcaData <- plotPCA(vsd, intgroup = "group_list", returnData = TRUE)

# Plot the PCA

ggplot(pcaData, aes(x = PC1, y = PC2, color = group_list, shape = group_list)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot of RNA-seq Data",
       x = "PC1",
       y = "PC2") +
  theme(legend.position = "right")

p <- ggplot(data=pcaData,aes(x=PC1,y=PC2,color=group_list,shape=group_list, label=name))+
  geom_point(size=3)+
  theme_bw()+theme(panel.grid=element_blank())+
  geom_text_repel() +
  xlab("PC1")+
  ylab("PC2")+
  theme(legend.position = c(0.9,0.9) + stat_ellipse(level = 0.9))

ggsave("pca_plot.png", width = 10, height = 8, units = "in")

# Plot the dendrogram

correlation_matrix <- cor(counts(dds, normalized =T), method = "spearman")
distance_matrix <- as.dist(1 - correlation_matrix)
hc <- hclust(distance_matrix, method = "average")
# plot(hc, main = "Hierarchical Clustering of NPC Samples", xlab = "", sub = "")
p <- pheatmap(correlation_matrix,
              cluster_rows = T,
              cluster_cols = T,
              clustering_distance_rows = distance_matrix,
              clustering_distance_cols = distance_matrix,
              main = "Heatmap of NPC Samples Correlation",
              color = colorRampPalette(c("white", "violet"))(30))
