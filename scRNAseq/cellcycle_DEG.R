library(dplyr)
library(magrittr)
library(patchwork)
library(Seurat)
library(presto)
library(Cairo)
library(ggplot2)
library(ggsci)
setwd('data/merge_data/scRNA/cluster_0.5')
C_WT <- readRDS("H9_ESC_cellcycle.rds")
C_24h <- readRDS("NPC_24h_cellcycle.rds")
C_44h <- readRDS("NPC_44h_cellcycle.rds")
C_46h <- readRDS("NPC_46h_cellcycle.rds")
C_48h <- readRDS("NPC_48h_cellcycle.rds")
C_50h <- readRDS("NPC_50h_cellcycle.rds")
C_52h <- readRDS("NPC_52h_cellcycle.rds")
C_54h <- readRDS("NPC_54h_cellcycle.rds")
C_72h <- readRDS("NPC_72h_cellcycle.rds")
C_d6 <- readRDS("NPC_D6_cellcycle.rds")

samples <- c(C_WT, C_24h, C_44h, C_46h, C_48h, C_50h, C_52h, C_54h, C_72h, C_d6)
early_samples <- c(C_WT, C_24h, C_44h, C_46h, C_48h, C_50h, C_52h, C_54h)
samples_name <- c('C_WT', 'C_24h', 'C_44h', 'C_46h', 'C_48h', 'C_50h', 'C_52h', 'C_54h', 'C_72h', 'C_d6')
early_samples_name <- c('C_WT', 'C_24h', 'C_44h', 'C_46h', 'C_48h', 'C_50h', 'C_52h', 'C_54h')
for (i in 1:10) {
  samples[[i]] <- RenameCells(samples[[i]], add.cell.id = samples_name[i])
}
names(samples) <- samples_name

scRNA <- merge(samples[[1]], samples[2:length(samples)])

scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
scRNA <- RunPCA(scRNA, verbose = F)
ElbowPlot(scRNA, ndims = 50)
table(scRNA$orig.ident)
scRNA <- FindNeighbors(scRNA, dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.5)
## UMAP
scRNA <- RunUMAP(scRNA, dims = 1:30)
p <- DimPlot(scRNA, reduction = "umap", label = T)

# Here set orig.idents as active.idents
Idents(scRNA)="orig.ident"
scRNA <- RenameIdents(scRNA, "WT"="ESC", "24h"="NPC_24h", 
                      '72h'="NPC_72h", 'D6'="NPC_D6")
# Must give Idents(scRNA) back to orig.ident, then Rename could be available
scRNA$orig.ident = Idents(scRNA)

p <- DimPlot(scRNA, group.by = "orig.ident", label = T, label.size=3)
ggsave("1.1_ident_cluster.tiff",p, width=10,height=8)

p <- FeaturePlot(scRNA,features = c('NANOG', 'POU5F1', 'PAX6', 'SOX2', 'VIM', 'NES'), reduction = "umap", label = T, ncol = 3) # ncol for amounts of fig in a row
ggsave("1.2_features.tiff",p, width = 10,height = 8)

Idents(scRNA)="orig.ident"
scRNA <- RenameIdents(scRNA, "44h"="NPC_44h", "46h"="NPC_46h", "48h"="NPC_48h", "50h"="NPC_50h",
                      "52h"="NPC_52h", "54h"="NPC_54h")
scRNA$orig.ident = Idents(scRNA)
levels(scRNA) <- c('ESC', "NPC_24h", "NPC_44h", "NPC_46h", "NPC_48h", "NPC_50h", "NPC_52h", "NPC_54h", "NPC_72h", "NPC_D6")

p <- DotPlot(scRNA, features = unique(c('NANOG', 'POU5F1', 'PAX6', 'SOX2', 'VIM', 'NES')), group.by = "orig.ident")+RotatedAxis()+scale_x_discrete("")+scale_y_discrete("")
ggsave("1.3_markers.tiff",p, width=8,height=6)


scRNA <- CellCycleScoring(scRNA, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = T)
scRNA <- RunPCA(scRNA, features = c(cc.genes$s.genes, cc.genes$g2m.genes))
levels(scRNA)

metadata <- scRNA@meta.data
# directly concatenate orig.ident and Phase to produce a new column as cell type
metadata$cell_type<- factor(paste0(metadata$orig.ident, '_', metadata$Phase))
scRNA@meta.data <- metadata 
Idents(scRNA) = 'cell_type'
levels(scRNA) <- c('ESC_G1','ESC_S','ESC_G2M',
                   'NPC_24h_G1','NPC_24h_S','NPC_24h_G2M',
                   'NPC_44h_G1','NPC_44h_S','NPC_44h_G2M',
                   'NPC_46h_G1','NPC_46h_S','NPC_46h_G2M',
                   'NPC_48h_G1','NPC_48h_S','NPC_48h_G2M',
                   'NPC_50h_G1','NPC_50h_S','NPC_50h_G2M',
                   'NPC_52h_G1','NPC_52h_S','NPC_52h_G2M',
                   'NPC_54h_G1','NPC_54h_S','NPC_54h_G2M',
                   'NPC_72h_G1','NPC_72h_S','NPC_72h_G2M',
                   'NPC_D6_G1','NPC_D6_S','NPC_D6_G2M')

### Differential gene expression
cluster1.markers <- FindMarkers(scRNA, ident.1 = "NPC_D6_G2M", ident.2 = "NPC_72h_G2M", min.pct = 0.1)  # adjust with needed groups
log2FC = 2
padj = 0.05
cluster1.markers$threshold="ns";
cluster1.markers[which(cluster1.markers$avg_log2FC  > log2FC & cluster1.markers$p_val_adj < padj),]$threshold="up";
cluster1.markers[which(cluster1.markers$avg_log2FC  < (-log2FC) & cluster1.markers$p_val_adj < padj),]$threshold="down";
cluster1.markers$threshold=factor(cluster1.markers$threshold, levels=c('down','ns','up'))
cluster1.markers$Label = ""
upordown <- subset(cluster1.markers, threshold=="up"|threshold=='down')   
upordown.gene <- c(as.character(rownames(upordown)))
cluster1.markers$Label[match(upordown.gene, rownames(cluster1.markers))] <- upordown.gene
upordown$Label[match(upordown.gene, rownames(upordown))] <- upordown.gene
write.csv(upordown,'2.1_D6.vs.D3_G2M.csv',row.names = TRUE)

color <- c("blue", "grey","red")
lb <- c()
  for (j in c(1:length(table(cluster1.markers$threshold)))) {
    lb_num <- paste0(names(table(cluster1.markers$threshold))[[j]], "(",table(cluster1.markers$threshold)[[j]],')')
    if (grepl("^ns", lb_num)) {
      lb_num <- "ns"
    }
    lb <- c(lb, lb_num)
  }
  lb_color <- color[names(table(cluster1.markers$threshold))]

g <- ggplot(data=cluster1.markers, aes(x=avg_log2FC, y=-log10(p_val_adj), color=threshold)) +
  geom_point(alpha=0.8, size=0.8) +
  geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
  geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
  xlab(bquote(Log[2]*FoldChange))+
  ylab(bquote(-Log[10]*italic(P.adj)))+
  labs(title='D6_G2M vs. D3_G2M')+
  theme_classic(base_size = 14) +
  scale_color_manual('', labels=lb, values=lb_color)+
  guides(color=guide_legend(override.aes = list(size=3, alpha=1)))
ggsave("2.1_D6.vs.D3_G2M_volcano.tiff",g, width=8,height=6)

ids=bitr(upordown$Label, 'SYMBOL','ENTREZID','org.Hs.eg.db')
upordown=merge(upordown,ids,by.x='Label',by.y='SYMBOL')
gcSample=split(upordown$ENTREZID,upordown$threshold)
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",   # Annotation database must correspond to species origin
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01)
p <- dotplot(xx) + scale_x_discrete(drop=FALSE)
ggsave("2.1_D6.vs.D3_G2M_GO.tiff",p, width=5,height=6)
