library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)

filenames=c("NPC_24h","NPC_44h","NPC_46h","NPC_48h","NPC_50h","NPC_52h","NPC_54h","NPC_72h","NPC_D6","H9_ESC")

#############scRNA_quality_control###################
for(i in 1:10){
  filename <- filenames[i]
  data <- Read10X(data.dir = paste0("data/merge_data/cellranger_arc/NPC_",filename,"/outs/filtered_feature_bc_matrix/"))
  # Initialize the Seurat object with the raw (non-normalized data).
  scRNA <- CreateSeuratObject(counts = data$`Gene Expression`, project = filename, min.cells = 3, min.features = 200)
  scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
  # Visualize QC metrics as a violin plot
  VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(filename = paste0("result/merge_data/quality_control/",filename,"/raw_feature_count_mt.png"),width = 15,height = 10)
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  # 
  # plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
  # plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  # plot1 + plot2
  # ggsave(filename = paste0("result/merge_data/quality_control/",filename,"/feature_count_and_mt_count.png"),width = 15,height = 10)
  scRNA <- subset(scRNA, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & percent.mt < 20)
  scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
  scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
  scRNA <- ScaleData(scRNA)
  scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
  scRNA <- RunUMAP(scRNA, dims = 1:30)
  scRNA <- FindNeighbors(scRNA, dims = 1:30)
  scRNA <- FindClusters(scRNA, resolution = 0.5)
  DimPlot(scRNA,label = T,repel = T,pt.size = 1.5)
  ggsave(paste0("result/merge_data/cluster/",filename,"_0.5.png"),width = 15,height = 10)
  VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "orig.ident")
  ggsave(filename = paste0("result/merge_data/quality_control/",filename,"/feature_count_and_mt_count_after_filter.png"),width = 15,height = 10)
  saveRDS(scRNA,paste0("data/merge_data/scRNA/cluster_0.5/",filename,".rds"))
}

#######################scRNA-seq merge####################################
npc44h <- readRDS("data/merge_data/scRNA/cluster_0.5/NPC_44h.rds")
npc46h <- readRDS("data/merge_data/scRNA/cluster_0.5/NPC_46h.rds")
npc48h <- readRDS("data/merge_data/scRNA/cluster_0.5/NPC_48h.rds")
npc50h <- readRDS("data/merge_data/scRNA/cluster_0.5/NPC_50h.rds")
npc52h <- readRDS("data/merge_data/scRNA/cluster_0.5/NPC_52h.rds")
npc54h <- readRDS("data/merge_data/scRNA/cluster_0.5/NPC_54h.rds")
npc24h <- readRDS("data/merge_data/scRNA/cluster_0.5/NPC_24h.rds")
npc72h <- readRDS("data/merge_data/scRNA/cluster_0.5/NPC_72h.rds")
WT <- readRDS("data/merge_data/scRNA/cluster_0.5/H9_ESC.rds")
D6 <- readRDS("data/merge_data/scRNA/cluster_0.5/D6.rds")
npc <- merge(npc24h,y=c(npc44h,npc46h,npc48h,npc50h,npc52h,npc54h,npc72h,WT,D6),
             add.cell.ids = c("npc24h","npc44h", "npc46h", "npc48h","npc50h","npc52h","npc54h","npc72h","WT","D6"),
             project = "NPC_RNA")
npc <- NormalizeData(npc, normalization.method = "LogNormalize", scale.factor = 10000)
npc <- FindVariableFeatures(npc, selection.method = "vst", nfeatures = 2000)
npc <- ScaleData(npc)
npc <- RunPCA(npc, features = VariableFeatures(object = npc))
npc <- RunUMAP(npc, dims = 1:30)
npc <- FindNeighbors(npc, dims = 1:30)
npc <- FindClusters(npc, resolution = 0.5)
saveRDS(npc,"data/merge_data/scRNA/cluster_0.5/npc.rds")
DimPlot(npc,label = T,repel = T,pt.size = 1.5, group.by = "orig.ident",raster = F)
ggsave("result/merge_data/annotation/merge_orig.ident_final.png",width = 15,height = 10)
