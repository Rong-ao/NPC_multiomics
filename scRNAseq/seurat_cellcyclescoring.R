library(Seurat)
library(ggplot2)
library(stringr)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
filenames=c("NPC_24h","NPC_44h","NPC_46h","NPC_48h","NPC_50h","NPC_52h","NPC_54h","NPC_72h","NPC_D6","H9_ESC")

#Predict cell cycle for merged samples
npc <- readRDS("data/merge_data/scRNA/cluster_0.5/npc.rds")
npc.cellcycle <- CellCycleScoring(npc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
npc.cellcycle$Phase<-factor(npc.cellcycle$Phase,levels = c("G1","S","G2M"))
DimPlot(npc.cellcycle,label = T,repel = T,pt.size = 1.5,raster = F,group.by = "Phase")
ggsave("result/merge_data/Phase/RNA_merge_phase.png",width = 15,height = 10)
DimPlot(npc.cellcycle,label = T,repel = T,pt.size = 1.5, cells.highlight = rownames(npc.cellcycle@meta.data)[which(npc.cellcycle$Phase=="G1")],cols.highlight="#F8766D",raster = F)
ggsave("result/merge_data/Phase/RNA_merge_phase_G1.png",width = 15,height = 10)
DimPlot(npc.cellcycle,label = T,repel = T,pt.size = 1.5,cells.highlight = rownames(npc.cellcycle@meta.data)[which(npc.cellcycle$Phase=="S")],cols.highlight = "#0cb702",raster = F)
ggsave("result/merge_data/Phase/RNA_merge_phase_S.png",width = 15,height = 10)
DimPlot(npc.cellcycle,label = T,repel = T,pt.size = 1.5,cells.highlight = rownames(npc.cellcycle@meta.data)[which(npc.cellcycle$Phase=="G2M")],cols.highlight="#00a9ff",raster = F)
ggsave("result/merge_data/Phase/RNA_merge_phase_G2M.png",width = 15,height = 10)
#Cell cycle proportion in all samples
npc.cellcycle$orig.ident <- factor(npc.cellcycle$orig.ident,levels = c("WT","24h","44h","46h","48h","50h","52h","54h","72h","D6"))
ggplot(npc.cellcycle@meta.data, aes(x = orig.ident, fill = Phase)) +
  geom_bar(width = 0.3, position = "fill")+
  theme_bw()+
  theme(axis.text.x = element_text(size=20),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("result/merge_data/Phase/merge_phase_percentage.png",width = 15,height = 10)

#Predict cell cycle for each sample
for(i in 1:10){
  filename <- filenames[i]
  scRNA <- readRDS(paste0("data/merge_data/scRNA/cluster_0.5/",filename,".rds"))
  scRNA$Phase<-"S"
  barcode_G1 <- as.data.frame(rownames(npc.cellcycle@meta.data)[which(npc.cellcycle$orig.ident==filename & npc.cellcycle$Phase=="G1")])
  colnames(barcode_G1)[1]<-"barcode"
  barcode_G1$barcode <- str_split_fixed(barcode_G1$barcode, "_", 2)
  barcode_G1 <- barcode_G1$barcode[,2]
  scRNA$Phase[which(rownames(scRNA@meta.data) %in% barcode_G1)]<-"G1"
  barcode_G2M <- as.data.frame(rownames(npc.cellcycle@meta.data)[which(npc.cellcycle$orig.ident==filename & npc.cellcycle$Phase=="G2M")])
  colnames(barcode_G2M)[1]<-"barcode"
  barcode_G2M$barcode <- str_split_fixed(barcode_G2M$barcode, "_", 2)
  barcode_G2M <- barcode_G2M$barcode[,2]
  scRNA$Phase[which(rownames(scRNA@meta.data) %in% barcode_G2M)]<-"G2M"
  scRNA$Phase<-factor(scRNA$Phase,levels = c("G1","S","G2M"))
  DimPlot(scRNA,label = T,repel = T,pt.size = 1.5,group.by = "Phase")
  ggsave(paste0("result/merge_data/Phase/",filename,"_phase.png"),width = 15,height = 10)
  saveRDS(scRNA,paste0("data/merge_data/scRNA/cluster_0.5/",filename,"_cellcycle.rds"))
}
