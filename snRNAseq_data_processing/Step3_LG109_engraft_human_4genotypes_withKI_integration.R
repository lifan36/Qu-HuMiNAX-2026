
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/athena/ganlab/scratch/lif4001/LG109_engraft/DF_2ndRound")

LG109_4_human <- readRDS(file = "LG109_4_human_singlets_PCA.rds")
LG109_5_human <- readRDS(file = "LG109_5_human_singlets_PCA.rds")
LG109_6_human <- readRDS(file = "LG109_6_human_singlets_PCA.rds")
LG109_139_human <- readRDS(file = "LG109_139_human_singlets_PCA.rds")

LG109_13_human <- readRDS(file = "LG109_13_human_singlets_PCA.rds")
LG109_41_human <- readRDS(file = "LG109_41_human_singlets_PCA.rds")
LG109_48_human <- readRDS(file = "LG109_48_human_singlets_PCA.rds")
LG109_75_76_human <- readRDS(file = "LG109_75_76_human_singlets_PCA.rds")
LG109_79_human <- readRDS(file = "LG109_79_human_singlets_PCA.rds")

LG109_71_human <- readRDS(file = "LG109_71_human_singlets_PCA.rds")
LG109_86_human <- readRDS(file = "LG109_86_human_singlets_PCA.rds")
LG109_90_human <- readRDS(file = "LG109_90_human_singlets_PCA.rds")

LG109_95_human <- readRDS(file = "LG109_95_human_singlets_PCA.rds")
LG109_97_human <- readRDS(file = "LG109_97_human_singlets_PCA.rds")
LG109_102_human <- readRDS(file = "LG109_102_human_singlets_PCA.rds")
LG109_103_human <- readRDS(file = "LG109_103_human_singlets_PCA.rds")

setwd("/athena/ganlab/scratch/lif4001/LG109_engraft/integration_human_4genotypes_withKI")

P301S_Ctrl <- c(LG109_4_human,LG109_5_human,LG109_6_human,LG109_139_human)
anchors_P301S_Ctrl <- FindIntegrationAnchors(object.list = P301S_Ctrl, dims = 1:30)
P301S_Ctrl_integrated <- IntegrateData(anchorset = anchors_P301S_Ctrl, dims = 1:30)
rm(LG109_4_human,LG109_5_human,LG109_6_human,LG109_139_human, P301S_Ctrl)

P301S_K18 <- c(LG109_13_human,LG109_41_human,LG109_48_human, LG109_75_76_human,LG109_79_human)
anchors_P301S_K18 <- FindIntegrationAnchors(object.list = P301S_K18, dims = 1:30)
P301S_K18_integrated <- IntegrateData(anchorset = anchors_P301S_K18, dims = 1:30)
rm(LG109_13_human,LG109_41_human,LG109_48_human, LG109_75_76_human,LG109_79_human, P301S_K18)

KI_P301S_K18 <- c(LG109_71_human,LG109_86_human, LG109_90_human)
anchors_KI_P301S_K18 <- FindIntegrationAnchors(object.list = KI_P301S_K18, dims = 1:30)
KI_P301S_K18_integrated <- IntegrateData(anchorset = anchors_KI_P301S_K18, dims = 1:30)
rm(LG109_71_human,LG109_86_human, LG109_90_human, KI_P301S_K18)

ISOKI_P301S_K18 <- c(LG109_95_human,LG109_97_human,LG109_102_human, LG109_103_human)
anchors_ISOKI_P301S_K18 <- FindIntegrationAnchors(object.list = ISOKI_P301S_K18, dims = 1:30)
ISOKI_P301S_K18_integrated <- IntegrateData(anchorset = anchors_ISOKI_P301S_K18, dims = 1:30)
rm(LG109_95_human,LG109_97_human,LG109_102_human, LG109_103_human, ISOKI_P301S_K18)

engraft_human <- c(P301S_Ctrl_integrated, P301S_K18_integrated, KI_P301S_K18_integrated, ISOKI_P301S_K18_integrated)
anchors_engraft_human <- FindIntegrationAnchors(object.list = engraft_human, dims = 1:30)
engraft_human_integrated <- IntegrateData(anchorset = anchors_engraft_human, dims = 1:30)
rm(P301S_Ctrl_integrated, P301S_K18_integrated, KI_P301S_K18_integrated, ISOKI_P301S_K18_integrated, engraft_human)

#saveRDS(engraft_human_integrated, file = "engraft_human_integrated.rds")

DefaultAssay(engraft_human_integrated) <- 'integrated'

# engraft_human_integrated <- NormalizeData(engraft_human_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# engraft_human_integrated <- FindVariableFeatures(engraft_human_integrated, selection.method = "vst", nfeatures = 3000)

engraft_human_integrated <- ScaleData(engraft_human_integrated, verbose = FALSE)
engraft_human_integrated <- RunPCA(engraft_human_integrated, features = VariableFeatures(object = engraft_human_integrated), verbose = FALSE)

engraft_human_integrated <- FindNeighbors(engraft_human_integrated, dims = 1:15)
engraft_human_integrated <- FindClusters(engraft_human_integrated, resolution = 0.1)
engraft_human_integrated <- RunUMAP(engraft_human_integrated, dims = 1: 15)

DefaultAssay(engraft_human_integrated) <- 'RNA'
engraft_human_integrated <- NormalizeData(engraft_human_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
engraft_human_integrated <- ScaleData(engraft_human_integrated, features = rownames(engraft_human_integrated))

engraft_human_integrated <- JoinLayers(engraft_human_integrated)
#saveRDS(engraft_human_integrated, file = 'engraft_human_integrated_PCA_0.1.rds')
#engraft_human_integrated <- readRDS(file = "engraft_human_integrated_PCA_0.1.rds")

engraft_human_integrated$Condition <- factor(x = engraft_human_integrated$Condition, levels = c("4R_P301S_Ctrl","4R_P301S_K18","KI_4R_P301S_K18","ISOKI_4R_P301S_K18"))
engraft_human_integrated$Sample_Name <- factor(x = engraft_human_integrated$Sample_Name, levels = c("4R_P301S_Ctrl_1","4R_P301S_Ctrl_2","4R_P301S_Ctrl_3","4R_P301S_Ctrl_4",
                                                                                                    "4R_P301S_K18_1","4R_P301S_K18_2","4R_P301S_K18_3","4R_P301S_K18_4","4R_P301S_K18_5",
                                                                                                    "KI_4R_P301S_K18_2","KI_4R_P301S_K18_3","KI_4R_P301S_K18_4",
                                                                                                    "ISOKI_4R_P301S_K18_1","ISOKI_4R_P301S_K18_2","ISOKI_4R_P301S_K18_3","ISOKI_4R_P301S_K18_4"))

pdf("engraft_human_QC_further_filter.pdf", width=9, height=4)
Idents(engraft_human_integrated) <- "Condition"
VlnPlot(object = engraft_human_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()


Idents(engraft_human_integrated) <- "Sample_Name"
pdf("engraft_human_QC_Sample_further_filter.pdf", width=18, height=4)

VlnPlot(object = engraft_human_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(engraft_human_integrated) <- "seurat_clusters"
pdf("engraft_human_integrated_umap.pdf", width=5, height=4)
DimPlot(engraft_human_integrated, reduction = 'umap', label = T)
dev.off()
pdf("engraft_human_integrated_umap_split_individual.pdf", width=8, height=8)
DimPlot(engraft_human_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
pdf("engraft_human_integrated_umap_split_Condition.pdf", width=10, height=2.5)
DimPlot(engraft_human_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()

write.csv(table(engraft_human_integrated$seurat_clusters, engraft_human_integrated$Sample_Name), "engraft_human_cell_counts_cluster_by_sample.csv")

DefaultAssay(engraft_human_integrated) <- 'RNA'
#Add marker genes
pdf("engraft_human_annotation_combine.pdf", width=20, height=5)
DotPlot(object = engraft_human_integrated, features = c("SYT1","SNAP25","GRIN1","SLC17A7","SLC17A6", "CAMK2A", "NRGN","GAD1", "GAD2","SEMA3C","LMX1B","PLP1", "MBP", "MOBP","AQP4","GFAP", 
                                                        "CD74","CSF1R","C3","PDGFRA","VCAN")) + RotatedAxis()
dev.off()


engraft_human_markers <- FindAllMarkers(engraft_human_integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1, test.use = "MAST")
write.csv(engraft_human_markers, "engraft_human_markers.csv")


saveRDS(engraft_human_integrated, file = 'engraft_human_integrated_PCA_0.1.rds')

#engraft_human_integrated <- readRDS(file = "engraft_human_integrated_PCA_0.1.rds")
#engraft_human_markers <- read.csv(file = "engraft_human_markers.csv", header=T,row.names =1)
#top5 <- engraft_human_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#top5$gene <- as.character(top5$gene)
#pdf("engraft_human_HeatMapTop5_0.1_new.pdf", width=24, height=16)
#DoHeatmap(engraft_human_integrated, features = top5$gene) + NoLegend()
#dev.off()




