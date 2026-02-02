
# ============================================================
# Cell type annotation workflow (control subset)
# Order and logic preserved; only added publication style comments
# Note: primary annotation was performed on integrated_snn_res.0.3
#       BEFORE subclustering cluster 8.
# ============================================================

## ============================================================
## Neighbors (dims 1:20), clustering (res 0.6), and UMAP (dims 1:20)
## ============================================================
s = readRDS("engraft_human_integrated_Annotation_5Genotypes_noKI.rds")
s <- FindNeighbors(s, dims = 1:20)
s.fn20 <- FindClusters(
  s,
  resolution   = 0.6,
  cluster.name = "new_cluster_res0.6",
  graph.name   = "integrated_snn"
)
s.fn20.umap20 <- RunUMAP(s.fn20, dims = 1:20)

saveRDS(s.fn20.umap20, file = "xeno.rds")


s.fn20.umap20<-xeno
xeno_wt <- subset(xeno, Condition %in% c("3R_WT_Ctrl", "4R_WT_Ctrl"))


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(tidyr)
  library(tibble)
  library(scales)
})


DefaultAssay(xeno_wt) <- "integrated"   # clustering performed on integrated assay

# Recompute PCA and SNN graph (may be missing/outdated after subsetting)
xeno_wt <- RunPCA(xeno_wt, verbose = FALSE)
xeno_wt <- FindNeighbors(xeno_wt, reduction = "pca", dims = 1:30)

# Graph name used by FindClusters (typically "integrated_snn")
graph_name <- paste0(DefaultAssay(xeno_wt), "_snn")

# Run clustering at multiple resolutions and store results in metadata
for (res in c(0.2, 0.3, 0.4, 0.5)) {
  xeno_wt <- FindClusters(
    xeno_wt,
    graph.name = graph_name,
    resolution = res,
    verbose    = FALSE
  )
}

# Confirm clustering columns created (e.g., integrated_snn_res.0.2/0.3/0.4/0.5)
grep(paste0(graph_name, "_res"), colnames(xeno_wt@meta.data), value = TRUE)

# Set res = 0.3 clusters as active identities for primary annotation
Idents(xeno_wt) <- paste0(graph_name, "_res.0.3")

# Save the object containing multi resolution cluster labels
saveRDS(xeno_wt, file = "xeno_wt_res_multi.rds")


Idents(xeno_wt) <- "integrated_snn_res.0.3"
DimPlot(xeno_wt)

# Confirm the clustering column exists (sanity check)
stopifnot("integrated_snn_res.0.3" %in% colnames(xeno_wt@meta.data))

# Use res = 0.3 cluster labels as identities for DEG discovery
Idents(xeno_wt) <- xeno_wt$integrated_snn_res.0.3

# Identify positive markers per cluster using Wilcoxon test
all_markers <- FindAllMarkers(
  xeno_wt,
  only.pos        = TRUE,
  min.pct         = 0.4,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

# Rank markers by detection specificity (pct.diff) then effect size (avg_log2FC)
all_markers <- all_markers %>%
  mutate(pct.diff = pct.1 - pct.2) %>%
  relocate(cluster, gene, avg_log2FC, pct.1, pct.2, pct.diff, p_val, p_val_adj)

# Select top 100 markers per cluster (pct.diff first, then avg_log2FC)
top100_per_cluster <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(pct.diff), desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 100) %>%
  ungroup()

# Export cluster marker list used for manual annotation
write.csv(top100_per_cluster,
          "xeno_wt_top100_markers_pctdiff_first.csv",
          row.names = FALSE)

# ------------------------------------------------------------
## Cluster 8 subclustering (performed AFTER initial cluster annotation)
# Purpose: resolve heterogeneity within cluster 8 using FindSubCluster
# ------------------------------------------------------------

# Make sure identities are set to the res = 0.3 clusters
Idents(xeno_wt) <- "integrated_snn_res.0.3"

# Confirm cluster "8" exists in the res = 0.3 clustering
table(xeno_wt$integrated_snn_res.0.3)

# Select an available SNN graph; prefer "integrated_snn" if present
available_graphs <- names(xeno_wt@graphs)
snn_graphs <- grep("_snn$", available_graphs, value = TRUE)

if ("integrated_snn" %in% snn_graphs) {
  graph_to_use <- "integrated_snn"
} else if (length(snn_graphs) > 0) {
  graph_to_use <- snn_graphs[1]
} else {
  # If no SNN graph exists, rebuild neighbors and create "integrated_snn"
  xeno_wt <- FindNeighbors(
    xeno_wt,
    reduction  = "pca",
    dims       = 1:30,
    graph.name = "integrated_snn"
  )
  graph_to_use <- "integrated_snn"
}

# Subcluster cluster 8 only; results stored in meta.data column "C8_subcluster"
xeno_wt <- FindSubCluster(
  xeno_wt,
  cluster         = 8,
  graph.name      = graph_to_use,
  subcluster.name = "C8_subcluster",
  resolution      = 0.2,
  algorithm       = 1
)

# Inspect subcluster labels (cluster 8 cells will have labels like 8_0, 8_1, 8_2)
table(xeno_wt$C8_subcluster, useNA = "ifany")

# ------------------------------------------------------------
# Order C8_subcluster factor levels for consistent label ordering
# Purpose: cosmetic ordering only; does not change assignment
# ------------------------------------------------------------
xeno_wt$C8_subcluster

labs <- as.character(xeno_wt$C8_subcluster)
u    <- unique(labs)

main <- as.integer(sub("_.*", "", u))                              # main cluster number
subn <- ifelse(grepl("_", u), as.integer(sub(".*_", "", u)), Inf)  # subcluster index for 8_*
ord  <- order(main, subn)

xeno_wt$C8_subcluster <- factor(labs, levels = u[ord])

# ------------------------------------------------------------
# Merge cluster 8 subcluster labels into a combined identity label
# Purpose: keep all clusters as 0/1/2/... while replacing cluster 8 with 8_0/8_1/8_2
# ------------------------------------------------------------
base_cluster_col <- "integrated_snn_res.0.3"
sub_col          <- "C8_subcluster"

stopifnot(base_cluster_col %in% colnames(xeno_wt@meta.data))
stopifnot(sub_col %in% colnames(xeno_wt@meta.data))

base_lab <- as.character(xeno_wt[[base_cluster_col]][, 1])
sub_lab  <- as.character(xeno_wt[[sub_col]][, 1])
is8      <- base_lab == "8" | base_lab == 8

full_lab <- base_lab
full_lab[is8 & !is.na(sub_lab) & grepl("^8_", sub_lab)] <-
  sub_lab[is8 & !is.na(sub_lab) & grepl("^8_", sub_lab)]
full_lab <- as.character(full_lab)

# Preserve original C8_subcluster column, then overwrite with combined labels
if (!"C8_subcluster_orig" %in% colnames(xeno_wt@meta.data)) {
  xeno_wt$C8_subcluster_orig <- xeno_wt[[sub_col]][, 1]
}
xeno_wt[[sub_col]] <- full_lab

# Use combined labels as identities for downstream marker discovery
Idents(xeno_wt) <- "C8_subcluster"
DefaultAssay(xeno_wt) <- "RNA"

# ------------------------------------------------------------
# H) Marker finding across combined identities (includes 8_0/8_1/8_2)
# Purpose: obtain marker genes for refined identities after subclustering cluster 8
# ------------------------------------------------------------
all_markers_sub <- FindAllMarkers(
  object          = xeno_wt,
  only.pos        = TRUE,
  min.pct         = 0.4,
  logfc.threshold = 0.25,
  test.use        = "wilcox",
  verbose         = TRUE
)

# Rank markers by pct_diff then avg_log2FC within each identity
all_markers_ranked <- all_markers_sub %>%
  dplyr::mutate(pct_diff = pct.1 - pct.2) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(dplyr::desc(pct_diff),
                 dplyr::desc(avg_log2FC),
                 p_val_adj,
                 .by_group = TRUE) %>%
  dplyr::ungroup()

write.csv(all_markers_ranked,
          "all_clusters_markers_ranked_by_pctdiff.csv",
          row.names = FALSE)

# Safety: ensure expected columns exist before formatting/ranking
stopifnot(all(c("cluster","gene","pct.1","pct.2","avg_log2FC","p_val","p_val_adj") %in%
                names(all_markers_sub)))

# Add pct.diff and tidy columns for export and reporting
all_markers_sub <- all_markers_sub %>%
  mutate(
    pct.1      = as.numeric(pct.1),
    pct.2      = as.numeric(pct.2),
    avg_log2FC = as.numeric(avg_log2FC),
    pct.diff   = pct.1 - pct.2
  ) %>%
  relocate(cluster, gene, avg_log2FC, pct.1, pct.2, pct.diff, p_val, p_val_adj)

# Select top 100 markers per identity (pct.1 then pct.diff then avg_log2FC)
top100_per_cluster_sub <- all_markers_sub %>%
  group_by(cluster) %>%
  dplyr::arrange(dplyr::desc(pct.1),
                 dplyr::desc(pct.diff),
                 dplyr::desc(avg_log2FC),
                 p_val_adj,
                 .by_group = TRUE) %>%
  slice_head(n = 100) %>%
  ungroup()


write.csv(top100_per_cluster_sub,
          "sub8_top100_by_pct1_then_pctdiff.csv",
          row.names = FALSE)



## ------------------------------------------------------------
#Top expressed genes within a specific subcluster (example: 8_0)
## ------------------------------------------------------------

stopifnot("C8_subcluster" %in% colnames(xeno_wt@meta.data))

assay_to_use <- if ("RNA" %in% Assays(xeno_wt)) "RNA" else DefaultAssay(xeno_wt)

cells_8_0 <- WhichCells(xeno_wt, expression = C8_subcluster == "8_0")
if (length(cells_8_0) == 0) {
  stop('No cells labeled "8_0". Available labels: ',
       paste(sort(unique(xeno_wt$C8_subcluster)), collapse = ", "))
}

mat <- GetAssayData(xeno_wt, assay = assay_to_use, layer  = "data")
avg_expr_vec <- Matrix::rowMeans(mat[, cells_8_0, drop = FALSE])

top100_avg <- tibble(
  gene     = names(avg_expr_vec),
  avg_expr = as.numeric(avg_expr_vec)
) %>%
  arrange(desc(avg_expr)) %>%
  slice_head(n = 100)

top100_avg_no_mito_ribo <- top100_avg %>%
  filter(!grepl("^MT-", gene, ignore.case = TRUE),
         !grepl("^RPL|^RPS", gene)) %>%
  slice_head(n = 100)

write.csv(top100_avg_no_mito_ribo, "C8_0_top100_by_avgExpr_noMitoRibo.csv",  row.names = FALSE)

#Based on the gene expression profile of each cluster, cell type annotation was given to each cluster



# 1) Map: cluster id -> cell type
map_ctrl <- c(
  "0"   = "AST",
  "1"   = "MG",
  "2"   = "Neu_IN_BCL11B",
  "3"   = "Neu_IN_Gly",
  "4"   = "Neu_IN_BCL11B",
  "5"   = "Neu_EN_RELN",
  "6"   = "AST",
  "7"   = "Neu_IN_Gly",
  "8_0" = "Neu_Early",
  "8_1" = "AST",
  "8_2" = "VLMC",
  "9"   = "Neu_IN_PROX1",
  "10"  = "Neu_EN_PIEZO2",
  "11"  = "Neu_EN_KCNC2",
  "12"  = "Neu_IN_Gly",
  "13"  = "Neu_IN_Gly",
  "14"  = "Neu_IN_RELN",
  "15"  = "Neu_EN_BCL11B",
  "16"  = "Neu_IN_PROX1",
  "17"  = "Neu_IN_Gly",
  "18"  = "Neu_EN_BCL11B",
  "19"  = "Neu_EN_PROX1",
  "20"  = "Neu_EN_KCNC2",
  "21"  = "MG",
  "22"  = "Neu_EN_PIEZO2",
  "23"  = "NPC_RG"
)

# 2) Apply mapping to metadata
key <- as.character(xeno_wt$C8_subcluster)
xeno_wt$control_celltype <- unname(map_ctrl[key])

# 3) Optional: make it an ordered factor and sanity check
xeno_wt$control_celltype <- factor(
  xeno_wt$control_celltype,
  levels = unique(unname(map_ctrl))
)

# Report any cluster labels without a mapping
unmapped <- setdiff(unique(key), names(map_ctrl))
if (length(unmapped) > 0) {
  message("Unmapped cluster labels: ", paste(unmapped, collapse = ", "))
}

# Quick check
table(xeno_wt$control_celltype, useNA = "ifany")

colnames(xeno_wt@meta.data)
View(xeno_wt)




###############################################################################
## Figure 1: Panel-by-panel code organization
## Required inputs:
##   - xeno_wt  (Seurat object)
##   - map_ctrl (named vector: names = C8_subcluster labels, values = control_celltype)
###############################################################################

library(Seurat)
library(ggplot2)

###############################################################################
## Panel f (and g): UMAP colored by control_celltype with % in legend
###############################################################################

## (f,g-1) Apply mapping to metadata
key <- as.character(xeno_wt$C8_subcluster)
xeno_wt$control_celltype <- unname(map_ctrl[key])

## sanity: unmapped clusters
unmapped <- setdiff(unique(key), names(map_ctrl))
if (length(unmapped) > 0) message("Unmapped cluster labels: ", paste(unmapped, collapse = ", "))

## (f,g-2) Set desired base order for control_celltype
## ----- Set desired_base levels for control_celltype -----
ct_all <- unique(as.character(xeno_wt$control_celltype))
ct_all <- ct_all[!is.na(ct_all)]

EN_levels  <- sort(ct_all[grepl("^Neu_EN_", ct_all)])
IN_levels  <- sort(ct_all[grepl("^Neu_IN_", ct_all)])
NE_levels  <- intersect("Neu_Early", ct_all)
AST_levels <- intersect("AST", ct_all)
MG_levels  <- intersect("MG", ct_all)
NPC_levels <- intersect("NPC_RG", ct_all)

## anything not captured above (eg VLMC) goes last
desired_base <- c(EN_levels, IN_levels, NE_levels, AST_levels, MG_levels, NPC_levels)
desired_base <- c(desired_base, setdiff(ct_all, desired_base))  # leftovers at end

xeno_wt$control_celltype <- factor(xeno_wt$control_celltype, levels = desired_base)


## ----- Make a 3R WT control object -----
xeno_3R <- subset(xeno_wt, subset = Condition == "3R_WT_Ctrl")


## ----- Build % labels, KEEP desired_base order -----
ct  <- as.character(xeno_3R$control_celltype)
pct <- prop.table(table(ct)) * 100

## IMPORTANT: base_present must follow desired_base order
base_present <- desired_base[desired_base %in% names(pct)]

lab_map <- setNames(
  sprintf("%s (%.1f%%)", base_present, as.numeric(pct[base_present])),
  base_present
)

vals <- unname(lab_map[ct])
names(vals) <- colnames(xeno_3R)

xeno_3R <- AddMetaData(xeno_3R, metadata = vals, col.name = "control_celltype_pct")

## IMPORTANT: levels must be the label strings (not desired_base)
xeno_3R$control_celltype_pct <- factor(
  xeno_3R$control_celltype_pct,
  levels = unname(lab_map)
)

levels(xeno_3R$control_celltype_pct)


## ----- Final plot for Panel f -----
p_f <- DimPlot(
  xeno_3R,
  group.by = "control_celltype_pct",
  label    = TRUE,
  repel    = TRUE,
  pt.size  = 1.2
) + NoAxes()

print(p_f)

#UMAP was exported out for minor editing in illustrator
#percentage of each cell type was exported out to be plotted in Prism



###############################################################################
## Panel h: NPC lineage pseudotime (non-MG cells) using velocity_pseudotime
###############################################################################

## subset: all cell types except MG, within 3R_WT_Ctrl
all_other_celltypes <- subset(
  xeno_wt,
  subset = Condition == "3R_WT_Ctrl" & control_celltype != "MG"
)

## pseudotime FeaturePlot (Panel h)
p_h <- FeaturePlot(
  all_other_celltypes,
  features  = "velocity_pseudotime",
  reduction = "umap",
  pt.size   = 1.1,
  raster    = FALSE
) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
  guides(color = guide_colorbar(title = "")) +
  NoAxes() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right"
  )

print(p_h)



###############################################################################
## Panel i: HPC lineage pseudotime (MG only) using latent_time
###############################################################################

## subset: MG only, within 3R_WT_Ctrl
MG_3R_only <- subset(
  xeno_wt,
  subset = control_celltype == "MG" & Condition == "3R_WT_Ctrl"
)

## latent time FeaturePlot (Panel i)
p_i <- FeaturePlot(
  MG_3R_only,
  features   = "latent_time",
  reduction  = "umap",
  pt.size    = 2.2,
  alpha      = 0.8,
  raster     = FALSE
) +
  scale_color_gradient(low = "blue", high = "#FFB90F") +
  guides(color = guide_colorbar(title = "")) +
  NoAxes() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right"
  )

print(p_i)



###############################################################################
## Panels j, k, l: Microglia marker expression (CSF1R, P2RY12, SALL1) in MG UMAP
###############################################################################

############################################################
# Panels j–l: MG marker FeaturePlots (match original style)
############################################################

DefaultAssay(xeno_wt) <- "RNA"

## MG subset (3R WT Ctrl)
MG_3R_only <- subset(
  xeno_wt,
  subset = control_celltype == "MG" & Condition == "3R_WT_Ctrl"
)

## Panel j: CSF1R
p_j <- FeaturePlot(
  MG_3R_only,
  features  = "CSF1R",
  reduction = "umap",
  pt.size   = 2.2,
  order     = TRUE
) +
  scale_color_gradient(low = "lightgrey", high = "red", na.value = "grey90") +
  guides(color = guide_colorbar(title = "")) +
  NoAxes() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right"
  )
print(p_j)


## Panel k: P2RY12
p_k <- FeaturePlot(
  MG_3R_only,
  features  = "P2RY12",
  reduction = "umap",
  pt.size   = 2.2,
  order     = TRUE
) +
  scale_color_gradient(low = "lightgrey", high = "red", na.value = "grey90") +
  guides(color = guide_colorbar(title = "")) +
  NoAxes() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right"
  )
print(p_k)


## Panel l: SALL1
p_l <- FeaturePlot(
  MG_3R_only,
  features  = "SALL1",
  reduction = "umap",
  pt.size   = 2.2,
  order     = TRUE
) +
  scale_color_gradient(low = "lightgrey", high = "red", na.value = "grey90") +
  guides(color = guide_colorbar(title = "")) +
  NoAxes() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right"
  )
print(p_l)


## Panal m: microglia signature comparison
library(Seurat)
library(ggplot2)
library(ggrepel)

setwd("/Users/lifan/Desktop/data_analysis/LG109_engraft/public_datasets/microglia profile comparison")

library(Seurat)

# ===============================
# Step 0: Load Seurat objects
# ===============================

obj1 <- readRDS("/Users/lifan/Desktop/data_analysis/Human_R47H_new/AD_Mayo_UPenn_integration/all_54/MG_only/all_54_MG.rds")
obj1 <- UpdateSeuratObject(obj1)

# extract non-AD samples only
Idents(obj1) <- "orig.ident"
obj1 <- subset(
  obj1,
  idents = c(
    "Non_WT_E2E3_F_1",
    "Non_WT_E2E3_F_2",
    "Non_WT_E2E3_M_1",
    "Non_WT_E2E3_M_2",
    "Non_WT_E2E4_M",
    "Non_WT_E3E3_M_2",
    "Non_WT_NA_F"
  )
)

obj2 <- readRDS("WQ_graft_microglia.rds")

# ===============================
# Step 1: Extract scaled RNA data and aggregate per sample
# ===============================
extract_sample_means <- function(obj) {
  expr <- GetAssayData(obj, assay = "RNA", layer = "scale.data")
  orig_id <- obj$orig.ident[colnames(expr)]
  
  sapply(unique(orig_id), function(s) {
    rowMeans(expr[, orig_id == s, drop = FALSE])
  }, simplify = "matrix")
}

sample_means1 <- extract_sample_means(obj1)
sample_means2 <- extract_sample_means(obj2)

# ===============================
# Step 2: Align genes across objects
# ===============================
common_genes <- intersect(
  rownames(sample_means1),
  rownames(sample_means2)
)

sample_means1 <- sample_means1[common_genes, , drop = FALSE]
sample_means2 <- sample_means2[common_genes, , drop = FALSE]

# ===============================
# Step 3: Combine scRNA-seq samples
# ===============================
sc_samples_all <- cbind(sample_means1, sample_means2)

dataset_labels <- c(
  rep("human_microglia (n=7)", ncol(sample_means1)),
  rep("WQ_engraft (n=3)", ncol(sample_means2))
)

# ===============================
# Step 4: PCA on samples
# ===============================
pca_res <- prcomp(t(sc_samples_all), scale. = FALSE)
pca_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

pca_df <- data.frame(
  Sample = colnames(sc_samples_all),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Dataset = dataset_labels
)


# ===============================
# Step 5: Add bulk RNA-seq samples
# ===============================

Ginroux <- read.csv("Ginroux_GSE241127_readcounts_genes.csv")
rownames(Ginroux) <- Ginroux$X
Ginroux$X <- NULL

# Ensure numeric matrix
bulk_counts <- apply(Ginroux, 2, as.numeric)
rownames(bulk_counts) <- rownames(Ginroux)

# CPM normalization + log1p
bulk_cpm <- t(t(bulk_counts) / colSums(bulk_counts)) * 1e6
bulk_log <- log1p(bulk_cpm)

# ===============================
# Project bulk into scRNA PCA
# ===============================

# Genes used by PCA (IN ORDER)
pca_genes_used <- rownames(pca_res$rotation)

# Initialize bulk matrix with all PCA genes
bulk_full <- matrix(
  0,
  nrow = length(pca_genes_used),
  ncol = ncol(bulk_log),
  dimnames = list(pca_genes_used, colnames(bulk_log))
)

# Fill overlapping genes
shared_genes <- intersect(pca_genes_used, rownames(bulk_log))
bulk_full[shared_genes, ] <- bulk_log[shared_genes, ]

# Transpose for prediction
bulk_newdata <- t(bulk_full)

# FINAL sanity check (must be TRUE)
stopifnot(all(colnames(bulk_newdata) == rownames(pca_res$rotation)))

# Project bulk samples
bulk_pca_proj <- predict(pca_res, newdata = bulk_newdata)

bulk_df <- data.frame(
  Sample  = colnames(bulk_log),
  PC1     = bulk_pca_proj[, 1],
  PC2     = bulk_pca_proj[, 2],
  Dataset = "Ginroux_bulk (n=6)"
)
#
# ===============================
# Step 6: Combine scRNA-seq + bulk
# ===============================
plot_df <- rbind(pca_df, bulk_df)

# ===============================
# Step 7: Plot PCA
# ===============================

library(ggplot2)

ggplot(plot_df, aes(x = PC1, y = PC2, color = Dataset, shape = Dataset)) +
  geom_point(size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Sample-level PCA",
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)")
  ) +
  scale_color_manual(values = c(
    "human_microglia (n=7)" = "#1f78b4",
    "WQ_engraft (n=3)"      = "#33a02c",
    "Ginroux_bulk (n=6)"   = "#ff7f00"
  )) +
  scale_shape_manual(values = c(
    "human_microglia (n=7)" = 16,
    "WQ_engraft (n=3)"      = 17,
    "Ginroux_bulk (n=6)"   = 8
  ))

ggsave(
  "microglia_profile_comparisons_PCA.pdf",
  plot = last_plot(),
  width = 7, height = 4.5, dpi = 600, limitsize = FALSE
)

ggsave(
  "microglia_profile_comparisons_PCA_large.pdf",
  plot = last_plot(),
  width = 12, height = 10, dpi = 600, limitsize = FALSE
)

ggsave(
  "microglia_profile_comparisons_PCA_figure.pdf",
  plot = last_plot(),
  width = 5, height = 3, dpi = 600, limitsize = FALSE
)


