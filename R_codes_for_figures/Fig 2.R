
# ============================================================
# Figure 2 snRNA-seq analysis and cell type annotation workflow
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


###add volocity 
# Packages
library(Seurat)
library(dplyr)
library(tibble)
library(readr)
xeno_wt <- xeno_wt_res_multi 
# 1) Read the mapping table produced by velocity
#    Expected columns (any subset is OK): cellid.ori, cellid.vc, velocity_pseudotime, latent_time
map_path <- "All_3RWT_4RWT_cellid2velocity_time.tsv"
c2t <- read.delim(map_path, stringsAsFactors = FALSE, check.names = FALSE)

# 2) Build a data frame from xeno_wt meta with a canonical "cell" column = Seurat cell names
meta <- xeno_wt@meta.data %>%
  rownames_to_column(var = "cell")

# 3) Choose the best ID columns to join by
#    Try to match one meta ID column to one mapping ID column by maximizing overlap
meta_id_candidates <- c("cell", "cellID_ori", "cellid.ori", "cellid", "barcode")
map_id_candidates  <- c("cellid.ori", "cellid.vc", "cellid", "barcode", "cell")

meta_id_candidates <- intersect(meta_id_candidates, names(meta))
map_id_candidates  <- intersect(map_id_candidates, names(c2t))

best_match <- NULL
best_overlap <- -1
best_meta_col <- NA_character_
best_map_col  <- NA_character_

for (mcol in meta_id_candidates) {
  mvals <- unique(meta[[mcol]])
  for (pcol in map_id_candidates) {
    pvals <- unique(c2t[[pcol]])
    ov <- length(intersect(mvals, pvals))
    if (ov > best_overlap) {
      best_overlap <- ov
      best_meta_col <- mcol
      best_map_col  <- pcol
    }
  }
}

if (best_overlap <= 0) {
  stop("Could not find any overlapping cell IDs between xeno_wt metadata and the mapping table. 
       Check that the mapping file columns (e.g., 'cellid.ori' or 'cellid.vc') match a column in xeno_wt@meta.data (e.g., 'cellID_ori') or the Seurat cell names.")
}

message(sprintf("Joining by meta column '%s' <-> mapping column '%s' (overlap = %d).",
                best_meta_col, best_map_col, best_overlap))

# 4) Prepare a minimal mapping table with consistent column names
need_cols <- intersect(c("velocity_pseudotime", "latent_time", "cellid.vc", "cellid.ori"), names(c2t))
map_min <- c2t %>%
  select(all_of(c(best_map_col, need_cols))) %>%
  rename(join_id = !!best_map_col)

# 5) Add the matching ID as 'join_id' to meta and left-join
meta_for_join <- meta %>%
  mutate(join_id = .data[[best_meta_col]])

merged <- meta_for_join %>%
  left_join(map_min, by = "join_id")

# 6) Add columns to Seurat metadata (row names must match Seurat cell names)
to_add <- merged %>%
  select(cell,
         velocity_pseudotime = any_of("velocity_pseudotime"),
         latent_time         = any_of("latent_time"),
         sampleid.cellid     = any_of("cellid.vc")) %>%
  column_to_rownames("cell")

# Keep only columns that exist after select
to_add <- to_add[, colnames(to_add)[colSums(!is.na(to_add)) + colSums(is.na(to_add)) > 0], drop = FALSE]

xeno_wt <- AddMetaData(xeno_wt, to_add)

# 7) Quick sanity check: show a few rows
check_cols <- c("cell", best_meta_col, "velocity_pseudotime", "latent_time", "sampleid.cellid")
preview <- xeno_wt@meta.data %>%
  rownames_to_column("cell") %>%
  select(any_of(check_cols)) %>%
  head(10)
print(preview)

# 8) Optional: visualize pseudotime on UMAP
# Make sure your default reduction is UMAP/PCA as needed
if ("velocity_pseudotime" %in% colnames(xeno_wt@meta.data)) {
  p1 <- FeaturePlot(xeno_wt, features = "velocity_pseudotime", reduction = "umap")
  print(p1)
}
if ("latent_time" %in% colnames(xeno_wt@meta.data)) {
  p2 <- FeaturePlot(xeno_wt, features = "latent_time", reduction = "umap")
  print(p2)
}


###############################################################################
## Panel B: NPC lineage pseudotime (non-MG cells) using velocity_pseudotime
###############################################################################

## subset: all cell types except MG, within 4R_WT_Ctrl
all_other_celltypes <- subset(
  xeno_wt,
  subset = Condition == "4R_WT_Ctrl" & control_celltype != "MG"
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
## Panel B: HPC lineage pseudotime (MG only) using latent_time
###############################################################################

## subset: MG only, within 4R_WT_Ctrl
MG_4R_only <- subset(
  xeno_wt,
  subset = control_celltype == "MG" & Condition == "4R_WT_Ctrl"
)

## latent time FeaturePlot (Panel i)
p_i <- FeaturePlot(
  MG_4R_only,
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



#4R 
## ----- Make a 4R WT control object -----
xeno_4R <- subset(xeno_wt, subset = Condition == "4R_WT_Ctrl")


## ----- Build % labels, KEEP desired_base order -----
ct  <- as.character(xeno_4R$control_celltype)
pct <- prop.table(table(ct)) * 100

## IMPORTANT: base_present must follow desired_base order
base_present <- desired_base[desired_base %in% names(pct)]

lab_map <- setNames(
  sprintf("%s (%.1f%%)", base_present, as.numeric(pct[base_present])),
  base_present
)

vals <- unname(lab_map[ct])
names(vals) <- colnames(xeno_4R)

xeno_4R <- AddMetaData(xeno_4R, metadata = vals, col.name = "control_celltype_pct")

## IMPORTANT: levels must be the label strings (not desired_base)
xeno_4R$control_celltype_pct <- factor(
  xeno_4R$control_celltype_pct,
  levels = unname(lab_map)
)

levels(xeno_4R$control_celltype_pct)


## ----- Final plot for Panel c -----
p_f <- DimPlot(
  xeno_4R,
  group.by = "control_celltype_pct",
  label    = TRUE,
  repel    = TRUE,
  pt.size  = 1.2
) + NoAxes()

print(p_f)

#UMAP was exported out for minor editing in Affinity and Illustrator





##separate NPC and VLMC
library(dplyr)
library(tibble)

## collapse detailed cell types into broad cell types
## keep NPC_RG and VLMC as separate categories
xeno_4R$broad_celltype <- case_when(
  grepl("^Neu", as.character(xeno_4R$control_celltype)) ~ "Neu",
  as.character(xeno_4R$control_celltype) == "AST" ~ "AST",
  as.character(xeno_4R$control_celltype) == "MG" ~ "MG",
  as.character(xeno_4R$control_celltype) == "NPC_RG" ~ "NPC_RG",
  as.character(xeno_4R$control_celltype) == "VLMC" ~ "VLMC",
  TRUE ~ NA_character_
)

## use character vectors to avoid unused factor levels
sample_vec <- as.character(xeno_4R$Sample_Name)
broad_vec  <- as.character(xeno_4R$broad_celltype)

## remove NA if any
keep <- !is.na(sample_vec) & !is.na(broad_vec)

sample_broad_counts <- table(
  sample_vec[keep],
  broad_vec[keep]
)

sample_broad_pct <- prop.table(sample_broad_counts, margin = 1) * 100

## force clean column order
celltype_order <- c("Neu", "AST", "MG", "NPC_RG", "VLMC")

sample_broad_pct_df <- as.data.frame.matrix(sample_broad_pct)

## add missing columns as 0 if any cell type is absent
missing_cols <- setdiff(celltype_order, colnames(sample_broad_pct_df))
sample_broad_pct_df[missing_cols] <- 0

sample_broad_pct_df <- sample_broad_pct_df[, celltype_order, drop = FALSE] %>%
  rownames_to_column("Sample_Name") %>%
  mutate(across(all_of(celltype_order), ~ round(.x, 2)))

print(sample_broad_pct_df)

write.csv(
  sample_broad_pct_df,
  file = "4R_WT_Ctrl_broad_celltype_percent_by_sample.csv",
  row.names = FALSE
)




#percentage of each cell type was exported out to be plotted in Prism




xeno_wt <- AAA_final_celltype_RELN 



## ----- Make a 4R WT control object -----
## Keep only 4R_WT_Ctrl samples based on Sample_Name
xeno_4R <- subset(
  xeno_wt,
  subset = grepl("^4R_WT_Ctrl", Sample_Name)
)

## Drop unused factor levels in metadata
xeno_4R$Sample_Name <- droplevels(factor(xeno_4R$Sample_Name))
xeno_4R$Condition <- droplevels(factor(xeno_4R$Condition))
xeno_4R$control_celltype <- droplevels(factor(xeno_4R$control_celltype))

## Also reset Seurat identities if needed
Idents(xeno_4R) <- droplevels(xeno_4R$control_celltype)

## Check
table(xeno_4R$Sample_Name)
table(xeno_4R$Condition)
unique(xeno_4R$control_celltype)

DefaultAssay(xeno_4R) <- "RNA"



#MG markers analysis

library(Seurat)
library(ggplot2)
library(patchwork)

# -----------------------------
# 1. Subset microglia
# -----------------------------
mg <- subset(
  xeno_4R,
  subset = control_celltype == "MG"
)

DefaultAssay(mg) <- "RNA"

# -----------------------------
# 2. Individual genes to show on UMAP
# -----------------------------

core_mg_markers <- c(
  "CX3CR1",
  "CSF1R"
)

in_vivo_maturation_markers <- c(
  "SALL1",
  "P2RY13",
  "GPR34",
  "SLCO2B1"
)

immature_prolif_markers <- c(
  "MKI67",
  "TOP2A"
)

individual_markers <- c(
  core_mg_markers,
  in_vivo_maturation_markers,
  immature_prolif_markers
)

# -----------------------------
# 3. Paper based mature human microglia score Fig 2E
# -----------------------------
# Top 30 enriched genes identified in mature human microglia from the paper PMID: 37500887


paper_mg_maturity_genes <- c(
  "SPP1",
  "CD74",
  "ACTB",
  "C3",
  "FTL",
  "FOS",
  "CSF1R",
  "B2M",
  "C1QC",
  "C1QB",
  "PSAP",
  "A2M",
  "ITM2B",
  "LAPTM5",
  "CTSB",
  "P2RY12",
  "C1QA",
  "SLCO2B1",
  "RGS1",
  "APOE",
  "CCL4L2",
  "RNASET2",
  "NEAT1",
  "CX3CR1",
  "DUSP1",
  "SAT1",
  "ZFP36",
  "CD81",
  "HLA-B",
  "HLA-DRA"
)

# Keep only genes present in your object
paper_mg_maturity_genes_use <- intersect(
  paper_mg_maturity_genes,
  rownames(mg)
)

missing_paper_maturity_genes <- setdiff(
  paper_mg_maturity_genes,
  rownames(mg)
)

print("Genes used for paper based MG maturity score:")
print(paper_mg_maturity_genes_use)

print("Genes missing from object:")
print(missing_paper_maturity_genes)

# Add module score
set.seed(42)

mg <- AddModuleScore(
  object = mg,
  features = list(paper_mg_maturity_genes_use),
  name = "Paper_MG_Maturity_Score"
)

# This creates:
# Paper_MG_Maturity_Score1

# -----------------------------
# 4. Final 3 by 3 panel order
# -----------------------------
features_to_plot <- c(
  "CX3CR1",
  "CSF1R",
  "SALL1",
  "P2RY13",
  "GPR34",
  "SLCO2B1",
  "MKI67",
  "TOP2A",
  "Paper_MG_Maturity_Score1"
)

features_to_plot <- features_to_plot[
  features_to_plot %in% c(rownames(mg), colnames(mg@meta.data))
]

# -----------------------------
# 5. Plot individual marker genes
# -----------------------------
gene_features <- setdiff(features_to_plot, "Paper_MG_Maturity_Score1")

gene_plots <- FeaturePlot(
  object = mg,
  features = gene_features,
  reduction = "umap",
  order = TRUE,
  pt.size = 1.5,
  combine = FALSE
)

gene_plots <- lapply(seq_along(gene_plots), function(i) {
  
  gene_plots[[i]] +
    scale_color_gradient(
      low = "lightgrey",
      high = "red",
      na.value = "grey90"
    ) +
    NoAxes() +
    ggtitle(gene_features[i]) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 20
      ),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "right",
      legend.title = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )
})

# -----------------------------
# 6. Plot paper based maturation score
# -----------------------------
score_plot <- FeaturePlot(
  object = mg,
  features = "Paper_MG_Maturity_Score1",
  reduction = "umap",
  order = TRUE,
  pt.size = 1.5,
  combine = FALSE
)[[1]] +
  scale_color_gradient(
    low = "lightgrey",
    high = "red",
    na.value = "grey90"
  ) +
  NoAxes() +
  ggtitle("MG maturity score") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 20
    ),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# -----------------------------
# 7. Combine into final 3 by 3 panel
# -----------------------------
final_plots <- c(gene_plots, list(score_plot))

p_mg_maturation_3x3 <- wrap_plots(
  final_plots,
  ncol = 3,
  nrow = 3
) +
  plot_annotation(
    title = "Human microglia identity and in vivo maturation",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )
  )

p_mg_maturation_3x3





##AST maturation 


library(Seurat)
library(ggplot2)
library(patchwork)
library(grid)

# -----------------------------
# 1. Subset astrocytes
# -----------------------------
ast <- subset(
  xeno_4R,
  subset = control_celltype == "AST"
)

DefaultAssay(ast) <- "RNA"




# -----------------------------
# 2.Mature ast markers PMID: 38418648 Fig 2J
# -----------------------------


wang_mature_ast_genes <- c(
  "SLC1A2",
  "SPARCL1",
  "SLC1A3",
  "GLUL",
  "CPE",
  "GJA1",
  "ATP1B2",
  "AQP4",
  "GPM6A",
  "ATP1A2",
  "GLUD1"
)





ast_maturation_genes<-
  wang_mature_ast_genes
ast_maturation_genes_use <- intersect(ast_maturation_genes, rownames(ast))
setdiff(ast_maturation_genes, rownames(ast))


genes_present <- intersect(wang_mature_ast_genes, rownames(ast))
genes_missing <- setdiff(wang_mature_ast_genes, rownames(ast))

genes_missing
length(genes_present)



# Fetch normalized expression
expr_mat <- FetchData(
  ast,
  vars = genes_present,
  layer = "data"
)

# Gene-level expression summary across all astrocytes
gene_expr_summary <- tibble(
  gene = colnames(expr_mat),
  mean_expr = colMeans(expr_mat, na.rm = TRUE),
  median_expr = apply(expr_mat, 2, median, na.rm = TRUE),
  pct_expressing = colMeans(expr_mat > 0, na.rm = TRUE) * 100,
  max_expr = apply(expr_mat, 2, max, na.rm = TRUE)
) %>%
  arrange(mean_expr)

gene_expr_summary


gene_expr_summary %>%
  arrange(pct_expressing)

low_detected_genes <- gene_expr_summary %>%
  filter(pct_expressing < 5 | mean_expr < 0.05) %>%
  arrange(pct_expressing)

low_detected_genes

# -----------------------------
# 3. Add astrocyte maturation score
# -----------------------------
set.seed(42)

ast <- AddModuleScore(
  object = ast,
  features = list(ast_maturation_genes_use),
  name = "AST_Maturation_Score"
)

# -----------------------------
# 4. Individual markers for 8-panel row
# -----------------------------
ast_identity_maturation_markers <- c(
  "ALDH1L1",
  "AQP4",
  "GJA1",
  "SLC1A2",
  "SPARCL1"
)

proliferation_markers <- c(
  "MKI67",
  "TOP2A"
)

gene_features <- c(
  ast_identity_maturation_markers,
  proliferation_markers
)

gene_features <- intersect(gene_features, rownames(ast))

# -----------------------------
# 5. Shared style
# -----------------------------
purple_high <- "#7B3F98"

small_legend <- guides(
  color = guide_colorbar(
    barheight = unit(2, "cm"),
    barwidth  = unit(.5, "cm"),
    ticks = FALSE
  )
)

panel_theme <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width = unit(0.25, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, -3, 0, -3),
    plot.margin = margin(1, 1, 1, 1),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# -----------------------------
# 6. Gene plots
# -----------------------------
gene_plots <- FeaturePlot(
  object = ast,
  features = gene_features,
  reduction = "umap",
  order = TRUE,
  pt.size = 1.5,
  combine = FALSE
)

gene_plots <- lapply(seq_along(gene_plots), function(i) {
  gene_plots[[i]] +
    scale_color_gradient(
      low = "lightgrey",
      high = purple_high,
      na.value = "grey90"
    ) +
    small_legend +
    NoAxes() +
    ggtitle(gene_features[i]) +
    panel_theme
})

# -----------------------------
# 7. Maturation score plot
# -----------------------------
score_plot <- FeaturePlot(
  object = ast,
  features = "AST_Maturation_Score1",
  reduction = "umap",
  order = TRUE,
  pt.size = 1.5,
  combine = FALSE
)[[1]] +
  scale_color_gradient(
    low = "lightgrey",
    high = purple_high,
    na.value = "grey90"
  ) +
  small_legend +
  NoAxes() +
  ggtitle("AST maturation score") +
  panel_theme

# -----------------------------
# 8. Combine into one row of 8
# -----------------------------
final_plots <- c(gene_plots, list(score_plot))

p_ast_maturation_8panel <- wrap_plots(
  final_plots,
  nrow = 1
) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(hjust = 0.5,  size = 14)
    )
  ) &
  theme(
    plot.margin = margin(1, 1, 1, 1)
  )

p_ast_maturation_8panel

