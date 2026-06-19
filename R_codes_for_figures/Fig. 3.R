
## ------------------------------------------------------------
## Transfer WT annotations (control_celltype / control_superclass)
## to xeno_5genotype (exact ID copy first, then anchor-based transfer)
## ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
})

## 0) Inputs
xeno_5genotype <- engraft_5Genotypes_noKI_annotated.by.6.Humandatasets_recluster_subclusters_Celltype.updated.assigned
ref <- xeno_wt
qry <- xeno_5genotype

stopifnot(inherits(ref, "Seurat"), inherits(qry, "Seurat"))

cols_to_transfer <- c("control_celltype", "control_superclass")

## 1) (Optional) Recompute PCA, neighbors, clusters for xeno_5genotype
##     Only do this if you truly need re-clustering in this object.
DefaultAssay(qry) <- "integrated"   # adjust if needed

qry <- RunPCA(qry, verbose = FALSE)
qry <- FindNeighbors(qry, reduction = "pca", dims = 1:30, verbose = FALSE)

graph_name <- paste0(DefaultAssay(qry), "_snn")
res_use <- 0.3

qry <- FindClusters(
  object     = qry,
  graph.name = graph_name,
  resolution = res_use,
  verbose    = FALSE
)

ident_col <- paste0(graph_name, "_res.", res_use)
Idents(qry) <- ident_col

# Optional quick visualization
DimPlot(qry, group.by = ident_col, label = TRUE, repel = TRUE)

## 2) Helper: robust cell name getter
get_cells <- function(obj) {
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    return(SeuratObject::Cells(obj))
  }
  rownames(obj@meta.data)
}

## 3) Ensure target metadata columns exist
for (cc in cols_to_transfer) {
  if (!cc %in% colnames(qry@meta.data)) {
    qry@meta.data[[cc]] <- NA_character_
  }
}

## 4) Strategy A: direct copy for overlapping cell IDs
ref_cells <- get_cells(ref)
qry_cells <- get_cells(qry)
overlap_cells <- intersect(ref_cells, qry_cells)

if (length(overlap_cells) > 0) {
  meta_src <- ref@meta.data[overlap_cells, cols_to_transfer, drop = FALSE]
  meta_src <- data.frame(
    lapply(meta_src, function(v) if (is.factor(v)) as.character(v) else v),
    row.names = overlap_cells,
    check.names = FALSE
  )
  qry@meta.data[overlap_cells, cols_to_transfer] <- meta_src
  message(sprintf("Exact ID transfer completed for %d overlapping cells.", length(overlap_cells)))
} else {
  message("No overlapping cell IDs found between reference and query; skipping exact copy.")
}

## 5) Strategy B: anchor-based label transfer for remaining NA cells
need_transfer <- is.na(qry@meta.data[["control_celltype"]]) | is.na(qry@meta.data[["control_superclass"]])

if (any(need_transfer)) {
  norm_method <- if ("SCT" %in% Assays(ref) && "SCT" %in% Assays(qry)) "SCT" else "LogNormalize"
  
  anchors <- FindTransferAnchors(
    reference            = ref,
    query                = qry,
    normalization.method = norm_method,
    dims                 = 1:30,
    verbose              = FALSE
  )
  
  pred_celltype <- TransferData(
    anchorset = anchors,
    refdata   = ref$control_celltype,
    dims      = 1:30,
    verbose   = FALSE
  )
  
  pred_superclass <- TransferData(
    anchorset = anchors,
    refdata   = ref$control_superclass,
    dims      = 1:30,
    verbose   = FALSE
  )
  
  idx_na_ct <- is.na(qry@meta.data[["control_celltype"]])
  idx_na_sc <- is.na(qry@meta.data[["control_superclass"]])
  
  qry@meta.data[idx_na_ct, "control_celltype"]   <- pred_celltype$predicted.id[idx_na_ct]
  qry@meta.data[idx_na_sc, "control_superclass"] <- pred_superclass$predicted.id[idx_na_sc]
  
  # Store scores for QC (useful for thresholding / filtering later)
  qry$control_celltype_score   <- pred_celltype$prediction.score.max
  qry$control_superclass_score <- pred_superclass$prediction.score.max
  
  message(sprintf("Anchor-based transfer filled labels for %d cells.", sum(need_transfer)))
} else {
  message("No missing labels detected; anchor-based transfer not needed.")
}

## 6) (Optional) Harmonize factor levels to match reference
lvl_ct <- if (is.factor(ref$control_celltype)) levels(ref$control_celltype) else unique(as.character(ref$control_celltype))
lvl_sc <- if (is.factor(ref$control_superclass)) levels(ref$control_superclass) else unique(as.character(ref$control_superclass))

qry$control_celltype   <- factor(as.character(qry$control_celltype),   levels = lvl_ct)
qry$control_superclass <- factor(as.character(qry$control_superclass), levels = lvl_sc)

## 7) Write back + sanity checks + plots
xeno_5genotype <- qry

cat("\ncontrol_celltype counts:\n")
print(table(xeno_5genotype$control_celltype, useNA = "ifany"))

cat("\ncontrol_superclass counts:\n")
print(table(xeno_5genotype$control_superclass, useNA = "ifany"))

DimPlot(xeno_5genotype, group.by = "control_celltype", label = TRUE, repel = TRUE) + NoAxes()
DimPlot(xeno_5genotype, group.by = "control_superclass", label = TRUE, repel = TRUE) + NoAxes()

## 8) Save
saveRDS(xeno_5genotype, file = "xeno_5genotype_final.rds")


#4R only seurat

xeno_4genotype_no3r_RELNetc_high_col <- subset(xeno_5genotype, subset = Condition != "3R_WT_Ctrl")



#Fig. 3I
library(Seurat)
library(dplyr)
library(ggplot2)

# Object
obj <- xeno_4genotype_no3r_RELNetc_high_col #seurat with all 4R conditions 

# Check condition and cell type columns
table(obj$Condition)
table(obj$control_celltype)

# 1. Count nuclei per condition
condition_count_df <- obj@meta.data %>%
  count(Condition, name = "N_nuclei") %>%
  arrange(Condition)

print(condition_count_df)

# Save as csv
write.csv(
  condition_count_df,
  file = "nuclei_count_by_condition.csv",
  row.names = FALSE
)


# count nuclei per condition and cell type
condition_celltype_count_df <- obj@meta.data %>%
  count(Condition, control_celltype, name = "N_nuclei") %>%
  arrange(Condition, control_celltype)

print(condition_celltype_count_df)

write.csv(
  condition_celltype_count_df,
  file = "nuclei_count_by_condition_and_celltype.csv",
  row.names = FALSE
)

#entire umap
DimPlot(
  obj,
  reduction = "umap",
  group.by = "control_celltype",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.25
) +
  NoAxes()

#only label neurons
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

obj <- xeno_4genotype_no3r_RELNetc_high_col

# Extract UMAP coordinates
umap_df <- as.data.frame(Embeddings(obj, reduction = "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

# Add metadata
umap_df$control_celltype <- obj$control_celltype
umap_df$is_neuron <- grepl("^Neu_", umap_df$control_celltype)

# Split neuron and non-neuron nuclei
non_neu_df <- umap_df %>%
  filter(!is_neuron)

neu_df <- umap_df %>%
  filter(is_neuron)

# Label positions for neuron subtypes only
neu_label_df <- neu_df %>%
  group_by(control_celltype) %>%
  summarize(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

# Make neuron subtype colors
neu_types <- sort(unique(neu_df$control_celltype))
neu_cols <- hue_pal()(length(neu_types))
names(neu_cols) <- neu_types

# Plot
# Label positions for neuron subtypes only
neu_label_df <- neu_df %>%
  group_by(control_celltype) %>%
  summarize(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label_clean = gsub("^Neu_", "", control_celltype)
  )

# Plot
p_neuron_subtype_umap <- ggplot() +
  geom_point(
    data = non_neu_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = "grey85",
    size = 0.25,
    alpha = 0.5
  ) +
  geom_point(
    data = neu_df,
    aes(x = UMAP_1, y = UMAP_2, color = control_celltype),
    size = 0.25,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = neu_label_df,
    aes(x = UMAP_1, y = UMAP_2, label = label_clean),
    size = 5,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = neu_cols) +
  NoAxes() +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) 

p_neuron_subtype_umap




#Fig. 3J
## ------------------------------------------------------------
## Normalized nuclei count (drop 3R first)
## Outcome per sample:
##   percent_neuron = 100 * (Neuron nuclei) / (Total human nuclei)
## Then summarize by Condition for plotting/statistics.
## ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

## 0) Input object
obj <- xeno_5genotype   
md <- obj@meta.data

## 1) Drop 3R at the beginning (keep only 4R conditions)
keep_levels <- c("4R_WT_Ctrl", "4R_WT_K18", "4R_P301S_Ctrl", "4R_P301S_K18")

stopifnot("Condition" %in% colnames(md))
present <- intersect(keep_levels, unique(md$Condition))
stopifnot(length(present) > 0)

obj4r <- subset(obj, subset = Condition %in% present)
md4r  <- obj4r@meta.data

## Optional: enforce plotting order
md4r$Condition <- factor(md4r$Condition, levels = keep_levels)

## 2) Define what counts as “Neuron nuclei”
## Preferred: use control_superclass if it has EN/IN/Neu_Early
## Fallback: use control_celltype with "Neu_" prefix
has_super <- "control_superclass" %in% colnames(md4r)
has_ct    <- "control_celltype"   %in% colnames(md4r)

stopifnot(has_super || has_ct)

is_neuron <- rep(FALSE, nrow(md4r))

if (has_super && any(md4r$control_superclass %in% c("EN", "IN", "Neu_Early"), na.rm = TRUE)) {
  is_neuron <- md4r$control_superclass %in% c("EN", "IN", "Neu_Early")
} else if (has_ct) {
  is_neuron <- str_detect(as.character(md4r$control_celltype), "^Neu")
} else if (has_super) {
  ## if your object already collapsed EN/IN into "Neu"
  is_neuron <- as.character(md4r$control_superclass) %in% c("Neu", "Neu_Early")
}

## 3) Per-sample normalized neuron nuclei percent
stopifnot("Sample_Name" %in% colnames(md4r))

per_sample <- md4r %>%
  transmute(
    Sample_Name = as.character(Sample_Name),
    Condition   = as.character(Condition),
    is_neuron   = is_neuron
  ) %>%
  group_by(Sample_Name, Condition) %>%
  summarise(
    total_human_nuclei = n(),
    neuron_nuclei      = sum(is_neuron, na.rm = TRUE),
    percent_neuron     = 100 * neuron_nuclei / total_human_nuclei,
    .groups = "drop"
  ) %>%
  mutate(
    Condition = factor(Condition, levels = keep_levels)
  )

print(per_sample)

## Save the exact numbers used for the dots in the bar plot
write.csv(per_sample, "normalized_neuron_percent_per_sample_4Ronly.csv", row.names = FALSE)

#Data were replotted in Prism



#Fig. 3K


neu_obj <- subset(
  x = obj_no3r,
  subset = broad_celltype == "Neuron"
)

#DEGs
## ---------- Neu total ----------
DefaultAssay(neu_obj) <- "RNA"
Idents(neu_obj) <- "Condition"
unique(neu_obj$Condition)


NEU_4R301k_vs_4R301ctrl <- FindMarkers(neu_obj, ident.1 = "4R_P301S_K18", ident.2 = "4R_P301S_Ctrl",
                                       logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
NEU_4R301k_vs_4R301ctrl <- subset(NEU_4R301k_vs_4R301ctrl, p_val_adj < 0.05)
NEU_4R301k_vs_4R301ctrl$gene <- rownames(NEU_4R301k_vs_4R301ctrl)
NEU_4R301k_vs_4R301ctrl <- NEU_4R301k_vs_4R301ctrl[order(abs(NEU_4R301k_vs_4R301ctrl$avg_log2FC), decreasing = TRUE),]
write.csv(NEU_4R301k_vs_4R301ctrl, file.path(out_dir, "NEU_4R301k_vs_4R301ctrl.csv"), row.names = FALSE)




## Fig 3K: pathfindR KEGG term–gene network + export upregulated genes

# =========================================================
# 0) Libraries
# =========================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(pathfindR)
  library(patchwork)  # for wrap_plots()
})

# =========================================================
# 1) Input and thresholds
# =========================================================
Neu_301kvsct   <- NEU_4R301k_vs_4R301ctrl
padj_cutoff    <- 0.05
logfc_abs_cut  <- 0.1

# =========================================================
# 2) Build pathfindR input (Gene.symbol, logFC, adj.P.Val)
# =========================================================
Neu_301k <- Neu_301kvsct %>%
  transmute(
    Gene.symbol = as.character(gene),
    logFC       = as.numeric(avg_log2FC),
    adj.P.Val   = as.numeric(p_val_adj)
  ) %>%
  filter(!is.na(Gene.symbol), !is.na(logFC), !is.na(adj.P.Val)) %>%
  arrange(adj.P.Val) %>%
  distinct(Gene.symbol, .keep_all = TRUE)   # keep best padj per gene if duplicated

RA_input_301k <- Neu_301k %>%
  filter(adj.P.Val < padj_cutoff, abs(logFC) > logfc_abs_cut)

if (nrow(RA_input_301k) == 0) stop("No genes passed the thresholds.")

# =========================================================
# 3) Run pathfindR (KEGG)
# =========================================================
set.seed(301)
RA_output_301k <- run_pathfindR(
  input            = RA_input_301k,
  gene_sets        = "KEGG",         # alternatives: "Reactome", "GO-BP", etc.
  convert2alias    = FALSE,          # gene names already in symbols (per your note)
  pin_name_path    = "Biogrid",
  p_val_threshold  = 0.05
)

write.csv(
  RA_output_301k,
  file = "RA_output_301k_pathfindR_KEGG.csv",
  row.names = FALSE
)

# =========================================================
# 4) Term–gene graphs (all terms, then selected terms)
# =========================================================
if (nrow(RA_output_301k) > 0) {
  term_gene_graph(result_df = RA_output_301k, use_description = TRUE)
}

sel1 <- c("Apoptosis", "Necroptosis", "Cellular senescence", "Ferroptosis")
selected_res_path1 <- RA_output_301k %>% filter(Term_Description %in% sel1)

if (nrow(selected_res_path1) > 0) {
  term_gene_graph(
    result_df        = selected_res_path1,
    use_description  = TRUE,
    node_colors      = c("#E5D7BF", "red", "blue")
  )
}

# =========================================================
# 5) Label control for ggrepel + add node labels with ggraph
# =========================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# Let ggrepel place all labels by default
ggplot2::update_geom_defaults("text_repel",  list(max.overlaps = Inf))
ggplot2::update_geom_defaults("label_repel", list(max.overlaps = Inf))

library(ggraph)   # for geom_node_text

sel1 <- c("Apoptosis", "Necroptosis", "Cellular senescence", "Ferroptosis")
selected_res_path1 <- RA_output_301k %>% dplyr::filter(Term_Description %in% sel1)


if (nrow(selected_res_path1) > 0) {
  p <- term_gene_graph(
    result_df       = selected_res_path1,
    use_description = TRUE,
    node_colors     = c("#E5D7BF", "red", "blue")
  )
  
  p <- p + ggraph::geom_node_text(
    aes(label = name),
    repel = TRUE,
    max.overlaps = Inf,
    size = 5,
    
    # --- in-between spacing (closer than before, but still readable) ---
    point.padding = grid::unit(0.6, "lines"),
    box.padding   = grid::unit(0.45, "lines"),
    force         = 12,
    force_pull    = 0.6,
    
    # --- connector line visible but subtle ---
    min.segment.length = 0,
    segment.alpha = 0.5,
    segment.size  = 0.3
  )
  
  print(p)
}

# =========================================================
# 6) Export UPREGULATED genes for selected pathways (4 terms) for spatial
# =========================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

up_genes_long <- selected_res_path1 %>%
  mutate(Up_regulated = as.character(Up_regulated)) %>%
  select(Term_Description, Up_regulated) %>%
  tidyr::separate_rows(Up_regulated, sep = "\\s*,\\s*") %>%   # split comma-separated genes
  filter(!is.na(Up_regulated), Up_regulated != "") %>%
  rename(Gene.symbol = Up_regulated) %>%
  distinct()

# OPTIONAL: attach logFC + adj.P.Val from your pathfindR input (RA_input_301k)
up_genes_long_stats <- up_genes_long %>%
  left_join(
    RA_input_301k %>% select(Gene.symbol, logFC, adj.P.Val),
    by = "Gene.symbol"
  ) %>%
  arrange(Term_Description, adj.P.Val)

# 1) A unique union list of upregulated genes across the 4 pathways
up_genes_union <- sort(unique(up_genes_long_stats$Gene.symbol))

# Print results
up_genes_union
up_genes_long_stats

# 2) (Optional) genes grouped by pathway (nice compact view)
up_genes_by_pathway <- up_genes_long_stats %>%
  group_by(Term_Description) %>%
  summarise(
    n_up     = n_distinct(Gene.symbol),
    Up_genes = paste(sort(unique(Gene.symbol)), collapse = ", "),
    .groups  = "drop"
  )

up_genes_by_pathway

# 3) Save outputs
write.csv(up_genes_long_stats, "sel1_upregulated_genes_long.csv", row.names = FALSE)
write.csv(up_genes_by_pathway, "sel1_upregulated_genes_by_pathway.csv", row.names = FALSE)
writeLines(up_genes_union, "sel1_upregulated_genes_union.txt")

# =========================================================
# 7) Export UPREGULATED genes for selected pathways (3 terms)
# =========================================================
sel1 <- c("Apoptosis", "Necroptosis", "Ferroptosis")
selected_res_path1 <- RA_output_301k %>% dplyr::filter(Term_Description %in% sel1)

up_genes_long <- selected_res_path1 %>%
  mutate(Up_regulated = as.character(Up_regulated)) %>%
  select(Term_Description, Up_regulated) %>%
  tidyr::separate_rows(Up_regulated, sep = "\\s*,\\s*") %>%   # split comma-separated genes
  filter(!is.na(Up_regulated), Up_regulated != "") %>%
  rename(Gene.symbol = Up_regulated) %>%
  distinct()

# OPTIONAL: attach logFC + adj.P.Val from your pathfindR input (RA_input_301k)
up_genes_long_stats <- up_genes_long %>%
  left_join(
    RA_input_301k %>% select(Gene.symbol, logFC, adj.P.Val),
    by = "Gene.symbol"
  ) %>%
  arrange(Term_Description, adj.P.Val)

# 1) A unique union list of upregulated genes across the 3 pathways
up_genes_union <- sort(unique(up_genes_long_stats$Gene.symbol))

# Print results
up_genes_union
up_genes_long_stats

# 2) (Optional) genes grouped by pathway (nice compact view)
up_genes_by_pathway <- up_genes_long_stats %>%
  group_by(Term_Description) %>%
  summarise(
    n_up     = n_distinct(Gene.symbol),
    Up_genes = paste(sort(unique(Gene.symbol)), collapse = ", "),
    .groups  = "drop"
  )

up_genes_by_pathway

# 3) Save outputs
write.csv(up_genes_long_stats, "sel1_upregulated_genes_long.csv", row.names = FALSE)
write.csv(up_genes_by_pathway, "sel1_upregulated_genes_by_pathway.csv", row.names = FALSE)
writeLines(up_genes_union, "celldeath_upregulated_genes_union.txt")


#genes were further uploaded to Loupe browser for spatial transcriptomics


#Fig 3 L


##Cell count for globel 

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

# ---- 0) Input + subset to the two 4R P301S conditions ----
keep_conditions <- c("4R_P301S_Ctrl", "4R_P301S_K18")

obj <- subset(
  xeno_5genotype,
  subset = Condition %in% keep_conditions
)



# ---- 1) Settings for "high" nuclei ----
q            <- 0.90  # top 10%
genes        <- c("SLC6A5", "RELN", "KCNC2", "SLC17A6")
assay_use    <- DefaultAssay(obj)
by_condition <- FALSE  # TRUE = top 10% within each Condition; FALSE = global top 10%

# ---- 2) Seurat v4/v5 compatible expression getter (prefers normalized data) ----
.get_mat <- function(obj, assay, layer) {
  out <- tryCatch(GetAssayData(obj, assay = assay, layer = layer), error = function(e) NULL)  # v5
  if (!is.null(out)) return(out)
  out <- tryCatch(GetAssayData(obj, assay = assay, slot  = layer), error = function(e) NULL)  # v4
  if (!is.null(out)) return(out)
  tryCatch(LayerData(obj, assay = assay, layer = layer), error = function(e) NULL)            # v5 alt
}

.get_expr_vec <- function(obj, gene, assay) {
  mat_data   <- .get_mat(obj, assay, "data")
  mat_counts <- .get_mat(obj, assay, "counts")
  mat <- if (!is.null(mat_data) && nrow(mat_data) > 0) mat_data else mat_counts
  if (is.null(mat) || !(gene %in% rownames(mat))) return(NULL)
  v <- mat[gene, , drop = TRUE]
  names(v) <- colnames(mat)
  v
}

# ---- 3) Identify high nuclei per gene, then summarize by Sample_Name (n + %) ----
meta <- obj@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  transmute(
    cell_id,
    Sample_Name = as.character(Sample_Name),
    Condition   = as.character(Condition)
  )

flag_list <- list()

for (g in genes) {
  expr <- .get_expr_vec(obj, g, assay_use)
  if (is.null(expr)) {
    message(sprintf("Skipping %s: not found in assay '%s'.", g, assay_use))
    next
  }
  
  df <- tibble(
    cell_id = names(expr),
    expr    = as.numeric(expr)
  ) %>%
    left_join(meta, by = "cell_id")
  
  if (isTRUE(by_condition)) {
    df <- df %>%
      group_by(Condition) %>%
      mutate(thr = quantile(expr[is.finite(expr)], probs = q, na.rm = TRUE)) %>%
      ungroup()
    high <- df$expr > df$thr
  } else {
    thr  <- quantile(df$expr[is.finite(df$expr)], probs = q, na.rm = TRUE)
    high <- df$expr > thr
  }
  
  high[!is.finite(df$expr)] <- FALSE
  
  flag_list[[g]] <- tibble(
    Sample_Name = df$Sample_Name,
    gene        = g,
    high        = high
  )
}



gene_high_by_sample <- bind_rows(flag_list) %>%
  group_by(Sample_Name, gene) %>%
  summarise(
    n_high   = sum(high, na.rm = TRUE),
    n_total  = dplyr::n(),
    pct_high = 100 * n_high / n_total,
    .groups  = "drop"
  ) %>%
  arrange(gene, Sample_Name)

print(gene_high_by_sample)



write.csv(gene_high_by_sample, "gene_high_by_sample.csv")







suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

# -----------------------------
# Inputs
# -----------------------------
keep_conditions <- c("4R_P301S_Ctrl", "4R_P301S_K18")

obj <- subset(
  xeno_5genotype,
  subset = Condition %in% keep_conditions
)


keep_conditions <- c("4R_P301S_Ctrl", "4R_P301S_K18")
genes <- c("SLC6A5", "RELN", "KCNC2", "SLC17A6")
q <- 0.90
by_condition <- FALSE  # FALSE = global top 10% across both conditions 
out_dir <- "DEG_high_markers"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Choose the object you want to (1) label and (2) run DEG on
# For neuron-only DEG, set this to neu_obj
obj_use <- obj

stopifnot(inherits(obj_use, "Seurat"))

# Keep only the two conditions (matches your logic)
obj_use <- subset(obj_use, subset = Condition %in% keep_conditions)

# Ensure Condition exists and is character (helps robust comparisons)
obj_use$Condition <- as.character(obj_use$Condition)

# -----------------------------
# Seurat v4/v5 compatible expression getter (same as yours)
# -----------------------------
.get_mat <- function(obj, assay, layer) {
  out <- tryCatch(GetAssayData(obj, assay = assay, layer = layer), error = function(e) NULL)  # v5
  if (!is.null(out)) return(out)
  out <- tryCatch(GetAssayData(obj, assay = assay, slot  = layer), error = function(e) NULL)  # v4
  if (!is.null(out)) return(out)
  tryCatch(LayerData(obj, assay = assay, layer = layer), error = function(e) NULL)            # v5 alt
}

.get_expr_vec <- function(obj, gene, assay) {
  mat_data   <- .get_mat(obj, assay, "data")
  mat_counts <- .get_mat(obj, assay, "counts")
  mat <- if (!is.null(mat_data) && nrow(mat_data) > 0) mat_data else mat_counts
  if (is.null(mat) || !(gene %in% rownames(mat))) return(NULL)
  v <- mat[gene, , drop = TRUE]
  names(v) <- colnames(mat)
  v
}

# -----------------------------
# 1) Add <GENE>_high metadata columns 
# -----------------------------
assay_use <- DefaultAssay(obj_use)

# store thresholds for record keeping
obj_use@misc$high_thresholds <- list(q = q, by_condition = by_condition, assay = assay_use, thresholds = list())

for (g in genes) {
  expr <- .get_expr_vec(obj_use, g, assay_use)
  if (is.null(expr)) {
    message(sprintf("Skipping %s: not found in assay '%s'.", g, assay_use))
    next
  }
  
  expr <- as.numeric(expr)
  names(expr) <- colnames(obj_use)
  
  # your logic: non-finite treated as FALSE
  finite_expr <- expr[is.finite(expr)]
  
  if (length(finite_expr) == 0) {
    message(sprintf("Skipping %s: no finite expression values.", g))
    next
  }
  
  if (isTRUE(by_condition)) {
    # threshold within each condition (optional mode)
    thr_vec <- rep(NA_real_, length(expr))
    names(thr_vec) <- names(expr)
    
    for (cond in keep_conditions) {
      cells_c <- colnames(obj_use)[obj_use$Condition == cond]
      v_c <- expr[cells_c]
      v_c <- v_c[is.finite(v_c)]
      if (length(v_c) == 0) next
      thr_c <- quantile(v_c, probs = q, na.rm = TRUE)
      thr_vec[cells_c] <- thr_c
    }
    
    high <- expr > thr_vec
    high[!is.finite(expr)] <- FALSE
    
    # store one threshold per condition for this gene
    obj_use@misc$high_thresholds$thresholds[[g]] <- tapply(
      expr,
      obj_use$Condition,
      function(v) quantile(v[is.finite(v)], probs = q, na.rm = TRUE)
    )
    
  } else {
    # your default: global threshold across both conditions combined
    thr <- quantile(finite_expr, probs = q, na.rm = TRUE)
    high <- expr > thr
    high[!is.finite(expr)] <- FALSE
    
    obj_use@misc$high_thresholds$thresholds[[g]] <- thr
  }
  
  flag_col <- paste0(g, "_high")
  obj_use[[flag_col]] <- factor(ifelse(high, "Yes", "No"), levels = c("No", "Yes"))
}

# quick sanity check
for (g in genes) {
  flag_col <- paste0(g, "_high")
  if (flag_col %in% colnames(obj_use@meta.data)) {
    cat("\n", flag_col, "\n")
    print(table(obj_use@meta.data[[flag_col]], useNA = "ifany"))
  }
}


for (g in genes) {
  flag_col <- paste0(g, "_high")
  cat("\n", flag_col, "\n")
  tab <- table(obj_use$Sample_Name, obj_use@meta.data[[flag_col]])
  print(tab)
  cat("\nPercent Yes within each Sample_Name:\n")
  print(round(100 * prop.table(tab, margin = 1)[, "Yes"], 2))
}

# -----------------------------
# 2) DEG within each gene-high subset: K18 vs Ctrl
# -----------------------------
DefaultAssay(obj_use) <- "RNA"
Idents(obj_use) <- "Condition"

run_deg_high <- function(obj, gene, ident.1 = "4R_P301S_K18", ident.2 = "4R_P301S_Ctrl",
                         logfc.threshold = 0.1, min.pct = 0.1, test.use = "MAST",
                         out_dir = ".", p_adj_cutoff = 0.05) {
  
  flag_col <- paste0(gene, "_high")
  if (!flag_col %in% colnames(obj@meta.data)) {
    message(sprintf("Skipping %s DEG: missing %s.", gene, flag_col))
    return(NULL)
  }
  
  cells_high <- rownames(obj@meta.data)[obj@meta.data[[flag_col]] == "Yes"]
  if (length(cells_high) < 20) {
    message(sprintf("Skipping %s DEG: too few high cells (n = %d).", gene, length(cells_high)))
    return(NULL)
  }
  
  obj_high <- subset(obj, cells = cells_high)
  Idents(obj_high) <- "Condition"
  
  # ensure both groups exist in this subset
  n1 <- sum(Idents(obj_high) == ident.1)
  n2 <- sum(Idents(obj_high) == ident.2)
  if (n1 < 10 || n2 < 10) {
    message(sprintf("Skipping %s DEG: not enough cells in both groups (K18 n=%d, Ctrl n=%d).", gene, n1, n2))
    return(NULL)
  }
  
  res <- FindMarkers(
    obj_high,
    ident.1 = ident.1,
    ident.2 = ident.2,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    min.pct = min.pct,
    only.pos = FALSE
  )
  
  if (nrow(res) == 0) {
    message(sprintf("No markers returned for %s (high subset).", gene))
    return(NULL)
  }
  
  # standardize FC column name across Seurat versions
  fc_col <- if ("avg_log2FC" %in% colnames(res)) "avg_log2FC" else if ("avg_logFC" %in% colnames(res)) "avg_logFC" else NULL
  if (is.null(fc_col)) fc_col <- colnames(res)[1]  # fallback, rarely needed
  
  res <- res %>%
    tibble::rownames_to_column("gene") %>%
    filter(p_val_adj < p_adj_cutoff) %>%
    arrange(desc(abs(.data[[fc_col]])))
  
  out_csv <- file.path(out_dir, sprintf("NEU_%s_high__4R301k_vs_4R301ctrl.csv", gene))
  write.csv(res, out_csv, row.names = FALSE)
  message("Wrote: ", normalizePath(out_csv))
  
  return(res)
}

deg_list <- lapply(genes, function(g) run_deg_high(obj_use, g, out_dir = out_dir))
names(deg_list) <- genes




# =============================================================================
# Global-high nuclei DEG (MAST) and overlap analysis
# =============================================================================

genes_high <- c("RELN","KCNC2","SLC17A6","SLC6A5")
out_dir    <- "DEG_Neu_high_by_gene"
logfc_thr  <- 0.1
min_pct    <- 0.1
cond_1     <- "4R_P301S_K18"
cond_2     <- "4R_P301S_Ctrl"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 1) Subset neurons: control_celltype starts with "Neu_"
neu_obj <- subset(
  x = obj,
  subset = grepl("^Neu_", control_celltype)
)
DefaultAssay(neu_obj)<-"RNA"

# Helper: run DE for a given <GENE>_high == "Yes"
run_deg_for_high <- function(seu, gene_name) {
  flag <- paste0(gene_name, "_high")
  if (!flag %in% colnames(seu@meta.data)) {
    warning(sprintf("Missing '%s' in metadata. Skipping %s.", flag, gene_name))
    return(invisible(NULL))
  }
  
  # Cells with GENE_high == "Yes" 
  cells_yes <- rownames(seu@meta.data)[seu@meta.data[[flag]] == "Yes"]
  if (length(cells_yes) < 20) {
    warning(sprintf("%s_high has only %d cells after neuron subset. Skipping.",
                    gene_name, length(cells_yes)))
    return(invisible(NULL))
  }
  
  seu_high <- subset(seu, cells = cells_yes)
  
  # Set identities to Condition and keep only the two groups
  Idents(seu_high) <- "Condition"
  keep_ids <- intersect(levels(Idents(seu_high)), c(cond_1, cond_2))
  if (!all(c(cond_1, cond_2) %in% keep_ids)) {
    warning(sprintf("%s_high: one/both conditions missing after filtering. Skipping.", gene_name))
    return(invisible(NULL))
  }
  seu_high <- subset(seu_high, idents = c(cond_1, cond_2))
  
  # DEG with MAST
  res <- FindMarkers(
    object          = seu_high,
    ident.1         = cond_1,
    ident.2         = cond_2,
    logfc.threshold = logfc_thr,
    test.use        = "MAST",
    min.pct         = min_pct,
    only.pos        = FALSE
  )
  
  if (nrow(res) == 0) {
    warning(sprintf("%s_high: FindMarkers returned 0 rows.", gene_name))
    return(invisible(NULL))
  }
  
  res <- res %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::arrange(dplyr::desc(abs(avg_log2FC)))
  
  out_file <- file.path(out_dir, sprintf("DEG_Neu_%shigh_4R301k_vs_4R301ctrl.csv", gene_name))
  write.csv(res, out_file, row.names = FALSE)
  message("Wrote: ", normalizePath(out_file))
  invisible(res)
}

deg_results <- lapply(genes_high, function(g) run_deg_for_high(neu_obj, g))
names(deg_results) <- genes_high




DEG_Neu_SLC17A6high_4R301k_vs_4R301ctrl<-as.data.frame(NEU_SLC17A6_high_4R301k_vs_4R301ctrl)
DEG_Neu_SLC6A5high_4R301k_vs_4R301ctrl<-as.data.frame(NEU_SLC6A5_high_4R301k_vs_4R301ctrl)
DEG_Neu_RELNhigh_4R301k_vs_4R301ctrl<-as.data.frame(NEU_RELN_high_4R301k_vs_4R301ctrl)
DEG_Neu_KCNC2high_4R301k_vs_4R301ctrl<-as.data.frame(NEU_KCNC2_high_4R301k_vs_4R301ctrl)



## =============================================================================
## DEG overlap (Venn) + common UP/DOWN export + GO BP enrichment (clusterProfiler)
## Inputs expected in memory:
##   DEG_Neu_SLC17A6high_4R301k_vs_4R301ctrl
##   DEG_Neu_SLC6A5high_4R301k_vs_4R301ctrl
##   DEG_Neu_RELNhigh_4R301k_vs_4R301ctrl
##   DEG_Neu_KCNC2high_4R301k_vs_4R301ctrl
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggVennDiagram)
})

## ----------------------------- 0) Output dir ---------------------------------
out_dir <- "C:/Users/wenhq/Documents/test"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ----------------------------- 1) Helpers ------------------------------------
# Ensure a "gene" column exists
.ensure_gene_col <- function(df) {
  stopifnot(is.data.frame(df))
  if (!"gene" %in% colnames(df)) df <- tibble::rownames_to_column(df, "gene")
  df %>% filter(!is.na(gene), gene != "")
}

# Extract UP and DOWN gene vectors based on avg_log2FC sign
get_up_down <- function(df) {
  df <- .ensure_gene_col(df) %>% filter(!is.na(avg_log2FC))
  list(
    up   = df %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique(),
    down = df %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()
  )
}

# Save a simple gene list
.write_gene_list <- function(genes, filename) {
  write.csv(
    data.frame(gene = unique(as.character(genes))),
    file.path(out_dir, filename),
    row.names = FALSE
  )
}

## ----------------------------- 2) Inputs -------------------------------------
deg_list <- list(
  SLC17A6high = DEG_Neu_SLC17A6high_4R301k_vs_4R301ctrl,
  SLC6A5high  = DEG_Neu_SLC6A5high_4R301k_vs_4R301ctrl,
  RELNhigh    = DEG_Neu_RELNhigh_4R301k_vs_4R301ctrl,
  KCNC2high   = DEG_Neu_KCNC2high_4R301k_vs_4R301ctrl
)

## ------------------ 3) Build UP/DOWN sets for Venn ---------------------------
sets <- lapply(deg_list, get_up_down)
gene_sets_up   <- lapply(sets, `[[`, "up")
gene_sets_down <- lapply(sets, `[[`, "down")

## ----------------------------- 4) Venn plots ---------------------------------
venn_up <- ggVennDiagram(gene_sets_up, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "indianred3") +
  labs(title = "Overlap of UP-regulated DEGs across neuronal groups") +
  theme(plot.title = element_text(size = 14, face = "bold"))
print(venn_up)

venn_down <- ggVennDiagram(gene_sets_down, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Overlap of DOWN-regulated DEGs across neuronal groups") +
  theme(plot.title = element_text(size = 14, face = "bold"))
print(venn_down)

#venn plots were exported out and combined using Affinity software




###DEG anlaysis 
suppressPackageStartupMessages(library(Seurat))
DefaultAssay(xeno_5genotype) <- "RNA"

# Subset by control_superclass
xeno_MG  <- subset(xeno_5genotype, control_superclass == "MG")
xeno_EN  <- subset(xeno_5genotype, control_superclass == "EN")
xeno_IN  <- subset(xeno_5genotype, control_superclass == "IN")
xeno_AST <- subset(xeno_5genotype, control_superclass == "AST")

out_dir <- "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type"

## ---------- EN ----------
Idents(xeno_EN) <- "Condition"
unique(xeno_EN$Condition)

EN_4r_vs_3r <- FindMarkers(xeno_EN, ident.1 = "4R_WT_Ctrl", ident.2 = "3R_WT_Ctrl",
                           logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
EN_4r_vs_3r <- subset(EN_4r_vs_3r, p_val_adj < 0.05)
EN_4r_vs_3r$gene <- rownames(EN_4r_vs_3r)
EN_4r_vs_3r <- EN_4r_vs_3r[order(abs(EN_4r_vs_3r$avg_log2FC), decreasing = TRUE),]
write.csv(EN_4r_vs_3r, file.path(out_dir, "EN_4r_vs_3r.csv"), row.names = FALSE)

EN_4R301_vs_4RWT_ctrl <- FindMarkers(xeno_EN, ident.1 = "4R_P301S_Ctrl", ident.2 = "4R_WT_Ctrl",
                                     logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
EN_4R301_vs_4RWT_ctrl <- subset(EN_4R301_vs_4RWT_ctrl, p_val_adj < 0.05)
EN_4R301_vs_4RWT_ctrl$gene <- rownames(EN_4R301_vs_4RWT_ctrl)
EN_4R301_vs_4RWT_ctrl <- EN_4R301_vs_4RWT_ctrl[order(abs(EN_4R301_vs_4RWT_ctrl$avg_log2FC), decreasing = TRUE),]
write.csv(EN_4R301_vs_4RWT_ctrl, file.path(out_dir, "EN_4R301_vs_4RWT_ctrl.csv"), row.names = FALSE)

EN_4Rwtk_vs_4RWTctrl <- FindMarkers(xeno_EN, ident.1 = "4R_WT_K18", ident.2 = "4R_WT_Ctrl",
                                    logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
EN_4Rwtk_vs_4RWTctrl <- subset(EN_4Rwtk_vs_4RWTctrl, p_val_adj < 0.05)
EN_4Rwtk_vs_4RWTctrl$gene <- rownames(EN_4Rwtk_vs_4RWTctrl)
EN_4Rwtk_vs_4RWTctrl <- EN_4Rwtk_vs_4RWTctrl[order(abs(EN_4Rwtk_vs_4RWTctrl$avg_log2FC), decreasing = TRUE),]
write.csv(EN_4Rwtk_vs_4RWTctrl, file.path(out_dir, "EN_4Rwtk_vs_4RWTctrl.csv"), row.names = FALSE)

EN_4R301k_vs_4R301ctrl <- FindMarkers(xeno_EN, ident.1 = "4R_P301S_K18", ident.2 = "4R_P301S_Ctrl",
                                      logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
EN_4R301k_vs_4R301ctrl <- subset(EN_4R301k_vs_4R301ctrl, p_val_adj < 0.05)
EN_4R301k_vs_4R301ctrl$gene <- rownames(EN_4R301k_vs_4R301ctrl)
EN_4R301k_vs_4R301ctrl <- EN_4R301k_vs_4R301ctrl[order(abs(EN_4R301k_vs_4R301ctrl$avg_log2FC), decreasing = TRUE),]
write.csv(EN_4R301k_vs_4R301ctrl, file.path(out_dir, "EN_4R301k_vs_4R301ctrl.csv"), row.names = FALSE)


## ---------- IN ----------
Idents(xeno_IN) <- "Condition"
unique(xeno_IN$Condition)

IN_4r_vs_3r <- FindMarkers(xeno_IN, ident.1 = "4R_WT_Ctrl", ident.2 = "3R_WT_Ctrl",
                           logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
IN_4r_vs_3r <- subset(IN_4r_vs_3r, p_val_adj < 0.05)
IN_4r_vs_3r$gene <- rownames(IN_4r_vs_3r)
IN_4r_vs_3r <- IN_4r_vs_3r[order(abs(IN_4r_vs_3r$avg_log2FC), decreasing = TRUE),]
write.csv(IN_4r_vs_3r, file.path(out_dir, "IN_4r_vs_3r.csv"), row.names = FALSE)

IN_4R301_vs_4RWT_ctrl <- FindMarkers(xeno_IN, ident.1 = "4R_P301S_Ctrl", ident.2 = "4R_WT_Ctrl",
                                     logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
IN_4R301_vs_4RWT_ctrl <- subset(IN_4R301_vs_4RWT_ctrl, p_val_adj < 0.05)
IN_4R301_vs_4RWT_ctrl$gene <- rownames(IN_4R301_vs_4RWT_ctrl)
IN_4R301_vs_4RWT_ctrl <- IN_4R301_vs_4RWT_ctrl[order(abs(IN_4R301_vs_4RWT_ctrl$avg_log2FC), decreasing = TRUE),]
write.csv(IN_4R301_vs_4RWT_ctrl, file.path(out_dir, "IN_4R301_vs_4RWT_ctrl.csv"), row.names = FALSE)

IN_4Rwtk_vs_4RWTctrl <- FindMarkers(xeno_IN, ident.1 = "4R_WT_K18", ident.2 = "4R_WT_Ctrl",
                                    logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
IN_4Rwtk_vs_4RWTctrl <- subset(IN_4Rwtk_vs_4RWTctrl, p_val_adj < 0.05)
IN_4Rwtk_vs_4RWTctrl$gene <- rownames(IN_4Rwtk_vs_4RWTctrl)
IN_4Rwtk_vs_4RWTctrl <- IN_4Rwtk_vs_4RWTctrl[order(abs(IN_4Rwtk_vs_4RWTctrl$avg_log2FC), decreasing = TRUE),]
write.csv(IN_4Rwtk_vs_4RWTctrl, file.path(out_dir, "IN_4Rwtk_vs_4RWTctrl.csv"), row.names = FALSE)

IN_4R301k_vs_4R301ctrl <- FindMarkers(xeno_IN, ident.1 = "4R_P301S_K18", ident.2 = "4R_P301S_Ctrl",
                                      logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
IN_4R301k_vs_4R301ctrl <- subset(IN_4R301k_vs_4R301ctrl, p_val_adj < 0.05)
IN_4R301k_vs_4R301ctrl$gene <- rownames(IN_4R301k_vs_4R301ctrl)
IN_4R301k_vs_4R301ctrl <- IN_4R301k_vs_4R301ctrl[order(abs(IN_4R301k_vs_4R301ctrl$avg_log2FC), decreasing = TRUE),]
write.csv(IN_4R301k_vs_4R301ctrl, file.path(out_dir, "IN_4R301k_vs_4R301ctrl.csv"), row.names = FALSE)

## ---------- MG ----------
Idents(xeno_MG) <- "Condition"
unique(xeno_MG$Condition)

MG_4r_vs_3r <- FindMarkers(xeno_MG, ident.1 = "4R_WT_Ctrl", ident.2 = "3R_WT_Ctrl",
                           logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
MG_4r_vs_3r <- subset(MG_4r_vs_3r, p_val_adj < 0.05)
MG_4r_vs_3r$gene <- rownames(MG_4r_vs_3r)
MG_4r_vs_3r <- MG_4r_vs_3r[order(abs(MG_4r_vs_3r$avg_log2FC), decreasing = TRUE),]
write.csv(MG_4r_vs_3r, file.path(out_dir, "MG_4r_vs_3r.csv"), row.names = FALSE)

MG_4R301_vs_4RWT_ctrl <- FindMarkers(xeno_MG, ident.1 = "4R_P301S_Ctrl", ident.2 = "4R_WT_Ctrl",
                                     logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
MG_4R301_vs_4RWT_ctrl <- subset(MG_4R301_vs_4RWT_ctrl, p_val_adj < 0.05)
MG_4R301_vs_4RWT_ctrl$gene <- rownames(MG_4R301_vs_4RWT_ctrl)
MG_4R301_vs_4RWT_ctrl <- MG_4R301_vs_4RWT_ctrl[order(abs(MG_4R301_vs_4RWT_ctrl$avg_log2FC), decreasing = TRUE),]
write.csv(MG_4R301_vs_4RWT_ctrl, file.path(out_dir, "MG_4R301_vs_4RWT_ctrl.csv"), row.names = FALSE)

MG_4Rwtk_vs_4RWTctrl <- FindMarkers(xeno_MG, ident.1 = "4R_WT_K18", ident.2 = "4R_WT_Ctrl",
                                    logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
MG_4Rwtk_vs_4RWTctrl <- subset(MG_4Rwtk_vs_4RWTctrl, p_val_adj < 0.05)
MG_4Rwtk_vs_4RWTctrl$gene <- rownames(MG_4Rwtk_vs_4RWTctrl)
MG_4Rwtk_vs_4RWTctrl <- MG_4Rwtk_vs_4RWTctrl[order(abs(MG_4Rwtk_vs_4RWTctrl$avg_log2FC), decreasing = TRUE),]
write.csv(MG_4Rwtk_vs_4RWTctrl, file.path(out_dir, "MG_4Rwtk_vs_4RWTctrl.csv"), row.names = FALSE)

MG_4R301k_vs_4R301ctrl <- FindMarkers(xeno_MG, ident.1 = "4R_P301S_K18", ident.2 = "4R_P301S_Ctrl",
                                      logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
MG_4R301k_vs_4R301ctrl <- subset(MG_4R301k_vs_4R301ctrl, p_val_adj < 0.05)
MG_4R301k_vs_4R301ctrl$gene <- rownames(MG_4R301k_vs_4R301ctrl)
MG_4R301k_vs_4R301ctrl <- MG_4R301k_vs_4R301ctrl[order(abs(MG_4R301k_vs_4R301ctrl$avg_log2FC), decreasing = TRUE),]
write.csv(MG_4R301k_vs_4R301ctrl, file.path(out_dir, "MG_4R301k_vs_4R301ctrl.csv"), row.names = FALSE)

## ---------- AST ----------
Idents(xeno_AST) <- "Condition"
unique(xeno_AST$Condition)

AST_4r_vs_3r <- FindMarkers(xeno_AST, ident.1 = "4R_WT_Ctrl", ident.2 = "3R_WT_Ctrl",
                            logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
AST_4r_vs_3r <- subset(AST_4r_vs_3r, p_val_adj < 0.05)
AST_4r_vs_3r$gene <- rownames(AST_4r_vs_3r)
AST_4r_vs_3r <- AST_4r_vs_3r[order(abs(AST_4r_vs_3r$avg_log2FC), decreasing = TRUE),]
write.csv(AST_4r_vs_3r, file.path(out_dir, "AST_4r_vs_3r.csv"), row.names = FALSE)

AST_4R301_vs_4RWT_ctrl <- FindMarkers(xeno_AST, ident.1 = "4R_P301S_Ctrl", ident.2 = "4R_WT_Ctrl",
                                      logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
AST_4R301_vs_4RWT_ctrl <- subset(AST_4R301_vs_4RWT_ctrl, p_val_adj < 0.05)
AST_4R301_vs_4RWT_ctrl$gene <- rownames(AST_4R301_vs_4RWT_ctrl)
AST_4R301_vs_4RWT_ctrl <- AST_4R301_vs_4RWT_ctrl[order(abs(AST_4R301_vs_4RWT_ctrl$avg_log2FC), decreasing = TRUE),]
write.csv(AST_4R301_vs_4RWT_ctrl, file.path(out_dir, "AST_4R301_vs_4RWT_ctrl.csv"), row.names = FALSE)

AST_4Rwtk_vs_4RWTctrl <- FindMarkers(xeno_AST, ident.1 = "4R_WT_K18", ident.2 = "4R_WT_Ctrl",
                                     logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
AST_4Rwtk_vs_4RWTctrl <- subset(AST_4Rwtk_vs_4RWTctrl, p_val_adj < 0.05)
AST_4Rwtk_vs_4RWTctrl$gene <- rownames(AST_4Rwtk_vs_4RWTctrl)
AST_4Rwtk_vs_4RWTctrl <- AST_4Rwtk_vs_4RWTctrl[order(abs(AST_4Rwtk_vs_4RWTctrl$avg_log2FC), decreasing = TRUE),]
write.csv(AST_4Rwtk_vs_4RWTctrl, file.path(out_dir, "AST_4Rwtk_vs_4RWTctrl.csv"), row.names = FALSE)

AST_4R301k_vs_4R301ctrl <- FindMarkers(xeno_AST, ident.1 = "4R_P301S_K18", ident.2 = "4R_P301S_Ctrl",
                                       logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
AST_4R301k_vs_4R301ctrl <- subset(AST_4R301k_vs_4R301ctrl, p_val_adj < 0.05)
AST_4R301k_vs_4R301ctrl$gene <- rownames(AST_4R301k_vs_4R301ctrl)
AST_4R301k_vs_4R301ctrl <- AST_4R301k_vs_4R301ctrl[order(abs(AST_4R301k_vs_4R301ctrl$avg_log2FC), decreasing = TRUE),]
write.csv(AST_4R301k_vs_4R301ctrl, file.path(out_dir, "AST_4R301k_vs_4R301ctrl.csv"), row.names = FALSE)



