## =============================================================================
## PGRNKI: Transfer / restore control_celltype + control_superclass labels
## Reorganized for publication (same logic, same steps, no new analysis)
## Objects used:
##   - xeno_5genotype
##   - engraft_human_integrated_Annotation_4Genotypes_withKI  -> xeno_wKI
##   - xeno_p301s_labeled_new (reference labels)              -> xeno_p301s_labeled
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

## -----------------------------------------------------------------------------
## 0) Inputs / object setup
## -----------------------------------------------------------------------------
# source objects
xeno_wKI <- engraft_human_integrated_Annotation_4Genotypes_withKI
# xeno_5genotype exists already

# keep only the two P301S conditions as the reference for label transfer
xeno_p301s <- subset(
  xeno_5genotype,
  subset = Condition %in% c("4R_P301S_K18", "4R_P301S_Ctrl")
)

# save reference subset (as you did)
saveRDS(xeno_p301s, "xeno_p301s_labeled_new.RDS")

# reference labeled object (as you defined)
xeno_p301s_labeled <- xeno_p301s_labeled_new

# optional quick checks (as in your script)
unique(xeno_5genotype$Condition)
unique(xeno_wKI$Condition)

## -----------------------------------------------------------------------------
## 1) (Optional) Recluster + UMAP on wKI for visualization
## -----------------------------------------------------------------------------
xeno_wKI <- FindClusters(
  xeno_wKI,
  resolution   = 0.3,
  cluster.name = "res0.3",
  graph.name   = "integrated_snn"
)

xeno_wKI <- RunUMAP(
  xeno_wKI,
  reduction      = "pca",
  dims           = 1:30,
  reduction.name = "umap_res0.3"
)

DimPlot(xeno_wKI, group.by = "res0.3", reduction = "umap_res0.3", label = TRUE, repel = TRUE)

## -----------------------------------------------------------------------------
## 2) Visualize existing labels in reference
## -----------------------------------------------------------------------------
DimPlot(xeno_p301s_labeled, group.by = "control_celltype")

## -----------------------------------------------------------------------------
## 3) Label transfer ONLY for cells that are NA AND are KI_/ISOKI_ conditions
##    (Seurat TransferData workflow)
## -----------------------------------------------------------------------------
needs_fill <- which(
  is.na(xeno_wKI$control_celltype) &
    grepl("^(KI_|ISOKI_)", xeno_wKI$Condition)
)
length(needs_fill)

if (length(needs_fill) > 0) {
  
  query_sub <- subset(xeno_wKI, cells = colnames(xeno_wKI)[needs_fill])
  
  # RPCA if PCA exists in both, else CCA
  use_rpca <- ("pca" %in% names(xeno_p301s_labeled@reductions)) &&
    ("pca" %in% names(query_sub@reductions))
  
  dims_use <- 1:30
  
  # reference labels as character (avoid factor leakage)
  ref_ct <- setNames(as.character(xeno_p301s_labeled$control_celltype),
                     colnames(xeno_p301s_labeled))
  ref_sc <- setNames(as.character(xeno_p301s_labeled$control_superclass),
                     colnames(xeno_p301s_labeled))
  
  anchors <- FindTransferAnchors(
    reference = xeno_p301s_labeled,
    query     = query_sub,
    dims      = dims_use,
    reduction = if (use_rpca) "rpca" else "cca"
  )
  
  pred_ct <- TransferData(anchorset = anchors, refdata = ref_ct, dims = dims_use)
  pred_sc <- TransferData(anchorset = anchors, refdata = ref_sc, dims = dims_use)
  
  # ensure target columns exist and are character
  if (!"control_celltype"   %in% colnames(xeno_wKI@meta.data)) xeno_wKI$control_celltype   <- NA_character_
  if (!"control_superclass" %in% colnames(xeno_wKI@meta.data)) xeno_wKI$control_superclass <- NA_character_
  xeno_wKI$control_celltype   <- as.character(xeno_wKI$control_celltype)
  xeno_wKI$control_superclass <- as.character(xeno_wKI$control_superclass)
  
  if (!"ct_score" %in% colnames(xeno_wKI@meta.data)) xeno_wKI$ct_score <- NA_real_
  if (!"sc_score" %in% colnames(xeno_wKI@meta.data)) xeno_wKI$sc_score <- NA_real_
  
  # write predictions back ONLY to subset cells
  cells_sub <- colnames(query_sub)
  xeno_wKI@meta.data[cells_sub, "control_celltype"]   <- as.character(pred_ct[cells_sub, "predicted.id"])
  xeno_wKI@meta.data[cells_sub, "control_superclass"] <- as.character(pred_sc[cells_sub, "predicted.id"])
  
  # store max scores
  if ("prediction.score.max" %in% colnames(pred_ct)) {
    xeno_wKI@meta.data[cells_sub, "ct_score"] <- pred_ct[cells_sub, "prediction.score.max"]
  } else {
    ct_cols <- grep("^prediction\\.score", colnames(pred_ct), value = TRUE)
    xeno_wKI@meta.data[cells_sub, "ct_score"] <- apply(pred_ct[cells_sub, ct_cols, drop = FALSE], 1, max)
  }
  
  if ("prediction.score.max" %in% colnames(pred_sc)) {
    xeno_wKI@meta.data[cells_sub, "sc_score"] <- pred_sc[cells_sub, "prediction.score.max"]
  } else {
    sc_cols <- grep("^prediction\\.score", colnames(pred_sc), value = TRUE)
    xeno_wKI@meta.data[cells_sub, "sc_score"] <- apply(pred_sc[cells_sub, sc_cols, drop = FALSE], 1, max)
  }
  
  # optional: recast to factors using reference levels
  xeno_wKI$control_celltype <- factor(
    xeno_wKI$control_celltype,
    levels = levels(xeno_p301s_labeled$control_celltype)
  )
  xeno_wKI$control_superclass <- factor(
    xeno_wKI$control_superclass,
    levels = levels(xeno_p301s_labeled$control_superclass)
  )
}

## -----------------------------------------------------------------------------
## 4) Sanity checks + visualization
## -----------------------------------------------------------------------------
unique(xeno_wKI$control_celltype)
table(xeno_wKI$Condition, is.na(xeno_wKI$control_celltype))


DimPlot(xeno_wKI, group.by = "control_celltype", reduction = "umap_res0.3", label = TRUE, repel = TRUE)


## -----------------------------------------------------------------------------
## 5) Restore P301S labels by exact mapping (Sample_Name + CoreID join)
##    (same mapping logic you used)
## -----------------------------------------------------------------------------
# diagnose NA distribution in P301S
table(xeno_wKI$Condition, is.na(xeno_wKI$control_celltype))

# helper: take the part BEFORE the first underscore (10x core id)
core_id <- function(x) sub("_.*$", "", x)

ref_keys <- xeno_p301s_labeled@meta.data %>%
  mutate(
    .cell_ref = rownames(.),
    CoreID    = core_id(.cell_ref),
    control_celltype   = as.character(control_celltype),
    control_superclass = as.character(control_superclass)
  ) %>%
  dplyr::select(Sample_Name, CoreID, .cell_ref, control_celltype, control_superclass)

wki_keys <- xeno_wKI@meta.data %>%
  mutate(
    .cell_wki = rownames(.),
    CoreID    = core_id(.cell_wki)
  ) %>%
  dplyr::select(Sample_Name, Condition, .cell_wki, CoreID)

m_p301s <- inner_join(
  wki_keys %>% filter(Condition %in% c("4R_P301S_Ctrl", "4R_P301S_K18")),
  ref_keys,
  by = c("Sample_Name", "CoreID")
)

# write back as character (avoid factor code issues)
if (!"control_celltype"   %in% colnames(xeno_wKI@meta.data)) xeno_wKI$control_celltype   <- NA_character_
if (!"control_superclass" %in% colnames(xeno_wKI@meta.data)) xeno_wKI$control_superclass <- NA_character_
xeno_wKI$control_celltype   <- as.character(xeno_wKI$control_celltype)
xeno_wKI$control_superclass <- as.character(xeno_wKI$control_superclass)

xeno_wKI@meta.data[m_p301s$.cell_wki, "control_celltype"]   <- m_p301s$control_celltype
xeno_wKI@meta.data[m_p301s$.cell_wki, "control_superclass"] <- m_p301s$control_superclass

# optional: factor levels using UNION so valid labels never become NA
ref_lvls_ct <- levels(xeno_p301s_labeled$control_celltype)
ref_lvls_sc <- levels(xeno_p301s_labeled$control_superclass)

xeno_wKI$control_celltype <- factor(
  xeno_wKI$control_celltype,
  levels = union(ref_lvls_ct, sort(unique(xeno_wKI$control_celltype)))
)
xeno_wKI$control_superclass <- factor(
  xeno_wKI$control_superclass,
  levels = union(ref_lvls_sc, sort(unique(xeno_wKI$control_superclass)))
)

# sanity checks
table(xeno_wKI$Condition, is.na(xeno_wKI$control_celltype))
unique(xeno_wKI$control_celltype)

## -----------------------------------------------------------------------------
## 6) Report overlap counts (your diagnostic block)
## -----------------------------------------------------------------------------
ref_keys2 <- xeno_p301s_labeled@meta.data %>%
  mutate(.cell_ref = rownames(.),
         CoreID    = core_id(.cell_ref)) %>%
  dplyr::select(Sample_Name, Condition, .cell_ref, CoreID,
                control_celltype, control_superclass)

wki_keys2 <- xeno_wKI@meta.data %>%
  mutate(.cell_wki = rownames(.),
         CoreID    = core_id(.cell_wki)) %>%
  filter(Condition %in% c("4R_P301S_K18", "4R_P301S_Ctrl")) %>%
  dplyr::select(Sample_Name, Condition, .cell_wki, CoreID)

matches <- inner_join(
  wki_keys2,
  ref_keys2 %>% dplyr::select(Sample_Name, CoreID, .cell_ref),
  by = c("Sample_Name", "CoreID")
)

cat("n_ref =", nrow(ref_keys2), "\n")
cat("n_wki (2 conds) =", nrow(wki_keys2), "\n")
cat("n_matches (by Sample_Name + CoreID) =", nrow(matches), "\n")

if (nrow(matches) > 0) {
  ref_lab <- ref_keys2 %>%
    dplyr::select(.cell_ref, control_celltype, control_superclass)
  
  match_lab <- matches %>%
    left_join(ref_lab, by = c(".cell_ref")) %>%
    dplyr::select(.cell_wki, control_celltype, control_superclass)
  
  if (!"control_celltype" %in% colnames(xeno_wKI@meta.data)) xeno_wKI$control_celltype <- NA
  if (!"control_superclass" %in% colnames(xeno_wKI@meta.data)) xeno_wKI$control_superclass <- NA
  
  xeno_wKI@meta.data[match_lab$.cell_wki, c("control_celltype", "control_superclass")] <-
    match_lab[, c("control_celltype", "control_superclass")]
}

print(table(xeno_wKI$Condition, is.na(xeno_wKI$control_celltype)))

## -----------------------------------------------------------------------------
## 7) Save final labeled object
## -----------------------------------------------------------------------------
saveRDS(xeno_wKI, "xeno_wKI_final_label.RDS")

#Fig 5F


# ============================================================
# Use the same UMAP reduction as your original KI plot
# ============================================================

reduction_use <- "umap_res0.3"

# Check available reductions if needed
Reductions(obj)

# ============================================================
# 3. Entire UMAP, all cell types labeled
# ============================================================

p_KI_entire_umap <- DimPlot(
  obj,
  reduction = reduction_use,
  group.by = "control_celltype",
  label = TRUE,
  repel = TRUE,
  pt.size = .5
) +
  NoAxes() +
  ggtitle("KI graft cell type annotation") +
  theme(
    legend.title = element_blank()
  )

p_KI_entire_umap

# ggsave(
#   filename = file.path(outdir, "KI_entire_umap_celltype_labeled.pdf"),
#   plot = p_KI_entire_umap,
#   width = 6,
#   height = 5
# )

# ============================================================
# 4. Neuron subtype UMAP
# Non-neuronal nuclei are gray
# Only neuron subtypes are colored and labeled
# Neu_ prefix is removed from labels
# ============================================================

# Extract UMAP coordinates using the correct KI UMAP reduction
umap_df <- as.data.frame(Embeddings(obj, reduction = reduction_use))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

# Add metadata
umap_df$cell <- rownames(umap_df)
umap_df$control_celltype <- obj$control_celltype
umap_df$Condition <- obj$Condition
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
  ) %>%
  mutate(
    label_clean = gsub("^Neu_", "", control_celltype)
  )

# Make neuron subtype colors
neu_types <- sort(unique(neu_df$control_celltype))
neu_cols <- hue_pal()(length(neu_types))
names(neu_cols) <- neu_types

# Plot neuron subtype UMAP
p_KI_neuron_subtype_umap <- ggplot() +
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
  coord_fixed() +
  ggtitle("KI neuron subtype representation") +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_KI_neuron_subtype_umap

ggsave(
  filename = file.path(outdir, "KI_neuron_subtype_umap_gray_non_neurons.pdf"),
  plot = p_KI_neuron_subtype_umap,
  width = 6,
  height = 5
)




##Fig 5G: Neu nuclei percentage (Neu = EN + IN + Neu_Early)
DefaultAssay(xeno_wKI)<-"RNA"
unique(xeno_wKI$control_superclass)


## Neu nuclei percentage (Neu = EN + IN + Neu_Early)
## Uses your Seurat object: xeno_wKI

suppressPackageStartupMessages({
  library(Seurat)
  library(dittoSeq)
  library(dplyr)
  library(ggplot2)
})



# -----------------------------
# 0) Add combined superclass (EN + IN + Neu_Early -> Neu)
# -----------------------------
obj <- xeno_wKI

obj$combined_superclass <- as.character(obj$control_superclass)
obj$combined_superclass[obj$combined_superclass %in% c("EN", "IN", "Neu_Early")] <- "Neu"
obj$combined_superclass <- factor(obj$combined_superclass)

# -----------------------------
# 1) Per-sample percentages (Sample_Name on x-axis)
# -----------------------------
df_raw <- dittoBarPlot(
  object   = obj,
  var      = "combined_superclass",
  group.by = "Sample_Name",
  scale    = "percent",
  data.out = TRUE
)$data
# expected columns include: label, grouping, count, label.count.total.per.facet, percent

# -----------------------------
# 2) Ensure percent is numeric (compute if missing)
# -----------------------------
if (!"percent" %in% names(df_raw)) {
  df_raw$percent <- with(df_raw, 100 * count / label.count.total.per.facet)
}
df_raw$percent <- suppressWarnings(as.numeric(df_raw$percent))

# -----------------------------
# 3) Keep Neu only + attach Condition
# -----------------------------
condition_map <- obj@meta.data %>%
  distinct(Sample_Name, Condition)

neu_df <- df_raw %>%
  filter(label == "Neu") %>%
  transmute(
    Sample  = grouping,
    Class   = label,
    Percent = percent
  ) %>%
  left_join(condition_map, by = c("Sample" = "Sample_Name")) %>%
  filter(!is.na(Condition))

# (Optional) enforce a specific condition order
cond_levels <- c("4R_P301S_Ctrl", "4R_P301S_K18",  "ISOKI_4R_P301S_K18", "KI_4R_P301S_K18")
neu_df$Condition <- factor(neu_df$Condition, levels = intersect(cond_levels, unique(neu_df$Condition)))

# -----------------------------
# 4) Condition means (bars)
# -----------------------------
cond_means <- neu_df %>%
  group_by(Condition) %>%
  summarize(Percent = mean(Percent, na.rm = TRUE), .groups = "drop")

cond_means$Condition <- factor(cond_means$Condition, levels = levels(neu_df$Condition))

# -----------------------------
# 5) Plot (bars = mean, dots = samples)
# -----------------------------
p_neu_bar <- ggplot() +
  geom_col(
    data  = cond_means,
    aes(x = Condition, y = Percent),
    width = 0.7,
    alpha = 0.7
  ) +
  geom_point(
    data = neu_df,
    aes(x = Condition, y = Percent),
    position = position_jitter(width = 0.12, height = 0),
    size = 2
  ) +
  labs(
    x = "Condition",
    y = "Neu (EN + IN + Neu_Early) nuclei percent",
    title = "Neu nuclei percentage by condition"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank())

p_neu_bar

# -----------------------------
# 6) Export CSVs and plotted in Prism
# -----------------------------
write.csv(
  neu_df[order(neu_df$Condition, neu_df$Sample), ],
  file = "Neu_sample_level_points.csv",
  row.names = FALSE
)


###DEGs
unique(xeno_wKI$Condition)

# separate different cell types
DefaultAssay(xeno_wKI) <- "RNA"

PGRN_KI_MG  <- subset(xeno_wKI, control_superclass == "MG",  invert = FALSE)
PGRN_KI_EN  <- subset(xeno_wKI, control_superclass == "EN",  invert = FALSE)
PGRN_KI_IN  <- subset(xeno_wKI, control_superclass == "IN",  invert = FALSE)
PGRN_KI_AST <- subset(xeno_wKI, control_superclass == "AST", invert = FALSE)

# EN_DEG list
Idents(PGRN_KI_EN) <- "Condition"
unique(PGRN_KI_EN$Condition)

EN_301kvsct <- FindMarkers(PGRN_KI_EN, ident.1 = "4R_P301S_K18", ident.2 = "4R_P301S_Ctrl",
                           logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
EN_301kvsct <- subset(EN_301kvsct, EN_301kvsct$p_val_adj < 0.05)
EN_301kvsct$gene <- row.names(EN_301kvsct)
EN_301kvsct <- EN_301kvsct[order(abs(EN_301kvsct$avg_log2FC), decreasing = TRUE),]
write.csv(EN_301kvsct, "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type/PGRN/EN_p301k_vs_ct.csv")

EN_PGRNKIisov301k <- FindMarkers(PGRN_KI_EN, ident.1 = "ISOKI_4R_P301S_K18", ident.2 = "4R_P301S_K18",
                                 logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
EN_PGRNKIisov301k <- subset(EN_PGRNKIisov301k, EN_PGRNKIisov301k$p_val_adj < 0.05)
EN_PGRNKIisov301k$gene <- row.names(EN_PGRNKIisov301k)
EN_PGRNKIisov301k <- EN_PGRNKIisov301k[order(abs(EN_PGRNKIisov301k$avg_log2FC), decreasing = TRUE),]
write.csv(EN_PGRNKIisov301k, "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type/PGRN/EN_PGRNKIiso_vs_301k.csv")

EN_PGRNKI2v301k <- FindMarkers(PGRN_KI_EN, ident.1 = "KI_4R_P301S_K18", ident.2 = "4R_P301S_K18",
                               logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
EN_PGRNKI2v301k <- subset(EN_PGRNKI2v301k, EN_PGRNKI2v301k$p_val_adj < 0.05)
EN_PGRNKI2v301k$gene <- row.names(EN_PGRNKI2v301k)
EN_PGRNKI2v301k <- EN_PGRNKI2v301k[order(abs(EN_PGRNKI2v301k$avg_log2FC), decreasing = TRUE),]
write.csv(EN_PGRNKI2v301k, "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type/PGRN/EN_PGRNKI#2_vs_301k.csv")

# IN_DEG list
Idents(PGRN_KI_IN) <- "Condition"
unique(PGRN_KI_IN$Condition)

IN_301kvsct <- FindMarkers(PGRN_KI_IN, ident.1 = "4R_P301S_K18", ident.2 = "4R_P301S_Ctrl",
                           logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
IN_301kvsct <- subset(IN_301kvsct, IN_301kvsct$p_val_adj < 0.05)
IN_301kvsct$gene <- row.names(IN_301kvsct)
IN_301kvsct <- IN_301kvsct[order(abs(IN_301kvsct$avg_log2FC), decreasing = TRUE),]
write.csv(IN_301kvsct, "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type/PGRN/IN_p301k_vs_ct.csv")

IN_PGRNKIisov301k <- FindMarkers(PGRN_KI_IN, ident.1 = "ISOKI_4R_P301S_K18", ident.2 = "4R_P301S_K18",
                                 logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
IN_PGRNKIisov301k <- subset(IN_PGRNKIisov301k, IN_PGRNKIisov301k$p_val_adj < 0.05)
IN_PGRNKIisov301k$gene <- row.names(IN_PGRNKIisov301k)
IN_PGRNKIisov301k <- IN_PGRNKIisov301k[order(abs(IN_PGRNKIisov301k$avg_log2FC), decreasing = TRUE),]
write.csv(IN_PGRNKIisov301k, "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type/PGRN/IN_PGRNKIiso_vs_301k.csv")

IN_PGRNKI2v301k <- FindMarkers(PGRN_KI_IN, ident.1 = "KI_4R_P301S_K18", ident.2 = "4R_P301S_K18",
                               logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
IN_PGRNKI2v301k <- subset(IN_PGRNKI2v301k, IN_PGRNKI2v301k$p_val_adj < 0.05)
IN_PGRNKI2v301k$gene <- row.names(IN_PGRNKI2v301k)
IN_PGRNKI2v301k <- IN_PGRNKI2v301k[order(abs(IN_PGRNKI2v301k$avg_log2FC), decreasing = TRUE),]
write.csv(IN_PGRNKI2v301k, "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type/PGRN/IN_PGRNKI#2_vs_301k.csv")

# MG_DEG list
Idents(PGRN_KI_MG) <- "Condition"
unique(PGRN_KI_MG$Condition)

MG_301kvsct <- FindMarkers(PGRN_KI_MG, ident.1 = "4R_P301S_K18", ident.2 = "4R_P301S_Ctrl",
                           logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
MG_301kvsct <- subset(MG_301kvsct, MG_301kvsct$p_val_adj < 0.05)
MG_301kvsct$gene <- row.names(MG_301kvsct)
MG_301kvsct <- MG_301kvsct[order(abs(MG_301kvsct$avg_log2FC), decreasing = TRUE),]
write.csv(MG_301kvsct, "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type/PGRN/MG_p301k_vs_ct.csv")

MG_PGRNKIisov301k <- FindMarkers(PGRN_KI_MG, ident.1 = "ISOKI_4R_P301S_K18", ident.2 = "4R_P301S_K18",
                                 logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
MG_PGRNKIisov301k <- subset(MG_PGRNKIisov301k, MG_PGRNKIisov301k$p_val_adj < 0.05)
MG_PGRNKIisov301k$gene <- row.names(MG_PGRNKIisov301k)
MG_PGRNKIisov301k <- MG_PGRNKIisov301k[order(abs(MG_PGRNKIisov301k$avg_log2FC), decreasing = TRUE),]
write.csv(MG_PGRNKIisov301k, "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type/PGRN/MG_PGRNKIiso_vs_301k.csv")

MG_PGRNKI2v301k <- FindMarkers(PGRN_KI_MG, ident.1 = "KI_4R_P301S_K18", ident.2 = "4R_P301S_K18",
                               logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
MG_PGRNKI2v301k <- subset(MG_PGRNKI2v301k, MG_PGRNKI2v301k$p_val_adj < 0.05)
MG_PGRNKI2v301k$gene <- row.names(MG_PGRNKI2v301k)
MG_PGRNKI2v301k <- MG_PGRNKI2v301k[order(abs(MG_PGRNKI2v301k$avg_log2FC), decreasing = TRUE),]
write.csv(MG_PGRNKI2v301k, "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type/PGRN/MG_PGRNKI#2_vs_301k.csv")

# AST_DEG list
Idents(PGRN_KI_AST) <- "Condition"
unique(PGRN_KI_AST$Condition)

AST_301kvsct <- FindMarkers(PGRN_KI_AST, ident.1 = "4R_P301S_K18", ident.2 = "4R_P301S_Ctrl",
                            logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
AST_301kvsct <- subset(AST_301kvsct, AST_301kvsct$p_val_adj < 0.05)
AST_301kvsct$gene <- row.names(AST_301kvsct)
AST_301kvsct <- AST_301kvsct[order(abs(AST_301kvsct$avg_log2FC), decreasing = TRUE),]
write.csv(AST_301kvsct, "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type/PGRN/AST_p301k_vs_ct.csv")

AST_PGRNKIisov301k <- FindMarkers(PGRN_KI_AST, ident.1 = "ISOKI_4R_P301S_K18", ident.2 = "4R_P301S_K18",
                                  logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
AST_PGRNKIisov301k <- subset(AST_PGRNKIisov301k, AST_PGRNKIisov301k$p_val_adj < 0.05)
AST_PGRNKIisov301k$gene <- row.names(AST_PGRNKIisov301k)
AST_PGRNKIisov301k <- AST_PGRNKIisov301k[order(abs(AST_PGRNKIisov301k$avg_log2FC), decreasing = TRUE),]
write.csv(AST_PGRNKIisov301k, "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type/PGRN/AST_PGRNKIiso_vs_301k.csv")

AST_PGRNKI2v301k <- FindMarkers(PGRN_KI_AST, ident.1 = "KI_4R_P301S_K18", ident.2 = "4R_P301S_K18",
                                logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE)
AST_PGRNKI2v301k <- subset(AST_PGRNKI2v301k, AST_PGRNKI2v301k$p_val_adj < 0.05)
AST_PGRNKI2v301k$gene <- row.names(AST_PGRNKI2v301k)
AST_PGRNKI2v301k <- AST_PGRNKI2v301k[order(abs(AST_PGRNKI2v301k$avg_log2FC), decreasing = TRUE),]
write.csv(AST_PGRNKI2v301k, "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type/PGRN/AST_PGRNKI#2_vs_301k.csv")

# ---- NEW: neuron-only subset (control_celltype starts with "Neu_") ----
PGRN_KI_NEU <- subset(xeno_wKI, subset = grepl("^Neu_", control_celltype))

Idents(PGRN_KI_NEU) <- "Condition"
unique(PGRN_KI_NEU$Condition)

# ---- NEW: output folder for neuron-only DEGs ----
out_dir_neu <- "G:/Other computers/My Laptop/sequencing analysis/MS/DEGs_final cell type/PGRN_neurons_only"
if (!dir.exists(out_dir_neu)) dir.create(out_dir_neu, recursive = TRUE, showWarnings = FALSE)

# 4R_P301S_K18 vs 4R_P301S_Ctrl
NEU_301kvsct <- FindMarkers(
  PGRN_KI_NEU, ident.1 = "4R_P301S_K18", ident.2 = "4R_P301S_Ctrl",
  logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE
)
NEU_301kvsct <- subset(NEU_301kvsct, NEU_301kvsct$p_val_adj < 0.05)
NEU_301kvsct$gene <- row.names(NEU_301kvsct)
NEU_301kvsct <- NEU_301kvsct[order(abs(NEU_301kvsct$avg_log2FC), decreasing = TRUE),]
write.csv(NEU_301kvsct, file.path(out_dir_neu, "NEU_p301k_vs_ct.csv"))

# ISOKI_4R_P301S_K18 vs 4R_P301S_K18
NEU_PGRNKIisov301k <- FindMarkers(
  PGRN_KI_NEU, ident.1 = "ISOKI_4R_P301S_K18", ident.2 = "4R_P301S_K18",
  logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE
)
NEU_PGRNKIisov301k <- subset(NEU_PGRNKIisov301k, NEU_PGRNKIisov301k$p_val_adj < 0.05)
NEU_PGRNKIisov301k$gene <- row.names(NEU_PGRNKIisov301k)
NEU_PGRNKIisov301k <- NEU_PGRNKIisov301k[order(abs(NEU_PGRNKIisov301k$avg_log2FC), decreasing = TRUE),]
write.csv(NEU_PGRNKIisov301k, file.path(out_dir_neu, "NEU_PGRNKIiso_vs_301k.csv"))

# KI_4R_P301S_K18 vs 4R_P301S_K18
NEU_PGRNKI2v301k <- FindMarkers(
  PGRN_KI_NEU, ident.1 = "KI_4R_P301S_K18", ident.2 = "4R_P301S_K18",
  logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = FALSE
)
NEU_PGRNKI2v301k <- subset(NEU_PGRNKI2v301k, NEU_PGRNKI2v301k$p_val_adj < 0.05)
NEU_PGRNKI2v301k$gene <- row.names(NEU_PGRNKI2v301k)
NEU_PGRNKI2v301k <- NEU_PGRNKI2v301k[order(abs(NEU_PGRNKI2v301k$avg_log2FC), decreasing = TRUE),]
write.csv(NEU_PGRNKI2v301k, file.path(out_dir_neu, "NEU_PGRNKI#2_vs_301k.csv"))





#Fig 5 H
##DEG analysis
NEU_PGRNKI_2_vs_301k <- as.data.frame(NEU_PGRNKI_2_vs_301k)
NEU_p301k_vs_ct <- as.data.frame(NEU_p301k_vs_ct)
NEU_PGRNKIiso_vs_301k <- as.data.frame(NEU_PGRNKIiso_vs_301k)

colnames(NEU_PGRNKI_2_vs_301k)

suppressPackageStartupMessages({
  library(dplyr)
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    install.packages("VennDiagram")
  }
  library(VennDiagram)
  library(grid)
})

## 1) Define DEG lists (FDR < 0.05) for each comparison
deg1 <- NEU_PGRNKI_2_vs_301k  %>% filter(p_val_adj < 0.05)
deg2 <- NEU_p301k_vs_ct       %>% filter(p_val_adj < 0.05)
deg3 <- NEU_PGRNKIiso_vs_301k %>% filter(p_val_adj < 0.05)

gene_list <- list(
  PGRNKI2_vs_301k    = unique(deg1$gene),
  p301k_vs_ct        = unique(deg2$gene),
  PGRNKIiso_vs_301k  = unique(deg3$gene)
)

## 2) Venn diagram of overlapping DEGs (all, regardless of direction)
grid.newpage()
venn.plot <- venn.diagram(
  x = gene_list,
  filename = NULL,
  fill = c("#FF6A6A", "#CAE1FF", "#C0FF3E"),
  alpha = 0.6,
  category.names = names(gene_list),
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = 0
)
grid.draw(venn.plot)

## 3) Up and down overlaps between NEU_PGRNKI_2_vs_301k and NEU_PGRNKIiso_vs_301k
# significant only
deg1_sig <- deg1
deg3_sig <- deg3

# up and down sets
up1   <- deg1_sig %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
down1 <- deg1_sig %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()

up3   <- deg3_sig %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
down3 <- deg3_sig %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()

# same direction overlaps
up_same   <- intersect(up1,   up3)
down_same <- intersect(down1, down3)

# all genes that move in the same direction in both datasets
same_dir_genes <- union(up_same, down_same)

## 4) Subset these genes from NEU_PGRNKIiso_vs_301k
NEU_PGRNKIiso_vs_301k_sameDir <- NEU_PGRNKIiso_vs_301k %>%
  filter(gene %in% same_dir_genes)



#Fig 3I-J
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})


# ============================================================
#  Helpers (defined if missing)  
# ============================================================
if (!exists("normalize_symbols")) {
  normalize_symbols <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    unique(x[!is.na(x) & nzchar(x)])
  }
}

if (!exists("make_updown")) {
  make_updown <- function(df, fdr = 0.05) {
    df <- df %>%
      mutate(gene = normalize_symbols(gene)) %>%
      filter(!is.na(gene), nzchar(gene)) %>%
      distinct(gene, .keep_all = TRUE)
    
    deg <- df %>% filter(p_val_adj < fdr)
    
    list(
      up   = normalize_symbols(deg$gene[deg$avg_log2FC > 0]),
      down = normalize_symbols(deg$gene[deg$avg_log2FC < 0])
    )
  }
}

if (!exists("do_go_bp")) {
  do_go_bp <- function(genes, p_cut = 0.05, q_cut = 0.2) {
    genes <- normalize_symbols(genes)
    if (length(genes) < 5) return(NULL)
    
    eg <- enrichGO(
      gene          = genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cut,
      qvalueCutoff  = q_cut,
      readable      = TRUE
    )
    if (is.null(eg) || nrow(as.data.frame(eg)) == 0) return(NULL)
    
    clusterProfiler::simplify(eg, cutoff = 0.6, by = "p.adjust", select_fun = min)
  }
}

if (!exists("go_to_df_all")) {
  go_to_df_all <- function(eg) {
    if (is.null(eg)) return(NULL)
    df <- as.data.frame(eg)
    if (nrow(df) == 0) return(NULL)
    df$score <- -log10(df$p.adjust)
    df$label <- sprintf(
      "%s / %s",
      df$Count,
      as.numeric(sub("/.*", "", df$BgRatio))
    )
    df
  }
}

if (!exists("get_rescued_overlap_rank_by_K18")) {
  get_rescued_overlap_rank_by_K18 <- function(df_k18, df_other, top_n = 15) {
    if (is.null(df_k18) || is.null(df_other) || nrow(df_k18) == 0 || nrow(df_other) == 0) return(NULL)
    
    ov <- dplyr::inner_join(
      df_k18   %>% dplyr::select(ID, Description, p.adjust, score, label),
      df_other %>% dplyr::select(ID, Description, p.adjust, score, label),
      by = "ID",
      suffix = c("_k18", "_other")
    )
    if (nrow(ov) == 0) return(NULL)
    
    ov %>%
      dplyr::mutate(Term = stringr::str_trunc(Description_k18, 55)) %>%
      dplyr::arrange(p.adjust_k18, dplyr::desc(score_k18)) %>%
      dplyr::slice_head(n = top_n)
  }
}

# ============================================================
# Plot settings + mirror plot function 
# ============================================================
scheme1 <- list(k18 = "lightpink1",     other = "darkseagreen2") # UpK18/DownRes
scheme2 <- list(k18 = "lightsteelblue1", other = "coral")        # DownK18/UpRes

plot_mirror_ranked_by_K18_signed <- function(ov, title,
                                             k18_sign = +1, other_sign = -1,
                                             k18_color = scheme1$k18, other_color = scheme1$other) {
  if (is.null(ov) || nrow(ov) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste(title, "(no rescued pathways)")))
  }
  ov$Term <- factor(ov$Term, levels = rev(ov$Term))
  
  long <- dplyr::bind_rows(
    ov %>% dplyr::transmute(
      Term,
      Source = "p301k vs ctrl",
      signed_score = k18_sign * score_k18,
      label = label_k18
    ),
    ov %>% dplyr::transmute(
      Term,
      Source = "PGRN OE vs 301k",
      signed_score = other_sign * score_other,
      label = label_other
    )
  )
  
  ggplot(long, aes(x = Term, y = signed_score)) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_col(aes(fill = Source), width = 0.75) +
    geom_text(
      aes(label = label),
      color = "black",
      size  = 5,
      hjust = ifelse(long$signed_score > 0, 1.02, -0.02)
    ) +
    coord_flip() +
    scale_fill_manual(values = c(
      "p301k vs ctrl"   = k18_color,
      "PGRN OE vs 301k" = other_color
    )) +
    labs(x = NULL, y = expression(-log[10]("FDR")), title = title) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      plot.title = element_text(size = 13, face = "bold"),
      legend.position = "none"
    )
}

# ============================================================
# Build DEG sets (exactly your logic)
# ============================================================
neu_k18_sets <- make_updown(NEU_p301k_vs_ct,               fdr = 0.05)
neu_iso_sets <- make_updown(NEU_PGRNKIiso_vs_301k_sameDir, fdr = 0.05)

# ============================================================
# GO BP enrichment (rank side = K18)
# ============================================================
NEU_k18_up_bp   <- do_go_bp(neu_k18_sets$up)
NEU_k18_down_bp <- do_go_bp(neu_k18_sets$down)

NEU_iso_up_bp   <- do_go_bp(neu_iso_sets$up)
NEU_iso_down_bp <- do_go_bp(neu_iso_sets$down)

# Convert to dfs (with scores/labels)
K18_UP_all <- go_to_df_all(NEU_k18_up_bp,   direction = "K18_UP")
K18_DN_all <- go_to_df_all(NEU_k18_down_bp, direction = "K18_DN")
ISO_UP_all <- go_to_df_all(NEU_iso_up_bp,   direction = "ISO_UP")
ISO_DN_all <- go_to_df_all(NEU_iso_down_bp, direction = "ISO_DN")

# ============================================================
#  Overlaps ranked by K18 (ISO sameDir)
#    Case A: K18 Up & Rescue Down   -> scheme1
#    Case B: K18 Down & Rescue Up   -> scheme2
# ============================================================
ov_NEU_UPK18_DNres_ISO <- get_rescued_overlap_rank_by_K18(K18_UP_all, ISO_DN_all, top_n = 10)
ov_NEU_DNK18_UPres_ISO <- get_rescued_overlap_rank_by_K18(K18_DN_all, ISO_UP_all, top_n = 10)

# ============================================================
#  Mirror plots
# ============================================================
p_NEU_ISO_rescue1 <- plot_mirror_ranked_by_K18_signed(
  ov_NEU_UPK18_DNres_ISO,
  "NEU mirror pathways (ISO sameDir): Up in p301k vs ctrl AND Down in PGRN OE vs 301k",
  k18_sign    = +1,
  other_sign  = -1,
  k18_color   = scheme1$k18,
  other_color = scheme1$other
)

p_NEU_ISO_rescue2 <- plot_mirror_ranked_by_K18_signed(
  ov_NEU_DNK18_UPres_ISO,
  "NEU mirror pathways (ISO sameDir): Down in p301k vs ctrl AND Up in PGRN OE vs 301k",
  k18_sign    = -1,
  other_sign  = +1,
  k18_color   = scheme2$k18,
  other_color = scheme2$other
)

# Show
p_NEU_ISO_rescue1
p_NEU_ISO_rescue2





## =============================================================================
## Panel K correlation: highlight putative noncoding RNAs
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)   # stat_cor()
})


df2 <- merge(
  NEU_p301k_vs_ct[, c("gene", "avg_log2FC")],
  NEU_PGRNKIiso_vs_301k_sameDir[, c("gene", "avg_log2FC")],
  by = "gene"
)

df2 <- df2 %>%
  mutate(
    label_to_display = ifelse(
      abs(avg_log2FC.x) > 0.25 & abs(avg_log2FC.y) > 0.25,
      gene,
      NA
    )
  )

## -----------------------------------------------------------------------------
## 1) Noncoding pattern helpers (exact same regex logic you used)
## -----------------------------------------------------------------------------
patterns_any <- function(g) {
  grepl("^ENSG\\d+", g) |
    grepl("-AS\\d*$",  g) |
    grepl("^LINC\\d+", g) |
    grepl("^MIR\\d+",  g) |
    grepl("^SNORD\\d+|^SNORA\\d+|^SCARNA\\d+", g) |
    grepl("^RNU\\d+",  g)
}

## -----------------------------------------------------------------------------
## 2) Count noncoding-like names in ALL genes (df2)
## -----------------------------------------------------------------------------
gene_tbl_all <- df2 %>%
  transmute(gene = as.character(gene)) %>%
  distinct()

N_total_all <- nrow(gene_tbl_all)

counts_all <- gene_tbl_all %>%
  summarise(
    ENSG_ID = sum(grepl("^ENSG\\d+", gene), na.rm = TRUE),
    AS      = sum(grepl("-AS\\d*$",  gene), na.rm = TRUE),
    LINC    = sum(grepl("^LINC\\d+", gene), na.rm = TRUE),
    MIR     = sum(grepl("^MIR\\d+",  gene), na.rm = TRUE),
    SNOR    = sum(grepl("^SNORD\\d+|^SNORA\\d+|^SCARNA\\d+", gene), na.rm = TRUE),
    RNU     = sum(grepl("^RNU\\d+",  gene), na.rm = TRUE)
  ) %>%
  pivot_longer(cols = everything(), names_to = "type", values_to = "n") %>%
  mutate(
    total_unique_genes = N_total_all,
    pct = 100 * n / N_total_all
  ) %>%
  arrange(desc(n))

counts_all

nc_pct_total <- gene_tbl_all %>%
  summarise(
    n_any_nc = sum(patterns_any(gene), na.rm = TRUE),
    total_unique_genes = n(),
    pct_any_nc = 100 * n_any_nc / total_unique_genes
  )

nc_pct_total

## -----------------------------------------------------------------------------
## 3) Recalculate noncoding-like counts for "rescued" (opposite direction) genes only
##    (same exact selection: product < 0)
## -----------------------------------------------------------------------------
df_opposite <- df2 %>%
  mutate(
    gene = as.character(gene),
    avg_log2FC.x = as.numeric(avg_log2FC.x),
    avg_log2FC.y = as.numeric(avg_log2FC.y)
  ) %>%
  filter(!is.na(gene), !is.na(avg_log2FC.x), !is.na(avg_log2FC.y)) %>%
  filter((avg_log2FC.x * avg_log2FC.y) < 0) %>%
  distinct(gene, .keep_all = TRUE)

gene_tbl_opposite <- df_opposite %>% distinct(gene)
N_total_opposite  <- nrow(gene_tbl_opposite)

counts_opposite <- gene_tbl_opposite %>%
  summarise(
    ENSG_ID = sum(grepl("^ENSG\\d+", gene)),
    AS      = sum(grepl("-AS\\d*$",  gene)),
    LINC    = sum(grepl("^LINC\\d+", gene)),
    MIR     = sum(grepl("^MIR\\d+",  gene)),
    SNOR    = sum(grepl("^SNORD\\d+|^SNORA\\d+|^SCARNA\\d+", gene)),
    RNU     = sum(grepl("^RNU\\d+",  gene))
  ) %>%
  pivot_longer(cols = everything(), names_to = "type", values_to = "n") %>%
  mutate(
    total_unique_genes = N_total_opposite,
    pct = round(100 * .data$n / .data$total_unique_genes, 2)
  ) %>%
  arrange(desc(n))

counts_opposite

n_any_opposite <- sum(patterns_any(gene_tbl_opposite$gene), na.rm = TRUE)

counts_opposite_tot <- counts_opposite %>%
  bind_rows(
    tibble(
      type = "TOTAL_any_nc",
      n = n_any_opposite,
      total_unique_genes = N_total_opposite,
      pct = round(100 * n_any_opposite / N_total_opposite, 2)
    )
  ) %>%
  arrange(desc(type == "TOTAL_any_nc"), desc(n))

counts_opposite_tot

## -----------------------------------------------------------------------------
## 4) Panel K plot: highlight putative noncoding RNAs (same plotting logic)
## -----------------------------------------------------------------------------
df2_plot <- df2 %>%
  mutate(
    gene = as.character(gene),
    nc_flag = patterns_any(gene)
  )

df_bg <- df2_plot %>% filter(!nc_flag)  # blue only
df_nc <- df2_plot %>% filter(nc_flag)   # orange only

p_m_nc <- ggplot() +
  geom_point(
    data = df_bg,
    aes(x = avg_log2FC.x, y = avg_log2FC.y),
    color = adjustcolor("#1E90FF", alpha.f = 0.5),
    size  = 5
  ) +
  geom_point(
    data = df_nc,
    aes(x = avg_log2FC.x, y = avg_log2FC.y),
    color = "coral1",
    size  = 5,
    alpha = 0.5
  ) +
  geom_smooth(
    data = df2_plot,
    aes(x = avg_log2FC.x, y = avg_log2FC.y),
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 0.5,
    fill = adjustcolor("grey", alpha.f = 0.3)
  ) +
  stat_cor(
    data = df2_plot,
    aes(x = avg_log2FC.x, y = avg_log2FC.y),
    method = "pearson"
  ) +
  geom_hline(yintercept =  0.1, linetype = "dashed", color = "grey10") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey10") +
  geom_vline(xintercept =  0.1, linetype = "dashed", color = "grey10") +
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey10") +
  labs(
    x = "log2FC (Seeded vs Ctrl)",
    y = "log2FC (PGRN-OE vs PGRN-WT)"
  ) +
  theme_classic()

p_m_nc





#Fig. 5L
#Requried input data: sp (splicing analysis result table; named as SnISOr_Seq_raw_data_related_to_5n in the supplement)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})


df2 <- sp %>%
  filter(
    grepl("^(EN|IN|NeuEarly|Neu_Early)", celltype),
    comparison %in% c("k18_vs_ctrl", "grnoe_vs_k18")
  ) %>%
  dplyr::select(exon, gene, gene_name, celltype, comparison, dpsi, FDR) %>%
  group_by(exon, gene, gene_name, celltype, comparison) %>%
  summarise(
    dpsi = mean(as.numeric(dpsi), na.rm = TRUE),
    FDR  = min(as.numeric(FDR),  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = comparison,
    values_from = c(dpsi, FDR),
    names_sep   = "__"
  ) %>%
  filter(!is.na(dpsi__k18_vs_ctrl), !is.na(dpsi__grnoe_vs_k18))



suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
})

# ---- settings that control the look/logic ----
dpsi_cut <- 0.10
fdr_cut  <- 0.05

stopifnot(all(c("dpsi__k18_vs_ctrl","dpsi__grnoe_vs_k18") %in% names(df2)))

has_fdr <- all(c("FDR__k18_vs_ctrl","FDR__grnoe_vs_k18") %in% names(df2))

df2 <- df2 %>%
  mutate(
    rescued_flag = (dpsi__k18_vs_ctrl * dpsi__grnoe_vs_k18) < 0,
    dot_group    = ifelse(rescued_flag, "Rescued", "Not rescued"),
    sig_both     = if (has_fdr) (!is.na(FDR__k18_vs_ctrl) & !is.na(FDR__grnoe_vs_k18) &
                                   FDR__k18_vs_ctrl < fdr_cut & FDR__grnoe_vs_k18 < fdr_cut) else FALSE
  )

# If your sp already only contains events passing a |dpsi| cutoff, you can remove this.
df2_plot <- df2 %>%
  filter(abs(dpsi__k18_vs_ctrl) >= dpsi_cut | abs(dpsi__grnoe_vs_k18) >= dpsi_cut)

df_sig <- df2_plot %>% filter(sig_both)

# ---------------- base scatter ----------------
p <- ggscatter(
  df2_plot,
  x = "dpsi__k18_vs_ctrl",
  y = "dpsi__grnoe_vs_k18",
  color   = "dot_group",
  palette = c("Not rescued"="grey75", "Rescued"="orchid3"),
  size    = 3.8,
  alpha   = 0.35,
  add     = "reg.line",
  add.params = list(color = "black", size = 0.7, fill = adjustcolor("grey40", alpha.f = 0.35)),
  conf.int   = TRUE,
  # correlation computed on rescued set to match the “rescued” bracket logic
  cor.coef   = TRUE,
  cor.method = "pearson",
  ylab = "ΔPSI (GRN-Iso vs. WT-Seeded condition)",
  xlab = "ΔPSI (WT MG Seeded vs. Non-seeded)",
  repel = FALSE
) +
  geom_hline(yintercept = c(-dpsi_cut, dpsi_cut), linetype = "dashed", color = "grey20") +
  geom_vline(xintercept = c(-dpsi_cut, dpsi_cut), linetype = "dashed", color = "grey20")

# ---------------- yellow overlay: BOTH FDR < 0.05 ----------------
p <- p +
  geom_point(
    data = df_sig,
    aes(x = dpsi__k18_vs_ctrl, y = dpsi__grnoe_vs_k18),
    inherit.aes = FALSE,
    shape = 21, fill = "#FFD700", color = "black",
    size = 3.8, stroke = 0.8, alpha = 0.9
  )

# ---------------- labels + outlines built ONLY from df_sig ----------------
# (n=) counts are from significant rescued set ONLY, like your panel
df_sig_rescued <- df_sig %>% filter(rescued_flag)

gene_counts_sig <- df_sig_rescued %>% count(gene_name, name = "n_exons_sig")
singletons_sig  <- gene_counts_sig %>% filter(n_exons_sig == 1) %>% pull(gene_name)
pairs_sig       <- gene_counts_sig %>% filter(n_exons_sig == 2) %>% pull(gene_name)
multi_gt2_sig   <- gene_counts_sig %>% filter(n_exons_sig >= 3) %>% pull(gene_name)

# pick one representative label point per multi-event gene (most significant by combined -log10 FDR)
label_pts_multi <- df_sig_rescued %>%
  filter(gene_name %in% c(pairs_sig, multi_gt2_sig)) %>%
  group_by(gene_name) %>%
  mutate(weight =
           pmax(-log10(pmax(FDR__k18_vs_ctrl,  1e-300)), 0) +
           pmax(-log10(pmax(FDR__grnoe_vs_k18, 1e-300)), 0)) %>%
  slice_max(weight, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(gene_counts_sig, by = "gene_name") %>%
  mutate(label_txt = paste0(gene_name, " (n=", n_exons_sig, ")"))

label_pts_single <- df_sig_rescued %>%
  filter(gene_name %in% singletons_sig) %>%
  distinct(gene_name, .keep_all = TRUE) %>%
  mutate(label_txt = gene_name)

# segment outline for n=2 genes
pair_segments <- df_sig_rescued %>%
  filter(gene_name %in% pairs_sig) %>%
  group_by(gene_name) %>%
  summarise(
    x1 = dplyr::first(dpsi__k18_vs_ctrl),
    y1 = dplyr::first(dpsi__grnoe_vs_k18),
    x2 = dplyr::last(dpsi__k18_vs_ctrl),
    y2 = dplyr::last(dpsi__grnoe_vs_k18),
    .groups = "drop"
  )
# ellipse outline for n>=3 genes (need ≥3 distinct points)
ellipse_groups <- df_sig_rescued %>%
  filter(gene_name %in% multi_gt2_sig) %>%
  group_by(gene_name) %>%
  mutate(n_pts = n_distinct(interaction(dpsi__k18_vs_ctrl, dpsi__grnoe_vs_k18, drop = TRUE))) %>%
  ungroup() %>%
  filter(n_pts >= 3)

p2 <- p +
  stat_ellipse(
    data = ellipse_groups,
    mapping = aes(x = dpsi__k18_vs_ctrl, y = dpsi__grnoe_vs_k18, group = gene_name),
    type = "norm", level = 0.80, geom = "polygon",
    fill = "grey20", alpha = 0.12,
    color = "black", linewidth = 1.2, linetype = "dotted",
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = pair_segments,
    mapping = aes(x = x1, y = y1, xend = x2, yend = y2),
    inherit.aes = FALSE,
    color = "grey50", linewidth = 3.2, alpha = 0.25, lineend = "round"
  ) +
  geom_segment(
    data = pair_segments,
    mapping = aes(x = x1, y = y1, xend = x2, yend = y2),
    inherit.aes = FALSE,
    color = "black", linewidth = 0.5, lineend = "round"
  ) +
  geom_text_repel(
    data = label_pts_single,
    aes(x = dpsi__k18_vs_ctrl, y = dpsi__grnoe_vs_k18, label = label_txt),
    size = 4.6, max.overlaps = Inf,
    force = 12, force_pull = 1.2,
    box.padding = unit(0.9, "lines"), point.padding = unit(0.5, "lines"),
    min.segment.length = 0, segment.alpha = 0.8, segment.size = 0.35,
    seed = 101
  ) +
  geom_text_repel(
    data = label_pts_multi,
    aes(x = dpsi__k18_vs_ctrl, y = dpsi__grnoe_vs_k18, label = label_txt),
    size = 4.8, max.overlaps = Inf,
    force = 14, force_pull = 1.3,
    box.padding = unit(1.0, "lines"), point.padding = unit(0.55, "lines"),
    min.segment.length = 0, segment.alpha = 0.85, segment.size = 0.35,
    seed = 102
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.22))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.22)))

p2



