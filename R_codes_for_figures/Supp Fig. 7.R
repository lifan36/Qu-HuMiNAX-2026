############################################################
# Extended Data QC (two Seurat objects)
# Panels a–g: xeno_5genotype_4Ronly
# Panels h–n: xeno_wKI_KIonly (GRN KI-only subset)
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})


############################################################
# QC function (used for ALL objects)
############################################################

make_qc_plots <- function(obj, assay = "RNA") {
  stopifnot(inherits(obj, "Seurat"))
  if (assay %in% Assays(obj)) DefaultAssay(obj) <- assay
  
  # auto pick group.by
  md <- obj@meta.data
  group_by <- NULL
  if ("Sample_Name" %in% colnames(md)) {
    group_by <- "Sample_Name"
  } else if ("orig.ident" %in% colnames(md)) {
    group_by <- "orig.ident"
  }
  
  feats <- rownames(obj)
  
  # percent.mt: human MT- and mouse mt-
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    mt_features <- grep("^(MT-|mt-)", feats, value = TRUE)
    if (length(mt_features) > 0) {
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mt_features)
    } else {
      obj[["percent.mt"]] <- 0
      message("No MT-/mt- genes found in rownames(obj); percent.mt set to 0.")
    }
  }
  
  # percent.ribo: RPL/RPS
  if (!"percent.ribo" %in% colnames(obj@meta.data)) {
    ribo_features <- grep("^RP[SL]", feats, value = TRUE)
    if (length(ribo_features) > 0) {
      obj[["percent.ribo"]] <- PercentageFeatureSet(obj, features = ribo_features)
    } else {
      obj[["percent.ribo"]] <- 0
      message("No RPL/RPS genes found; percent.ribo set to 0.")
    }
  }
  
  qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
  qc_features <- qc_features[qc_features %in% colnames(obj@meta.data)]
  
  # Violin plots
  if (!is.null(group_by)) {
    p_vln <- VlnPlot(
      obj,
      features = qc_features,
      group.by = group_by,
      pt.size  = 0,
      ncol     = 4
    ) + NoLegend()
  } else {
    p_vln <- VlnPlot(
      obj,
      features = qc_features,
      pt.size  = 0,
      ncol     = 4
    ) + NoLegend()
  }
  
  # Scatter plots
  p_sc1 <- FeatureScatter(obj, feature1 = "nCount_RNA",   feature2 = "nFeature_RNA")
  p_sc2 <- FeatureScatter(obj, feature1 = "nCount_RNA",   feature2 = "percent.mt")
  p_sc3 <- FeatureScatter(obj, feature1 = "nFeature_RNA", feature2 = "percent.mt")
  p_scatter <- (p_sc1 | p_sc2 | p_sc3)
  
  list(seurat = obj, violin = p_vln, scatter = p_scatter, group_by = group_by)
}


############################################################
# PART 1: QC for xeno_5genotype_4Ronly (Panels a–g)
############################################################

DefaultAssay(xeno_5genotype_4Ronly) <- "RNA"

# ---- KEEP correct levels: Sample_Name ordering ----
# 1) Extract group name (remove _1, _2, etc.)
xeno_5genotype_4Ronly$Sample_Group <- sub("_[0-9]+$", "", xeno_5genotype_4Ronly$Sample_Name)

# 2) Extract replicate number
xeno_5genotype_4Ronly$Replicate <- as.integer(sub(".*_([0-9]+)$", "\\1", xeno_5genotype_4Ronly$Sample_Name))

# 3) Define the biological order you want
group_order <- c("4R_WT_Tau_seeds", "4R_P301S_Ctrl", "4R_P301S_Tau_seeds")

# 4) Build ordered Sample_Name levels (WT reps first, then P301S ctrl, then P301S seeds)
meta_tmp <- xeno_5genotype_4Ronly@meta.data

ordered_levels <- meta_tmp$Sample_Name[order(
  match(meta_tmp$Sample_Group, group_order),
  meta_tmp$Replicate
)]
ordered_levels <- unique(ordered_levels)

# 5) Apply ordering
xeno_5genotype_4Ronly$Sample_Name <- factor(
  xeno_5genotype_4Ronly$Sample_Name,
  levels = ordered_levels
)

# Check order
levels(xeno_5genotype_4Ronly$Sample_Name)

# Set identities (matches your workflow)
Idents(xeno_5genotype_4Ronly) <- "Sample_Name"


# ---- RUN QC ----
qc_xeno_5genotype_4Ronly <- make_qc_plots(xeno_5genotype_4Ronly)

qc_xeno_5genotype_4Ronly$group_by
qc_xeno_5genotype_4Ronly$violin   # Panels a–d
qc_xeno_5genotype_4Ronly$scatter  # Panels e–g




############################################################
# PART 2: QC for GRN (xeno_wKI_KIonly) (Panels h–n)
############################################################

DefaultAssay(xeno_wKI_final_label) <- "RNA"

# Confirm sample names
unique(xeno_wKI_final_label$Sample_Name)

# =========================================================
# 1) Subset xeno_wKI_final_label to KI + ISOKI ONLY
# =========================================================
xeno_wKI_KIonly <- subset(
  xeno_wKI_final_label,
  subset = grepl("^(KI_|ISOKI_)", Sample_Name)
)

unique(xeno_wKI_KIonly$Sample_Name)

# =========================================================
# 2) Rename K18 -> Tau_seeds (Sample_Name + Condition)
# =========================================================
xeno_wKI_KIonly$Sample_Name <- gsub("_K18_", "_Tau_seeds_", xeno_wKI_KIonly$Sample_Name)

if ("Condition" %in% colnames(xeno_wKI_KIonly@meta.data)) {
  xeno_wKI_KIonly$Condition <- gsub("_K18", "_Tau_seeds", xeno_wKI_KIonly$Condition)
  xeno_wKI_KIonly$Condition <- gsub("K18", "Tau_seeds", xeno_wKI_KIonly$Condition)
}

unique(xeno_wKI_KIonly$Sample_Name)
if ("Condition" %in% colnames(xeno_wKI_KIonly@meta.data)) unique(xeno_wKI_KIonly$Condition)

# =========================================================
# 3) Order Sample_Name nicely: iso first, then GRN_2
# =========================================================
xeno_wKI_KIonly$Sample_Name <- as.character(xeno_wKI_KIonly$Sample_Name)

# Rename sample names (ISOKI -> HPC_GRN_iso, KI -> HPC_GRN_2)
xeno_wKI_KIonly$Sample_Name <- gsub("^ISOKI_4R_P301S_Tau_seeds_", "HPC_GRN_iso_", xeno_wKI_KIonly$Sample_Name)
xeno_wKI_KIonly$Sample_Name <- gsub("^KI_4R_P301S_Tau_seeds_",    "HPC_GRN_2_",   xeno_wKI_KIonly$Sample_Name)

# Extract group + replicate number
xeno_wKI_KIonly$Sample_Group <- sub("_[0-9]+$", "", xeno_wKI_KIonly$Sample_Name)
xeno_wKI_KIonly$Replicate    <- as.integer(sub(".*_([0-9]+)$", "\\1", xeno_wKI_KIonly$Sample_Name))

# Desired order: iso first, then GRN_2
group_order <- c("HPC_GRN_iso", "HPC_GRN_2")

meta_tmp <- xeno_wKI_KIonly@meta.data

ordered_levels <- meta_tmp$Sample_Name[order(
  match(meta_tmp$Sample_Group, group_order),
  meta_tmp$Replicate
)]
ordered_levels <- unique(ordered_levels)

# Apply ordering + set identities
xeno_wKI_KIonly$Sample_Name <- factor(xeno_wKI_KIonly$Sample_Name, levels = ordered_levels)
Idents(xeno_wKI_KIonly) <- "Sample_Name"

# Check final order
levels(xeno_wKI_KIonly$Sample_Name)

# =========================================================
# 4) Run QC
# =========================================================
DefaultAssay(xeno_wKI_KIonly) <- "RNA"

qc_xeno_wKI_KIonly <- make_qc_plots(xeno_wKI_KIonly)

qc_xeno_wKI_KIonly$group_by
qc_xeno_wKI_KIonly$violin   # Panels h–k
qc_xeno_wKI_KIonly$scatter  # Panels l–n

table(xeno_wKI_KIonly$Condition)



