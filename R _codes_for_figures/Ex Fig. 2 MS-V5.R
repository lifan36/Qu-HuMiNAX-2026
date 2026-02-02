
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



unmatched_cells <- rownames(xeno_wt)[is.na(xeno_wt$velocity_pseudotime)]

# Count and inspect
cat("Number of unmatched cells:", length(unmatched_cells), "\n")
head(unmatched_cells)  # show first few

###cell cycle

library(Seurat)
#cell cycle expression matrix was downloaded from https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "C:/Users/wenhq/Documents/snrna_rds/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes

g2m.genes <- cc.genes$g2m.genes


xeno_wt <- CellCycleScoring(xeno_wt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
# head(xeno_wt[[]])
# 
# RidgePlot(xeno_wt, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

DimPlot(xeno_wt, group.by = "Phase", label = TRUE, repel = TRUE) + 
  ggtitle("Cell Cycle Phase")

unique(xeno_wt$Phase)

table(xeno_wt$Phase)


############################################################
# Extended Data Fig. 2 (xeno_wt)
# Panels a–k in figure order
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

DefaultAssay(xeno_wt) <- "RNA"

Idents(xeno_wt) <- "Sample_Name"
############################################################
# Panels a–g: QC violin plots + scatter plots
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
  
  # Violin plots (Panels a–d)
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
  
  # Scatter plots (Panels e–g)
  p_sc1 <- FeatureScatter(obj, feature1 = "nCount_RNA",   feature2 = "nFeature_RNA")
  p_sc2 <- FeatureScatter(obj, feature1 = "nCount_RNA",   feature2 = "percent.mt")
  p_sc3 <- FeatureScatter(obj, feature1 = "nFeature_RNA", feature2 = "percent.mt")
  p_scatter <- (p_sc1 | p_sc2 | p_sc3)
  
  list(seurat = obj, violin = p_vln, scatter = p_scatter, group_by = group_by)
}

qc_AAA <- make_qc_plots(xeno_wt)

qc_AAA$group_by     # check grouping column
qc_AAA$violin       # Panels a–d
qc_AAA$scatter      # Panels e–g


############################################################
# Panel h: UMAP colored by Condition (3R vs 4R WT Ctrl)
############################################################

#4R dots were a little more transparent to make 3R more visible
DimPlot(xeno_wt, group.by = "Condition", label = F)

tplot <- DimPlot(xeno_wt, group.by = "Condition", repel = TRUE )

tplot[[1]]$layers[[1]]$aes_params$alpha =  ifelse ( xeno_wt@meta.data$Condition == "3R_WT_Ctrl" ,1,0.3)
tplot


############################################################
# Panel i: Cell cycle phase on UMAP
############################################################

# Seurat cell cycle genes
###cell cycle
library(Seurat)

# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "C:/Users/wenhq/Documents/snrna_rds/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

xeno_wt <- CellCycleScoring(xeno_wt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(xeno_wt, group.by = "Phase", label = F) + 
  ggtitle("Cell Cycle Phase")

unique(xeno_wt$Phase)
table(xeno_wt$Phase)




############################################################
# Panel j: UMAP colored by cluster identity
############################################################

# Order by numerical order: 0,1,2,...,8_0,8_1,8_2,...
labs <- as.character(xeno_wt$C8_subcluster)
u    <- unique(labs)

main <- as.integer(sub("_.*", "", u))                              # 0,1,2,...,8,...
subn <- ifelse(grepl("_", u), as.integer(sub(".*_", "", u)), Inf)  # 0,1,2 for 8_*, Inf otherwise
ord  <- order(main, subn)

xeno_wt$C8_subcluster <- factor(labs, levels = u[ord])

p_j <- DimPlot(
  xeno_wt,
  group.by = "C8_subcluster",
  label    = TRUE,
  repel    = TRUE
) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_blank()
  )

print(p_j)

############################################################
# Panel k: DotPlot (markers) with desired cluster order + y labels
############################################################

# Marker set (same as yours)
top_markers <- c(
  "CX3CR1","P2RY12","CSF1R","SALL1",                 # MG
  "COL1A1","COL1A2","DCN",                           # VLMC
  "MKI67","NES","PAX6","SOX1",                       # NPC
  "GFAP","AQP4","ALDH1L1",                           # Astro
  "RBFOX3",
  "GAD1","GAD2",                                     # Inhibitory
  "SLC6A5","GATA3",                                  # Glycinergic
  "RELN","SLC17A6",                                  # Excitatory
  "BCL11B","PROX1","PIEZO2"                          # Subtypes
)

# 1) Desired cluster order
desired <- c(
  "1","21",
  "8_2",
  "23",
  "0","6","8_1",
  "8_0",
  "12","13","17","3","7",
  "2","4",
  "16","9",
  "14",
  "15","18",
  "19",
  "11","20","5",
  "10","22"
)

labs    <- as.character(xeno_wt$C8_subcluster)
present <- desired[desired %in% labs]
extras  <- setdiff(unique(labs), present)
xeno_wt$C8_subcluster <- factor(labs, levels = c(present, extras))

# mapping from cluster -> cell type (your table)
ct_by_cluster <- c(
  "1"="MG","21"="MG",
  "8_2"="VLMC",
  "23"="NPC",
  "0"="AST","6"="AST","8_1"="AST",
  "8_0"="Neu_Early",
  "12"="Neu_IN_Gly","13"="Neu_IN_Gly","17"="Neu_IN_Gly","3"="Neu_IN_Gly","7"="Neu_IN_Gly",
  "2"="Neu_IN_BCL11B","4"="Neu_IN_BCL11B",
  "16"="Neu_IN_PROX1","9"="Neu_IN_PROX1",
  "14"="Neu_IN_RELN",
  "15"="Neu_EN_BCL11B","18"="Neu_EN_BCL11B",
  "19"="Neu_EN_PROX1",
  "11"="Neu_EN_RELN","20"="Neu_EN_RELN","5"="Neu_EN_RELN",
  "10"="Neu_EN_PIEZO2","22"="Neu_EN_PIEZO2"
)

# build y-axis labels using your existing factor order
levs <- levels(xeno_wt$C8_subcluster)
y_labels <- setNames(
  ifelse(
    levs %in% names(ct_by_cluster),
    paste0(levs, " | ", ct_by_cluster[levs]),   # e.g. "1 | MG"
    levs
  ),
  levs
)

DefaultAssay(xeno_wt) <- "RNA"


DotPlot(
  object    = xeno_wt,
  features  = top_markers,
  group.by  = "C8_subcluster",
  scale     = TRUE,
  dot.scale = 6
) +
  scale_color_gradientn(colors = c("darkblue", "gray98", "red")) +
  scale_size(range = c(2, 8)) +
  scale_y_discrete(labels = y_labels) +
  theme_minimal() +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 14),
    plot.title  = element_text(size = 16, face = "bold")
  ) +
  guides(
    size  = guide_legend(override.aes = list(size = 6)),
    color = guide_colorbar(barheight = unit(40, "mm"))
  ) +
  labs(title = "", x = "", y = "")




