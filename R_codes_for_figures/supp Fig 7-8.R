unique(xeno_4genotype_no3r_RELNetc_high_col$Condition)


##seurat object
xeno_5genotype_final

#subset 4R only condition
xeno_4genotype_no3r_RELNetc_high_col <- subset(xeno_5genotype_final, Condition != "3R_WT_Ctrl")


library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(grid)

############################################################
## Use 4R object
############################################################

obj <- xeno_4genotype_no3r_RELNetc_high_col

DefaultAssay(obj) <- "RNA"

############################################################
## Set desired_base levels for control_celltype
############################################################

ct_all <- unique(as.character(obj$control_celltype))
ct_all <- ct_all[!is.na(ct_all)]

EN_levels  <- sort(ct_all[grepl("^Neu_EN_", ct_all)])
IN_levels  <- sort(ct_all[grepl("^Neu_IN_", ct_all)])
NE_levels  <- intersect("Neu_Early", ct_all)
AST_levels <- intersect("AST", ct_all)
MG_levels  <- intersect("MG", ct_all)
NPC_levels <- intersect("NPC_RG", ct_all)

desired_base <- c(EN_levels, IN_levels, NE_levels, AST_levels, MG_levels, NPC_levels)
desired_base <- c(desired_base, setdiff(ct_all, desired_base))

obj$control_celltype <- factor(
  obj$control_celltype,
  levels = desired_base
)

obj$Condition <- factor(
  obj$Condition,
  levels = c("4R_WT_Ctrl", "4R_WT_K18", "4R_P301S_Ctrl", "4R_P301S_K18")
)


############################################################
## Helper function: make cell type percentage UMAP
## Remove Neu_ prefix from neuron subtype labels,
## but keep Neu_Early unchanged
## Larger labels, stronger repel, no bold text
############################################################

library(ggrepel)

make_pct_umap <- function(seu, plot_title = NULL) {
  
  seu$Sample_Name <- droplevels(factor(seu$Sample_Name))
  seu$Condition <- droplevels(factor(seu$Condition))
  seu$control_celltype <- droplevels(factor(seu$control_celltype))
  
  Idents(seu) <- droplevels(seu$control_celltype)
  
  ############################################################
  ## Build percentage labels
  ############################################################
  
  ct <- as.character(seu$control_celltype)
  pct <- prop.table(table(ct)) * 100
  
  base_present <- desired_base[desired_base %in% names(pct)]
  
  label_names <- base_present
  
  label_names <- ifelse(
    grepl("^Neu_", label_names) & label_names != "Neu_Early",
    sub("^Neu_", "", label_names),
    label_names
  )
  
  lab_map <- setNames(
    sprintf("%s (%.1f%%)", label_names, as.numeric(pct[base_present])),
    base_present
  )
  
  vals <- unname(lab_map[ct])
  names(vals) <- colnames(seu)
  
  seu <- AddMetaData(
    seu,
    metadata = vals,
    col.name = "control_celltype_pct"
  )
  
  seu$control_celltype_pct <- factor(
    seu$control_celltype_pct,
    levels = unname(lab_map)
  )
  
  ############################################################
  ## Build label positions using median UMAP location
  ############################################################
  
  umap_df <- as.data.frame(Embeddings(seu, reduction = "umap"))
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
  
  umap_df$control_celltype_pct <- seu$control_celltype_pct
  
  label_df <- umap_df %>%
    group_by(control_celltype_pct) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      .groups = "drop"
    )
  
  ############################################################
  ## Plot UMAP
  ############################################################
  
  p <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "control_celltype_pct",
    label = FALSE,
    repel = FALSE,
    pt.size = 1.2
  ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        label = control_celltype_pct
      ),
      size = 5,
      fontface = "plain",
      color = "black",
      box.padding = 0.6,
      point.padding = 0.4,
      force = 2,
      force_pull = 0.5,
      max.overlaps = Inf,
      min.segment.length = 0
    ) +
    NoAxes() +
    ggtitle(plot_title) +
    theme(
      text = element_text(face = "plain"),
      plot.title = element_text(
        hjust = 0.5,
        face = "plain",
        size = 13
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        face = "plain",
        size = 8
      )
    )
  
  return(p)
}
############################################################
## 1. Entire UMAP: all four conditions included
############################################################

xeno_4R_all <- subset(
  obj,
  subset = Condition %in% c(
    "4R_WT_Ctrl",
    "4R_WT_K18",
    "4R_P301S_Ctrl",
    "4R_P301S_K18"
  )
)

p_f_4R_all <- make_pct_umap(
  xeno_4R_all,
  plot_title = "All 4R conditions"
)

# print(p_f_4R_all)

############################################################
## 2. Individual UMAPs
## Include only these three conditions individually
## Do not make individual UMAP for 4R_WT_Ctrl
############################################################

xeno_4R_WT_K18 <- subset(
  obj,
  subset = Condition == "4R_WT_K18"
)

xeno_4R_P301S_Ctrl <- subset(
  obj,
  subset = Condition == "4R_P301S_Ctrl"
)

xeno_4R_P301S_K18 <- subset(
  obj,
  subset = Condition == "4R_P301S_K18"
)

p_f_4R_WT_K18 <- make_pct_umap(
  xeno_4R_WT_K18,
  plot_title = "4R_WT_K18"
)

p_f_4R_P301S_Ctrl <- make_pct_umap(
  xeno_4R_P301S_Ctrl,
  plot_title = "4R_P301S_Ctrl"
)

p_f_4R_P301S_K18 <- make_pct_umap(
  xeno_4R_P301S_K18,
  plot_title = "4R_P301S_K18"
)

# print(p_f_4R_WT_K18)
# print(p_f_4R_P301S_Ctrl)
# print(p_f_4R_P301S_K18)

############################################################
## 3. Combined panel: one row of 4 UMAPs
## Includes all-condition UMAP plus 3 individual condition UMAPs
############################################################

p_f_4R_all <- p_f_4R_all +
  ggtitle("All 4R conditions")

p_f_4R_WT_K18 <- p_f_4R_WT_K18 +
  ggtitle("4R WT + Tau seeds")

p_f_4R_P301S_Ctrl <- p_f_4R_P301S_Ctrl +
  ggtitle("4R P301S control")

p_f_4R_P301S_K18 <- p_f_4R_P301S_K18 +
  ggtitle("4R P301S + Tau seeds")

p_f_4R_combined <- 
  p_f_4R_all |
  p_f_4R_WT_K18 |
  p_f_4R_P301S_Ctrl |
  p_f_4R_P301S_K18

p_f_4R_combined <- p_f_4R_combined +
  plot_layout(nrow = 1, widths = c(1, 1, 1, 1)) &
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18
    ),
    legend.title = element_blank()
  )

print(p_f_4R_combined)



## ------------------------------------------------------------
## Normalized nuclei percent for ALL subclusters (control_celltype)
## Same pipeline as your “% neuron nuclei / total human nuclei”
## but computed per subcluster instead of “Neuron” aggregate.
##
## Output:
##   1) per-sample long table: percent of each subcluster per Sample_Name
##   2) condition summary table: mean +/- SEM across samples (for bars)
##   3) optional faceted bar+dot plot (one panel per subcluster)
## ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

## 0) Input object
obj <- xeno_5genotype

stopifnot(inherits(obj, "Seurat"))

## 1) Drop 3R at the beginning (keep only 4R conditions)
keep_levels <- c("4R_WT_Ctrl", "4R_WT_K18", "4R_P301S_Ctrl", "4R_P301S_K18")

stopifnot(all(c("Condition", "Sample_Name", "control_celltype") %in% colnames(obj@meta.data)))

present <- intersect(keep_levels, unique(obj$Condition))
stopifnot(length(present) > 0)

obj4r <- subset(obj, subset = Condition %in% present)
md4r  <- obj4r@meta.data

md4r$Condition <- factor(as.character(md4r$Condition), levels = keep_levels)
md4r$Sample_Name <- as.character(md4r$Sample_Name)

## 2) Define subcluster order (use your list; keep only those present)
ct_order <- c(
  "MG","VLMC","NPC_RG","AST","Neu_Early",
  "Neu_IN_Gly","Neu_IN_BCL11B","Neu_IN_RELN","Neu_IN_PROX1",
  "Neu_EN_RELN","Neu_EN_BCL11B","Neu_EN_KCNC2","Neu_EN_PROX1","Neu_EN_PIEZO2"
)
present_ct <- unique(as.character(md4r$control_celltype))
ct_levels  <- intersect(ct_order, present_ct)

md4r$control_celltype <- as.character(md4r$control_celltype)
md4r$control_celltype[is.na(md4r$control_celltype) | md4r$control_celltype == ""] <- "Unknown"
md4r$control_celltype <- factor(md4r$control_celltype, levels = c(ct_levels, setdiff(unique(md4r$control_celltype), ct_levels)))

## 3) Per-sample normalized percent for each subcluster
##    percent_subcluster = 100 * (nuclei in subcluster) / (total human nuclei in that Sample_Name)
per_sample_ct <- md4r %>%
  transmute(
    Sample_Name      = as.character(Sample_Name),
    Condition        = as.character(Condition),
    control_celltype = as.character(control_celltype)
  ) %>%
  group_by(Sample_Name, Condition) %>%
  mutate(total_human_nuclei = n()) %>%
  ungroup() %>%
  count(Sample_Name, Condition, control_celltype, total_human_nuclei, name = "subcluster_nuclei") %>%
  mutate(
    percent_subcluster = 100 * subcluster_nuclei / total_human_nuclei,
    Condition = factor(Condition, levels = keep_levels),
    control_celltype = factor(control_celltype, levels = levels(md4r$control_celltype))
  ) %>%
  arrange(control_celltype, Condition, Sample_Name)

print(per_sample_ct)

write.csv(per_sample_ct,
          file = "normalized_subcluster_percent_per_sample_4Ronly.csv",
          row.names = FALSE)

#Data were replotted in PRISM








#Ex. Fig 7 m-n DEGs_Vulnerble neuron subpopulation_4R_P301S seeded vs 4R_P301S ctrl


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


## ---------------------- 5) Intersections across all four ---------------------
deg_up_intersect_all_full   <- Reduce(intersect, gene_sets_up)
deg_down_intersect_all_full <- Reduce(intersect, gene_sets_down)

cat("Common UP genes (all four):",   length(deg_up_intersect_all_full),   "\n")
cat("Common DOWN genes (all four):", length(deg_down_intersect_all_full), "\n")

.write_gene_list(deg_up_intersect_all_full,   "intersect_all_four_UP_full.csv")
.write_gene_list(deg_down_intersect_all_full, "intersect_all_four_DOWN_full.csv")


## =============================================================================
## 6) GO BP enrichment on intersected UP/DOWN genes
## =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# intersected gene lists you already made
up_genes   <- unique(as.character(deg_up_intersect_all))
down_genes <- unique(as.character(deg_down_intersect_all))

length(up_genes)
length(down_genes)

## ---- Enrichment: BP only ----
do_go_bp <- function(genes, p_cut = 0.05, q_cut = 0.2, simplify_cutoff = 0.6) {
  genes <- unique(as.character(genes))
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
  
  clusterProfiler::simplify(eg, cutoff = simplify_cutoff, by = "p.adjust", select_fun = min)
}

## ---- Format for plotting (same as your old fmt_df, BP only) ----
fmt_bp_df <- function(x, direction, top_n = 10, trunc_width = 55) {
  if (is.null(x)) return(NULL)
  df <- as.data.frame(x)
  if (nrow(df) == 0) return(NULL)
  
  df$term_size <- as.numeric(sub("/.*", "", df$BgRatio))
  df$score     <- -log10(df$p.adjust)
  df$label     <- sprintf("%s / %s", df$Count, df$term_size)
  df$Direction <- direction
  
  df %>%
    arrange(p.adjust, desc(score)) %>%
    slice_head(n = top_n) %>%
    mutate(Description = str_trunc(Description, trunc_width))
}

## ---- Plot (same layout as before; single facet row for BP) ----
plot_bp <- function(df, title, fill_color) {
  if (is.null(df) || nrow(df) == 0)
    return(ggplot() + theme_void() + ggtitle(paste(title, "(no significant BP terms)")))
  
  ggplot(df, aes(x = reorder(Description, score), y = score)) +
    geom_col(fill = fill_color) +
    geom_text(aes(label = label), hjust = 1.02, size = 3, color = "white") +
    coord_flip() +
    labs(x = NULL, y = expression(-log[10]("FDR")), title = title) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          plot.title = element_text(size = 13, face = "bold"))
}

## ---- Run BP enrichment ----
up_bp <- do_go_bp(up_genes)
dn_bp <- do_go_bp(down_genes)

UP_BP   <- fmt_bp_df(up_bp, "Up",   top_n = 15)
DOWN_BP <- fmt_bp_df(dn_bp, "Down", top_n = 15)

p_up   <- plot_bp(UP_BP,   "GO BP (Up genes, intersect)",   "#E84A5F")
p_down <- plot_bp(DOWN_BP, "GO BP (Down genes, intersect)", "#4E79A7")

p_up
p_down

