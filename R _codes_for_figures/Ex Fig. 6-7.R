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

