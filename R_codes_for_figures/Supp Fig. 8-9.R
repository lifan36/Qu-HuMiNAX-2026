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











##Fig. s9O


#Nat Commun. 2023 Nov 2;14:6801. doi: 10.1038/s41467-023-42626-3

psp_sndrnaseq<-as.data.frame(psp_snrnaseq)

MG_4R301k_vs_4R301ctrl<-as.data.frame(MG_4R301k_vs_4R301ctrl)

AST_4R301k_vs_4R301ctrl<- as.data.frame(AST_4R301k_vs_4R301ctrl)



EN_4R301k_vs_4R301ctrl<-as.data.frame(EN_4R301k_vs_4R301ctrl)

IN_4R301k_vs_4R301ctrl<-as.data.frame(IN_4R301k_vs_4R301ctrl)



MG_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1<-as.data.frame(MG_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1)



IN_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1<-as.data.frame(IN_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1)


EN_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1<-as.data.frame(EN_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1)



AST_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1<-as.data.frame(AST_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1)







library(dplyr)

deg_count_table <- psp_snrnaseq %>%
  filter(
    !is.na(`MAST-logFC`),
    !is.na(`MAST-FDR`),
    abs(`MAST-logFC`) > 0.1,
    `MAST-FDR` < 0.05
  ) %>%
  distinct(
    cell_type,
    GeneName
  ) %>%
  count(
    cell_type,
    name = "number_of_DEGs"
  ) %>%
  arrange(
    desc(number_of_DEGs)
  )

print(deg_count_table)




#from human psp brains
head(psp_sndrnaseq)
unique(psp_snrnaseq$cell_type)

#from Huminax study
head(MG_4R301k_vs_4R301ctrl)
head(EN_4R301k_vs_4R301ctrl)
head(IN_4R301k_vs_4R301ctrl)
head(AST_4R301k_vs_4R301ctrl)

#from mouse ps19 study
head(MG_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1)

head(IN_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1)

head(EN_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1)

head(AST_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1)





library(gprofiler2)
library(dplyr)
library(purrr)
library(readr)

# Mouse PS19 DEG tables
mouse_tables <- list(
  EN  = EN_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1,
  IN  = IN_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1,
  AST = AST_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1,
  MG  = MG_E3WT_PS19_vs_E3WT_DEGs_log2FC_0_1
)

# Collect all unique mouse gene symbols
mouse_genes <- mouse_tables %>%
  map(~ as.character(.x$X)) %>%
  unlist(use.names = FALSE) %>%
  trimws() %>%
  unique()

mouse_genes <- mouse_genes[
  !is.na(mouse_genes) &
    mouse_genes != ""
]

# Mouse symbols -> human ortholog symbols
orthologs <- gprofiler2::gorth(
  query = mouse_genes,
  source_organism = "mmusculus",
  target_organism = "hsapiens",
  mthreshold = Inf,
  filter_na = TRUE
)

# Simple mapping table
mouse_to_human <- orthologs %>%
  transmute(
    mouse_gene = trimws(as.character(input)),
    human_gene = toupper(
      trimws(as.character(ortholog_name))
    )
  ) %>%
  filter(
    !is.na(mouse_gene),
    mouse_gene != "",
    !is.na(human_gene),
    human_gene != ""
  ) %>%
  distinct()

head(mouse_to_human)

cat(
  "Input mouse symbols:",
  length(mouse_genes),
  "\nMouse symbols with human orthologs:",
  n_distinct(mouse_to_human$mouse_gene),
  "\nUnique human ortholog symbols:",
  n_distinct(mouse_to_human$human_gene),
  "\n"
)



# ============================================================
# Prepare each mouse DEG table using mapped human symbols
# ============================================================

prepare_mouse_gorth <- function(
    df,
    cell_group,
    ortholog_map,
    fdr_cutoff = 0.05
) {
  
  mapped <- df %>%
    transmute(
      mouse_gene = trimws(as.character(X)),
      mouse_logFC = as.numeric(avg_log2FC),
      mouse_pvalue = as.numeric(p_val),
      mouse_FDR = as.numeric(p_val_adj)
    ) %>%
    filter(
      !is.na(mouse_gene),
      mouse_gene != "",
      is.finite(mouse_logFC),
      !is.na(mouse_FDR),
      mouse_FDR < fdr_cutoff
    ) %>%
    inner_join(
      ortholog_map,
      by = "mouse_gene"
    )
  
  # Resolve the mapped results at the human-gene level.
  # One mouse gene may map to multiple human genes,
  # and multiple mouse genes may map to one human gene.
  mapped %>%
    group_by(human_gene) %>%
    summarise(
      any_up = any(
        mouse_logFC > 0,
        na.rm = TRUE
      ),
      
      any_down = any(
        mouse_logFC < 0,
        na.rm = TRUE
      ),
      
      mouse_genes = paste(
        sort(unique(mouse_gene)),
        collapse = ";"
      ),
      
      n_mouse_genes = n_distinct(mouse_gene),
      
      min_mouse_FDR = min(
        mouse_FDR,
        na.rm = TRUE
      ),
      
      max_positive_logFC = ifelse(
        any(mouse_logFC > 0),
        max(
          mouse_logFC[mouse_logFC > 0],
          na.rm = TRUE
        ),
        NA_real_
      ),
      
      min_negative_logFC = ifelse(
        any(mouse_logFC < 0),
        min(
          mouse_logFC[mouse_logFC < 0],
          na.rm = TRUE
        ),
        NA_real_
      ),
      
      .groups = "drop"
    ) %>%
    mutate(
      direction = case_when(
        any_up & !any_down ~ "Up",
        !any_up & any_down ~ "Down",
        any_up & any_down ~ "Mixed",
        TRUE ~ "Zero"
      ),
      
      gene = human_gene,
      dataset = "PS19_mouse",
      cell_group = cell_group
    )
}


mouse_prepared <- imap(
  mouse_tables,
  ~ prepare_mouse_gorth(
    df = .x,
    cell_group = .y,
    ortholog_map = mouse_to_human,
    fdr_cutoff = 0.05
  )
)


# Check the actual mapped PS19 counts
mouse_mapping_count_summary <- imap_dfr(
  mouse_prepared,
  function(df, cell_group) {
    
    tibble(
      cell_group = cell_group,
      mapped_human_orthologs = n_distinct(df$gene),
      mapped_up = n_distinct(
        df$gene[df$direction == "Up"]
      ),
      mapped_down = n_distinct(
        df$gene[df$direction == "Down"]
      ),
      mapped_mixed = n_distinct(
        df$gene[df$direction == "Mixed"]
      )
    )
  }
)

print(mouse_mapping_count_summary)

head(mouse_prepared)


# ============================================================
# Convert PSP gene symbols:
#
# Original PSP symbol
#       ↓
# GENCODE v32 Ensembl gene ID
#       ↓
# GENCODE v50 gene symbol
#
# Final R object:
#   psp_snrnaseq_updated
#
# Final CSV:
#   PSP_gene_name_conversion/psp_snrnaseq_updated.csv
#
# Important:
#   - Every original PSP DEG row is retained.
#   - No rows are deduplicated.
#   - Ambiguous or unmapped symbols retain their original name.
#   - GeneName contains the updated symbol whenever conversion
#     is uniquely resolved.
# ============================================================

library(rtracklayer)
library(dplyr)
library(readr)
library(tibble)

# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------

output_dir <- "PSP_gene_name_conversion"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

old_release <- 32
new_release <- 50

old_gtf_file <- file.path(
  output_dir,
  paste0(
    "gencode.v",
    old_release,
    ".annotation.gtf.gz"
  )
)

new_gtf_file <- file.path(
  output_dir,
  paste0(
    "gencode.v",
    new_release,
    ".annotation.gtf.gz"
  )
)

old_gtf_url <- paste0(
  "https://ftp.ebi.ac.uk/pub/databases/gencode/",
  "Gencode_human/release_",
  old_release,
  "/gencode.v",
  old_release,
  ".annotation.gtf.gz"
)

new_gtf_url <- paste0(
  "https://ftp.ebi.ac.uk/pub/databases/gencode/",
  "Gencode_human/release_",
  new_release,
  "/gencode.v",
  new_release,
  ".annotation.gtf.gz"
)

# ------------------------------------------------------------
# 2. Download annotations
# ------------------------------------------------------------

if (!file.exists(old_gtf_file)) {
  download.file(
    url = old_gtf_url,
    destfile = old_gtf_file,
    mode = "wb"
  )
}

if (!file.exists(new_gtf_file)) {
  download.file(
    url = new_gtf_url,
    destfile = new_gtf_file,
    mode = "wb"
  )
}

if (
  !file.exists(old_gtf_file) ||
  file.info(old_gtf_file)$size == 0
) {
  stop(
    "GENCODE v",
    old_release,
    " annotation was not downloaded correctly."
  )
}

if (
  !file.exists(new_gtf_file) ||
  file.info(new_gtf_file)$size == 0
) {
  stop(
    "GENCODE v",
    new_release,
    " annotation was not downloaded correctly."
  )
}

# ------------------------------------------------------------
# 3. Import GENCODE annotations
# ------------------------------------------------------------

old_gtf <- rtracklayer::import(
  old_gtf_file
)

new_gtf <- rtracklayer::import(
  new_gtf_file
)

# ------------------------------------------------------------
# 4. GENCODE v32:
#    original symbol -> Ensembl gene ID
# ------------------------------------------------------------

old_gene_map <- as.data.frame(old_gtf) %>%
  dplyr::filter(
    type == "gene"
  ) %>%
  dplyr::transmute(
    GeneName_original = trimws(
      as.character(gene_name)
    ),
    
    GeneName_original_upper = toupper(
      trimws(
        as.character(gene_name)
      )
    ),
    
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      as.character(gene_id)
    )
  ) %>%
  dplyr::filter(
    !is.na(GeneName_original),
    GeneName_original != "",
    !is.na(ensembl_gene_id),
    ensembl_gene_id != ""
  ) %>%
  dplyr::distinct()

# ------------------------------------------------------------
# 5. GENCODE v50:
#    Ensembl gene ID -> updated symbol
# ------------------------------------------------------------

new_gene_map <- as.data.frame(new_gtf) %>%
  dplyr::filter(
    type == "gene"
  ) %>%
  dplyr::transmute(
    ensembl_gene_id = sub(
      "\\..*$",
      "",
      as.character(gene_id)
    ),
    
    GeneName_current = trimws(
      as.character(gene_name)
    )
  ) %>%
  dplyr::filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != "",
    !is.na(GeneName_current),
    GeneName_current != ""
  ) %>%
  dplyr::distinct()

# Summarize current annotation so there is exactly one row
# for each Ensembl gene ID.

new_gene_map_summary <- new_gene_map %>%
  dplyr::group_by(
    ensembl_gene_id
  ) %>%
  dplyr::summarise(
    n_current_symbols = dplyr::n_distinct(
      GeneName_current
    ),
    
    current_symbols = paste(
      sort(
        unique(GeneName_current)
      ),
      collapse = ";"
    ),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 6. Unique PSP symbols
# ------------------------------------------------------------

all_psp_symbols <- psp_snrnaseq %>%
  dplyr::transmute(
    GeneName_original = trimws(
      as.character(GeneName)
    ),
    
    GeneName_original_upper = toupper(
      trimws(
        as.character(GeneName)
      )
    )
  ) %>%
  dplyr::filter(
    !is.na(GeneName_original),
    GeneName_original != ""
  ) %>%
  dplyr::distinct()

# ------------------------------------------------------------
# 7. Exact v32 symbol mapping
# ------------------------------------------------------------

old_exact_mapping <- old_gene_map %>%
  dplyr::group_by(
    GeneName_original
  ) %>%
  dplyr::summarise(
    n_exact_ensembl_ids = dplyr::n_distinct(
      ensembl_gene_id
    ),
    
    exact_ensembl_ids = paste(
      sort(
        unique(ensembl_gene_id)
      ),
      collapse = ";"
    ),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 8. Case-insensitive v32 rescue
#
# Used only when there is no exact match and the uppercase
# symbol identifies exactly one Ensembl gene.
# ------------------------------------------------------------

old_case_insensitive_mapping <- old_gene_map %>%
  dplyr::group_by(
    GeneName_original_upper
  ) %>%
  dplyr::summarise(
    n_case_insensitive_ensembl_ids =
      dplyr::n_distinct(
        ensembl_gene_id
      ),
    
    case_insensitive_ensembl_ids = paste(
      sort(
        unique(ensembl_gene_id)
      ),
      collapse = ";"
    ),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 9. Select one old Ensembl ID when mapping is unambiguous
# ------------------------------------------------------------

psp_old_symbol_mapping <- all_psp_symbols %>%
  dplyr::left_join(
    old_exact_mapping,
    by = "GeneName_original"
  ) %>%
  dplyr::left_join(
    old_case_insensitive_mapping,
    by = "GeneName_original_upper"
  ) %>%
  dplyr::mutate(
    n_exact_ensembl_ids = dplyr::coalesce(
      n_exact_ensembl_ids,
      0L
    ),
    
    n_case_insensitive_ensembl_ids =
      dplyr::coalesce(
        n_case_insensitive_ensembl_ids,
        0L
      ),
    
    ensembl_gene_id = dplyr::case_when(
      n_exact_ensembl_ids == 1 ~
        exact_ensembl_ids,
      
      n_exact_ensembl_ids == 0 &
        n_case_insensitive_ensembl_ids == 1 ~
        case_insensitive_ensembl_ids,
      
      TRUE ~
        NA_character_
    ),
    
    old_mapping_status = dplyr::case_when(
      n_exact_ensembl_ids == 1 ~
        "Exact GENCODE v32 symbol match",
      
      n_exact_ensembl_ids > 1 ~
        "Ambiguous exact GENCODE v32 mapping",
      
      n_exact_ensembl_ids == 0 &
        n_case_insensitive_ensembl_ids == 1 ~
        "Case-insensitive GENCODE v32 rescue",
      
      n_exact_ensembl_ids == 0 &
        n_case_insensitive_ensembl_ids > 1 ~
        "Ambiguous case-insensitive GENCODE v32 mapping",
      
      TRUE ~
        "Not found in GENCODE v32"
    )
  ) %>%
  dplyr::select(
    GeneName_original,
    ensembl_gene_id,
    old_mapping_status
  )

# ------------------------------------------------------------
# 10. Convert v32 Ensembl IDs to v50 symbols
# ------------------------------------------------------------

psp_symbol_conversion <- psp_old_symbol_mapping %>%
  dplyr::left_join(
    new_gene_map_summary,
    by = "ensembl_gene_id"
  ) %>%
  dplyr::mutate(
    n_current_symbols = dplyr::coalesce(
      n_current_symbols,
      0L
    ),
    
    GeneName_current = dplyr::case_when(
      !is.na(ensembl_gene_id) &
        n_current_symbols == 1 ~
        current_symbols,
      
      TRUE ~
        GeneName_original
    ),
    
    mapping_status = dplyr::case_when(
      old_mapping_status ==
        "Not found in GENCODE v32" ~
        "Not found in GENCODE v32",
      
      grepl(
        "^Ambiguous",
        old_mapping_status
      ) ~
        old_mapping_status,
      
      !is.na(ensembl_gene_id) &
        n_current_symbols == 0 ~
        paste0(
          "Ensembl gene absent from GENCODE v",
          new_release
        ),
      
      !is.na(ensembl_gene_id) &
        n_current_symbols > 1 ~
        paste0(
          "Ambiguous GENCODE v",
          new_release,
          " symbol"
        ),
      
      !is.na(ensembl_gene_id) &
        n_current_symbols == 1 &
        GeneName_original == GeneName_current ~
        "Already current",
      
      !is.na(ensembl_gene_id) &
        n_current_symbols == 1 &
        GeneName_original != GeneName_current ~
        "Updated through Ensembl",
      
      TRUE ~
        "Unresolved"
    ),
    
    symbol_changed =
      mapping_status == "Updated through Ensembl"
  ) %>%
  dplyr::select(
    GeneName_original,
    ensembl_gene_id,
    GeneName_current,
    old_mapping_status,
    mapping_status,
    symbol_changed
  )

# Confirm exactly one mapping row exists per original symbol.

if (
  anyDuplicated(
    psp_symbol_conversion$GeneName_original
  ) > 0
) {
  stop(
    "The conversion table contains duplicate original symbols."
  )
}

# ------------------------------------------------------------
# 11. Create the final updated PSP table
#
# No rows are removed or deduplicated.
# GeneName becomes the updated symbol.
# ------------------------------------------------------------

original_row_count <- nrow(
  psp_snrnaseq
)

psp_snrnaseq_updated <- psp_snrnaseq %>%
  dplyr::mutate(
    PSP_original_row = dplyr::row_number(),
    
    GeneName_original = trimws(
      as.character(GeneName)
    )
  ) %>%
  dplyr::left_join(
    psp_symbol_conversion,
    by = "GeneName_original"
  ) %>%
  dplyr::mutate(
    GeneName_current = dplyr::coalesce(
      GeneName_current,
      GeneName_original
    ),
    
    mapping_status = dplyr::coalesce(
      mapping_status,
      "Original symbol retained"
    ),
    
    symbol_changed = dplyr::coalesce(
      symbol_changed,
      FALSE
    ),
    
    GeneName = GeneName_current
  ) %>%
  dplyr::arrange(
    PSP_original_row
  ) %>%
  dplyr::select(
    -PSP_original_row
  )

# Confirm that the conversion did not add or remove DEG rows.

if (
  nrow(psp_snrnaseq_updated) !=
  original_row_count
) {
  stop(
    "Row count changed during gene-symbol conversion."
  )
}

# ------------------------------------------------------------
# 12. Print a simple conversion summary
# ------------------------------------------------------------

conversion_summary <- psp_symbol_conversion %>%
  dplyr::count(
    mapping_status,
    name = "n_unique_symbols"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_unique_symbols)
  )

print(
  conversion_summary,
  n = Inf
)

cat(
  "\nOriginal PSP rows:",
  nrow(psp_snrnaseq),
  "\nUpdated PSP rows:",
  nrow(psp_snrnaseq_updated),
  "\nChanged symbols:",
  sum(
    psp_symbol_conversion$symbol_changed,
    na.rm = TRUE
  ),
  "\n"
)



head(mouse_prepared)
head(psp_snrnaseq_updated)

# ============================================================
# DOWNSTREAM ANALYSIS ONLY
#
# Directional DEG comparison:
#   HuMiNAX vs human PSP vs mouse PS19
#
# Required existing objects:
#   1. psp_snrnaseq_updated
#   2. mouse_prepared
#   3. EN_4R301k_vs_4R301ctrl
#   4. IN_4R301k_vs_4R301ctrl
#   5. AST_4R301k_vs_4R301ctrl
#   6. MG_4R301k_vs_4R301ctrl
#
# Cell-type matching:
#   EN  = HuMiNAX EN  vs PSP exc_neu vs PS19 EN
#   IN  = HuMiNAX IN  vs PSP inh_neu vs PS19 IN
#   AST = HuMiNAX AST vs PSP ast     vs PS19 AST
#   MG  = HuMiNAX MG  vs PSP mic     vs PS19 MG
#
# Significance:
#   HuMiNAX: p_val_adj < 0.05
#   PSP:     MAST-FDR < 0.05
#   PS19:    already filtered when mouse_prepared was created
#
# No additional log2FC cutoff is applied.
#
# Outputs:
#   - Eight directional Venn diagrams
#   - One combined Venn figure
#   - PSP directional recapitulation heatmap
#   - Directional gene-set and overlap CSV files
# ============================================================


# ------------------------------------------------------------
# 1. Packages and settings
# ------------------------------------------------------------

library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(readr)
library(ggplot2)
library(ggVennDiagram)
library(patchwork)
library(scales)

output_dir <- "HuMiNAX_PSP_PS19_directional_analysis"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

fdr_cutoff <- 0.05

cell_groups <- c(
  "EN",
  "IN",
  "AST",
  "MG"
)

dataset_order <- c(
  "HuMiNAX",
  "Human PSP",
  "Mouse PS19"
)


# ------------------------------------------------------------
# 2. Confirm that the required objects exist
# ------------------------------------------------------------

required_object_names <- c(
  "psp_snrnaseq_updated",
  "mouse_prepared",
  "EN_4R301k_vs_4R301ctrl",
  "IN_4R301k_vs_4R301ctrl",
  "AST_4R301k_vs_4R301ctrl",
  "MG_4R301k_vs_4R301ctrl"
)

missing_objects <- required_object_names[
  !vapply(
    required_object_names,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0) {
  stop(
    "These required objects are missing: ",
    paste(
      missing_objects,
      collapse = ", "
    )
  )
}

if (
  !is.list(mouse_prepared) ||
  !all(cell_groups %in% names(mouse_prepared))
) {
  stop(
    "mouse_prepared must be a named list containing ",
    "EN, IN, AST, and MG."
  )
}


# ------------------------------------------------------------
# 3. Connect the existing input objects
# ------------------------------------------------------------

huminax_tables <- list(
  EN  = EN_4R301k_vs_4R301ctrl,
  IN  = IN_4R301k_vs_4R301ctrl,
  AST = AST_4R301k_vs_4R301ctrl,
  MG  = MG_4R301k_vs_4R301ctrl
)

psp_cell_types <- c(
  EN  = "exc_neu",
  IN  = "inh_neu",
  AST = "ast",
  MG  = "mic"
)


# ------------------------------------------------------------
# 4. Confirm the required columns
# ------------------------------------------------------------

check_columns <- function(
    df,
    required_columns,
    object_name
) {
  
  missing_columns <- setdiff(
    required_columns,
    colnames(df)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      object_name,
      " is missing these columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
}

purrr::iwalk(
  huminax_tables,
  function(df, cell_group) {
    check_columns(
      df,
      c(
        "gene",
        "avg_log2FC",
        "p_val_adj"
      ),
      paste0(
        "HuMiNAX ",
        cell_group,
        " table"
      )
    )
  }
)

check_columns(
  psp_snrnaseq_updated,
  c(
    "GeneName",
    "cell_type",
    "MAST-logFC",
    "MAST-FDR"
  ),
  "psp_snrnaseq_updated"
)

purrr::iwalk(
  mouse_prepared[cell_groups],
  function(df, cell_group) {
    check_columns(
      df,
      c(
        "gene",
        "direction"
      ),
      paste0(
        "mouse_prepared$",
        cell_group
      )
    )
  }
)


# ------------------------------------------------------------
# 5. Gene-symbol cleaning
#
# All comparisons are performed with uppercase human symbols.
# ------------------------------------------------------------

clean_gene <- function(x) {
  
  output <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  output[
    is.na(output) |
      output == ""
  ] <- NA_character_
  
  output
}


# ------------------------------------------------------------
# 6. Prepare HuMiNAX or PSP DEG rows
#
# This function:
#   1. Keeps FDR < 0.05.
#   2. Does not impose another logFC cutoff.
#   3. Assigns Up or Down from the logFC sign.
#   4. Marks a gene Mixed if significant rows have both signs.
# ------------------------------------------------------------

summarize_deg_directions <- function(
    df,
    gene_column,
    logfc_column,
    fdr_column
) {
  
  df %>%
    transmute(
      gene = clean_gene(
        .data[[gene_column]]
      ),
      
      logFC = suppressWarnings(
        as.numeric(
          .data[[logfc_column]]
        )
      ),
      
      FDR = suppressWarnings(
        as.numeric(
          .data[[fdr_column]]
        )
      )
    ) %>%
    filter(
      !is.na(gene),
      is.finite(logFC),
      is.finite(FDR),
      FDR < fdr_cutoff
    ) %>%
    group_by(gene) %>%
    summarise(
      any_up = any(
        logFC > 0,
        na.rm = TRUE
      ),
      
      any_down = any(
        logFC < 0,
        na.rm = TRUE
      ),
      
      minimum_FDR = min(
        FDR,
        na.rm = TRUE
      ),
      
      n_significant_rows = n(),
      
      .groups = "drop"
    ) %>%
    mutate(
      direction = case_when(
        any_up & !any_down ~ "Up",
        !any_up & any_down ~ "Down",
        any_up & any_down ~ "Mixed",
        TRUE ~ "Zero"
      )
    )
}


# ------------------------------------------------------------
# 7. Prepare HuMiNAX directional genes
# ------------------------------------------------------------

huminax_direction_summary <- purrr::imap(
  huminax_tables,
  function(df, cell_group) {
    
    summarize_deg_directions(
      df = df,
      gene_column = "gene",
      logfc_column = "avg_log2FC",
      fdr_column = "p_val_adj"
    )
  }
)


# ------------------------------------------------------------
# 8. Prepare PSP directional genes
#
# PSP contains multiple clusters within each broad cell type.
# A gene significant in opposite directions across clusters is
# classified as Mixed and excluded from directional comparisons.
# ------------------------------------------------------------

psp_direction_summary <- purrr::imap(
  psp_cell_types,
  function(psp_cell_type, cell_group) {
    
    psp_subset <- psp_snrnaseq_updated %>%
      filter(
        as.character(cell_type) ==
          psp_cell_type
      )
    
    summarize_deg_directions(
      df = psp_subset,
      gene_column = "GeneName",
      logfc_column = "MAST-logFC",
      fdr_column = "MAST-FDR"
    )
  }
)


# ------------------------------------------------------------
# 9. Use the already processed PS19 outcome
#
# mouse_prepared already contains:
#   gene      = mapped human ortholog
#   direction = Up, Down, or Mixed
#
# This section does not remap mouse genes.
# ------------------------------------------------------------

summarize_mouse_outcome <- function(df) {
  
  df %>%
    transmute(
      gene = clean_gene(gene),
      
      original_direction = case_when(
        tolower(
          trimws(
            as.character(direction)
          )
        ) == "up" ~ "Up",
        
        tolower(
          trimws(
            as.character(direction)
          )
        ) == "down" ~ "Down",
        
        tolower(
          trimws(
            as.character(direction)
          )
        ) == "mixed" ~ "Mixed",
        
        TRUE ~ "Zero"
      )
    ) %>%
    filter(
      !is.na(gene)
    ) %>%
    group_by(gene) %>%
    summarise(
      any_up = any(
        original_direction %in%
          c(
            "Up",
            "Mixed"
          )
      ),
      
      any_down = any(
        original_direction %in%
          c(
            "Down",
            "Mixed"
          )
      ),
      
      .groups = "drop"
    ) %>%
    mutate(
      direction = case_when(
        any_up & !any_down ~ "Up",
        !any_up & any_down ~ "Down",
        any_up & any_down ~ "Mixed",
        TRUE ~ "Zero"
      )
    )
}

mouse_direction_summary <- purrr::imap(
  mouse_prepared[cell_groups],
  function(df, cell_group) {
    summarize_mouse_outcome(df)
  }
)


# ------------------------------------------------------------
# 10. Combine the three processed datasets
# ------------------------------------------------------------

all_direction_summaries <- bind_rows(
  
  purrr::imap_dfr(
    huminax_direction_summary,
    function(df, cell_group) {
      
      df %>%
        mutate(
          dataset = "HuMiNAX",
          cell_group = cell_group
        )
    }
  ),
  
  purrr::imap_dfr(
    psp_direction_summary,
    function(df, cell_group) {
      
      df %>%
        mutate(
          dataset = "Human PSP",
          cell_group = cell_group
        )
    }
  ),
  
  purrr::imap_dfr(
    mouse_direction_summary,
    function(df, cell_group) {
      
      df %>%
        mutate(
          dataset = "Mouse PS19",
          cell_group = cell_group
        )
    }
  )
) %>%
  select(
    cell_group,
    dataset,
    gene,
    direction,
    everything()
  )


# Genes used in the directional plots

directional_gene_table <- all_direction_summaries %>%
  filter(
    direction %in%
      c(
        "Up",
        "Down"
      )
  ) %>%
  distinct(
    cell_group,
    dataset,
    direction,
    gene
  )


# Mixed-direction genes excluded from the plots

mixed_direction_gene_table <- all_direction_summaries %>%
  filter(
    direction == "Mixed"
  ) %>%
  distinct(
    cell_group,
    dataset,
    gene,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# 11. Directional DEG counts
# ------------------------------------------------------------

directional_deg_counts <- directional_gene_table %>%
  count(
    cell_group,
    dataset,
    direction,
    name = "n_directional_DEGs"
  ) %>%
  complete(
    cell_group = cell_groups,
    dataset = dataset_order,
    direction = c(
      "Up",
      "Down"
    ),
    fill = list(
      n_directional_DEGs = 0
    )
  ) %>%
  arrange(
    factor(
      cell_group,
      levels = cell_groups
    ),
    factor(
      dataset,
      levels = dataset_order
    ),
    factor(
      direction,
      levels = c(
        "Up",
        "Down"
      )
    )
  )

print(
  directional_deg_counts,
  n = Inf
)


# ------------------------------------------------------------
# 12. Helper to retrieve one directional gene set
# ------------------------------------------------------------

get_gene_set <- function(
    cell_group,
    dataset_name,
    direction_name
) {
  
  directional_gene_table %>%
    filter(
      cell_group == !!cell_group,
      dataset == dataset_name,
      direction == direction_name
    ) %>%
    pull(gene) %>%
    unique() %>%
    sort()
}


# ------------------------------------------------------------
# 13. Create one directional three-way Venn diagram
# ------------------------------------------------------------

make_directional_venn <- function(
    cell_group,
    direction_name
) {
  
  huminax_genes <- get_gene_set(
    cell_group,
    "HuMiNAX",
    direction_name
  )
  
  psp_genes <- get_gene_set(
    cell_group,
    "Human PSP",
    direction_name
  )
  
  ps19_genes <- get_gene_set(
    cell_group,
    "Mouse PS19",
    direction_name
  )
  
  venn_input <- list(
    "HuMiNAX" = huminax_genes,
    "Human PSP" = psp_genes,
    "Mouse PS19" = ps19_genes
  )
  
  high_fill <- if (
    direction_name == "Up"
  ) {
    "#B33A3A"
  } else {
    "#3F6F9F"
  }
  
  direction_title <- if (
    direction_name == "Up"
  ) {
    "upregulated DEGs"
  } else {
    "downregulated DEGs"
  }
  
  ggVennDiagram(
    venn_input,
    label_alpha = 0,
    edge_size = 0.8,
    set_size = 4
  ) +
    scale_fill_gradient(
      low = "white",
      high = high_fill
    ) +
    labs(
      title = paste0(
        cell_group,
        ": ",
        direction_title
      ),
      subtitle = paste0(
        "Same-direction significant genes; ",
        "mixed-direction genes excluded"
      )
    ) +
    theme_void() +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 13
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 8.5
      ),
      legend.position = "none",
      plot.margin = margin(
        10,
        25,
        10,
        25
      )
    )
}


# ------------------------------------------------------------
# 14. Generate and save the eight Venn diagrams
#
# Order:
#   Top row    = EN, IN, AST, MG up
#   Bottom row = EN, IN, AST, MG down
# ------------- -----------------------------------------------

venn_plots <- list()

for (
  direction_name in c(
    "Up",
    "Down"
  )
) {
  
  for (
    cell_group in cell_groups
  ) {
    
    plot_name <- paste(
      cell_group,
      direction_name,
      sep = "_"
    )
    
    current_plot <- make_directional_venn(
      cell_group = cell_group,
      direction_name = direction_name
    )
    
    venn_plots[[plot_name]] <- current_plot
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0(
          cell_group,
          "_",
          tolower(direction_name),
          "_directional_venn.png"
        )
      ),
      plot = current_plot,
      width = 6.5,
      height = 5.5,
      dpi = 300,
      bg = "white"
    )
  }
}

combined_directional_venn <- wrap_plots(
  plotlist = venn_plots,
  ncol = 4,
  nrow = 2
) +
  plot_annotation(
    title = paste0(
      "Directional DEG overlap: ",
      "HuMiNAX vs human PSP vs mouse PS19"
    ),
    subtitle = paste0(
      "Top row: upregulated DEGs; ",
      "bottom row: downregulated DEGs"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 18
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 11
      )
    )
  )

print(
  combined_directional_venn
)

# ggsave(
#   filename = file.path(
#     output_dir,
#     "combined_directional_venn.png"
#   ),
#   plot = combined_directional_venn,
#   width = 18,
#   height = 10,
#   dpi = 300,
#   bg = "white"
# )
# 
# ggsave(
#   filename = file.path(
#     output_dir,
#     "combined_directional_venn.pdf"
#   ),
#   plot = combined_directional_venn,
#   width = 18,
#   height = 10,
#   bg = "white"
# )


# ------------------------------------------------------------
# 15. Export the genes assigned to each Venn region
# ------------------------------------------------------------

make_venn_membership_table <- function(
    cell_group,
    direction_name
) {
  
  huminax_genes <- get_gene_set(
    cell_group,
    "HuMiNAX",
    direction_name
  )
  
  psp_genes <- get_gene_set(
    cell_group,
    "Human PSP",
    direction_name
  )
  
  ps19_genes <- get_gene_set(
    cell_group,
    "Mouse PS19",
    direction_name
  )
  
  all_genes <- sort(
    unique(
      c(
        huminax_genes,
        psp_genes,
        ps19_genes
      )
    )
  )
  
  tibble(
    cell_group = rep(
      cell_group,
      length(all_genes)
    ),
    
    direction = rep(
      direction_name,
      length(all_genes)
    ),
    
    gene = all_genes,
    
    in_HuMiNAX =
      all_genes %in%
      huminax_genes,
    
    in_human_PSP =
      all_genes %in%
      psp_genes,
    
    in_mouse_PS19 =
      all_genes %in%
      ps19_genes
  ) %>%
    mutate(
      venn_region = case_when(
        in_HuMiNAX &
          in_human_PSP &
          in_mouse_PS19 ~
          "HuMiNAX + human PSP + mouse PS19",
        
        in_HuMiNAX &
          in_human_PSP ~
          "HuMiNAX + human PSP",
        
        in_HuMiNAX &
          in_mouse_PS19 ~
          "HuMiNAX + mouse PS19",
        
        in_human_PSP &
          in_mouse_PS19 ~
          "Human PSP + mouse PS19",
        
        in_HuMiNAX ~
          "HuMiNAX only",
        
        in_human_PSP ~
          "Human PSP only",
        
        in_mouse_PS19 ~
          "Mouse PS19 only"
      )
    )
}

directional_venn_membership <- purrr::map_dfr(
  cell_groups,
  function(cell_group) {
    
    purrr::map_dfr(
      c(
        "Up",
        "Down"
      ),
      function(direction_name) {
        
        make_venn_membership_table(
          cell_group = cell_group,
          direction_name = direction_name
        )
      }
    )
  }
)


# ------------------------------------------------------------
# 16. Calculate directional recapitulation of human PSP
#
# For every tile:
#
#   fraction recapitulated =
#       model genes overlapping PSP in the same direction
#       -------------------------------------------------
#       total PSP genes in that direction
#
# Comparisons:
#   HuMiNAX vs human PSP
#   Mouse PS19 vs human PSP
# ------------------------------------------------------------

model_comparisons <- tibble(
  model_dataset = c(
    "HuMiNAX",
    "Mouse PS19"
  ),
  
  model_label = c(
    "HuMiNAX",
    "PS19"
  )
)

psp_directional_heatmap_df <- purrr::map_dfr(
  cell_groups,
  function(cell_group) {
    
    purrr::map_dfr(
      c(
        "Up",
        "Down"
      ),
      function(direction_name) {
        
        psp_genes <- get_gene_set(
          cell_group,
          "Human PSP",
          direction_name
        )
        
        purrr::map_dfr(
          seq_len(
            nrow(model_comparisons)
          ),
          function(model_index) {
            
            model_dataset <-
              model_comparisons$model_dataset[
                model_index
              ]
            
            model_label <-
              model_comparisons$model_label[
                model_index
              ]
            
            model_genes <- get_gene_set(
              cell_group,
              model_dataset,
              direction_name
            )
            
            overlapping_genes <- intersect(
              psp_genes,
              model_genes
            )
            
            n_psp_genes <- length(
              psp_genes
            )
            
            n_overlap <- length(
              overlapping_genes
            )
            
            fraction_recapitulated <- if (
              n_psp_genes == 0
            ) {
              NA_real_
            } else {
              n_overlap /
                n_psp_genes
            }
            
            tibble(
              cell_group = cell_group,
              model_dataset = model_dataset,
              direction = direction_name,
              
              comparison = paste0(
                model_label,
                ": PSP ",
                tolower(direction_name)
              ),
              
              n_model_directional_DEGs =
                length(model_genes),
              
              n_PSP_directional_DEGs =
                n_psp_genes,
              
              n_same_direction_overlap =
                n_overlap,
              
              fraction_recapitulated =
                fraction_recapitulated,
              
              label = ifelse(
                is.na(
                  fraction_recapitulated
                ),
                "NA",
                paste0(
                  n_overlap,
                  "/",
                  n_psp_genes,
                  "\n",
                  scales::percent(
                    fraction_recapitulated,
                    accuracy = 1
                  )
                )
              )
            )
          }
        )
      }
    )
  }
)


# Set the desired heatmap order

comparison_levels <- c(
  "HuMiNAX: PSP up",
  "HuMiNAX: PSP down",
  "PS19: PSP up",
  "PS19: PSP down"
)

psp_directional_heatmap_df <-
  psp_directional_heatmap_df %>%
  mutate(
    comparison = factor(
      comparison,
      levels = comparison_levels
    ),
    
    cell_group = factor(
      cell_group,
      levels = rev(
        cell_groups
      )
    )
  )




# ------------------------------------------------------------
# 17. Make the PSP directional recapitulation heatmap
#
# Upregulated PSP DEGs:
#   white -> dark red
#
# Downregulated PSP DEGs:
#   white -> dark blue
# ------------------------------------------------------------

# Install once if needed:
# install.packages("ggnewscale")

library(ggnewscale)

# Separate the upregulated and downregulated heatmap data

psp_heatmap_up <- psp_directional_heatmap_df %>%
  filter(
    direction == "Up"
  )

psp_heatmap_down <- psp_directional_heatmap_df %>%
  filter(
    direction == "Down"
  )


# Determine a separate maximum for each color scale

up_heatmap_max <- max(
  psp_heatmap_up$fraction_recapitulated,
  na.rm = TRUE
)

down_heatmap_max <- max(
  psp_heatmap_down$fraction_recapitulated,
  na.rm = TRUE
)

if (
  !is.finite(up_heatmap_max) ||
  up_heatmap_max <= 0
) {
  up_heatmap_max <- 1
}

if (
  !is.finite(down_heatmap_max) ||
  down_heatmap_max <= 0
) {
  down_heatmap_max <- 1
}


# Make the heatmap

psp_directional_heatmap <- ggplot(
  psp_directional_heatmap_df,
  aes(
    x = comparison,
    y = cell_group
  )
) +
  
  # ----------------------------------------------------------
# Upregulated PSP DEG recapitulation
# ----------------------------------------------------------

geom_tile(
  data = psp_heatmap_up,
  aes(
    fill = fraction_recapitulated
  ),
  color = "white",
  linewidth = 1
) +
  
  scale_fill_gradient(
    name = paste0(
      "Upregulated PSP DEGs\n",
      "recapitulated"
    ),
    low = "#FFFFF0",
    high = adjustcolor("#8B1A1A", alpha.f = 0.9),
    labels = scales::percent_format(
      accuracy = 1
    ),
    limits = c(
      0,
      up_heatmap_max
    ),
    oob = scales::squish,
    na.value = "grey90"
  ) +
  
  # Start a second independent fill scale
  
  ggnewscale::new_scale_fill() +
  
  # ----------------------------------------------------------
# Downregulated PSP DEG recapitulation
# ----------------------------------------------------------

geom_tile(
  data = psp_heatmap_down,
  aes(
    fill = fraction_recapitulated
  ),
  color = "#FFFFF0",
  linewidth = 1
) +
  
  scale_fill_gradient(
    name = paste0(
      "Downregulated PSP DEGs\n",
      "recapitulated"
    ),
    low = "#FFFFF0",
    high = adjustcolor("#08306B", alpha.f = 0.8),
    labels = scales::percent_format(
      accuracy = 1
    ),
    limits = c(
      0,
      down_heatmap_max
    ),
    oob = scales::squish,
    na.value = "grey70"
  ) +
  
  # ----------------------------------------------------------
# Numbers and percentages printed on all tiles
# ----------------------------------------------------------

geom_text(
  aes(
    label = label
  ),
  size = 4.0,
) +
  
  labs(
    title = "Recapitulation of human PSP directional DEGs",
    subtitle = paste0(
      "Each tile shows same-direction overlap / ",
      "total PSP up- or downregulated DEGs"
    ),
    x = NULL,
    y = NULL
  ) +
  
  theme_bw(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      size = 15
    ),
    
    plot.subtitle = element_text(
      size = 11
    ),
    
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    
    axis.text.y = element_text(
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 10
    )
  )


# Display the heatmap

print(
  psp_directional_heatmap
)










