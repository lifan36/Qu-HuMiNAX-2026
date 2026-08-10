## =============================================================================
## MAJIQ v3 splicing analysis (Ctrl vs KD)
## Generates:
##   (1) Ctrl PSI (median) vs KD PSI (median) scatter with highlighted events
##   (2) GO BP enrichment barplot for genes from significant splicing events
## Logic preserved from your original script.
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

## =============================================================================
## 1) Input and basic filtering (BH FDR)
## =============================================================================
HNRNPK_MAJIQv3_HET_USE <- as.data.frame(HNRNPK_MAJIQv3_HET_USE)
df <- HNRNPK_MAJIQv3_HET_USE

# BH-adjusted FDR (kept exactly as you did)
df$FDR <- p.adjust(df$`ttest-approximate_pvalue`, method = "BH")
sig <- subset(df, FDR < 0.05)

## =============================================================================
## 2) Extract median PSI columns, compute delta PSI, define highlights
## =============================================================================

# Find the median-quantile PSI columns (robust to separators in column names)
ctrl_col <- grep("Ctrl.*raw.*psi.*quantile.*0\\.500", names(df), value = TRUE)[1]
kd_col   <- grep("KD.*raw.*psi.*quantile.*0\\.500",   names(df), value = TRUE)[1]
stopifnot(!is.na(ctrl_col), !is.na(kd_col))

# Add PSI and delta PSI (KD - Ctrl)
df$CtrlPsi  <- as.numeric(df[[ctrl_col]])
df$KDPsi    <- as.numeric(df[[kd_col]])
df$deltaPsi <- df$KDPsi - df$CtrlPsi

# Highlight events: p < 0.05 AND |dPSI| > 0.1 (using raw p, same as your code)
df$highlight <- !is.na(df$`ttest-approximate_pvalue`) &
  df$`ttest-approximate_pvalue` < 0.05 &
  !is.na(df$deltaPsi) &
  abs(df$deltaPsi) > 0.1

## =============================================================================
## 3) Figure 6e: PSI scatter (with highlighted events)
## =============================================================================
dot_size <- 2

p_scatter <- ggplot(df, aes(x = CtrlPsi, y = KDPsi)) +
  geom_point(color = "grey70", alpha = 0.6, size = dot_size, na.rm = TRUE) +
  geom_point(
    data = subset(df, highlight),
    color = "palevioletred3",
    alpha = 0.8,
    size = dot_size,
    na.rm = TRUE
  ) +
  geom_abline(slope = 1, intercept = 0, color = "palevioletred3", linetype = "dashed") +
  theme_classic() +
  labs(
    x = "Ctrl PSI (median)",
    y = "KD PSI (median)",
    title = "Ctrl vs KD splicing (MAJIQ v3)",
    subtitle = "p < 0.05 and |dPSI| > 0.1"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0))

p_scatter

## =============================================================================
## 4) Export significant splicing events (p < 0.05 AND |dPSI| > 0.1)
## =============================================================================
sig_events <- df %>%
  filter(
    !is.na(`ttest-approximate_pvalue`),
    `ttest-approximate_pvalue` < 0.05,
    !is.na(deltaPsi),
    abs(deltaPsi) > 0.1
  ) %>%
  arrange(`ttest-approximate_pvalue`, desc(abs(deltaPsi)))

# Quick check
nrow(sig_events)

# Write CSV (same filename as your code)
write_csv(sig_events, "HNRNPK_MAJIQ_significant_events_p0.05_dpsi0.1.csv")


# =============================================================================
# DESeq2 differential expression (KD vs Ctrl)
# Input: count_gene_names (data.frame) with a first column "gene_name"
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

# ----------------------------
# 1) User parameters
# ----------------------------
FC    <- 1.4   # kept as in your script (not used downstream)
FDR   <- 0.05  # kept as in your script (not used downstream)
alpha <- 0.05  # used in DESeq2::results()

# ----------------------------
# 2) Prepare count matrix
# ----------------------------
# Starting from converted/filtered counts where gene IDs were converted to Ensembl,
# and a "gene_name" column exists.
raw_counts <- as.data.frame(count_gene_names)

raw_counts$gene_name <- as.character(raw_counts$gene_name)
raw_counts$gene_name <- make.unique(raw_counts$gene_name)
rownames(raw_counts) <- raw_counts$gene_name

# Remove the gene_name column (counts only)
raw_counts <- raw_counts[, -1, drop = FALSE]

# Optional quick check
head(raw_counts)

# ----------------------------
# 3) Sample metadata (colData)
# ----------------------------
col_data <- data.frame(
  sample = c("Ctrl_1", "Ctrl_2", "Ctrl_3", "KD_1", "KD_2", "KD_3"),
  groups = c("Ctrl",   "Ctrl",   "Ctrl",   "KD",   "KD",   "KD"),
  row.names = c("Ctrl_1", "Ctrl_2", "Ctrl_3", "KD_1", "KD_2", "KD_3")
)

# Ensure count matrix columns match colData rows (same samples, same order)
stopifnot(all(colnames(raw_counts) %in% rownames(col_data)))
col_data <- col_data[colnames(raw_counts), , drop = FALSE]

# ----------------------------
# 4) Run DESeq2
# ----------------------------
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData   = col_data,
  design    = ~ groups
)

dds <- DESeq2::DESeq(dds)

# ----------------------------
# 5) Extract results (KD - Ctrl)
# ----------------------------
res <- DESeq2::results(
  dds,
  contrast = c("groups", "KD", "Ctrl"),
  independentFiltering = TRUE,
  alpha = alpha
)

res_df <- as.data.frame(res)
res_clean <- res_df[complete.cases(res_df), ]


head(res_clean)  # top (lowest row order; not necessarily top by padj)
tail(res_clean)

nrow(res_clean)

# ----------------------------
# 6) Save outputs
# ----------------------------
out_file <- "DESeq2_lgFC0.5_moreinfo.csv"
write.csv(
  as.data.frame(res_clean),
  file = out_file,
  row.names = TRUE
)





####GO anlaysis ex. Fig. 12 a-b
library(clusterProfiler)
library(org.Hs.eg.db)
library(genekitr)

# Step 1: Subset upregulated and downregulated genes
res_up <- res_clean[res_clean$log2FoldChange > 0.5 & res_clean$padj < 0.05, ]
res_down <- res_clean[res_clean$log2FoldChange < -0.5 & res_clean$padj < 0.05, ]

# Step 2: Run GO enrichment
GO_results_up <- enrichGO(
  gene = row.names(res_up),
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

GO_results_down <- enrichGO(
  gene = row.names(res_down),
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)


# Step 3: Advanced comparative visualization with genekitr

up_go <- as.data.frame(GO_results_up) %>%
  filter(ONTOLOGY == "BP") %>%
  arrange(p.adjust) %>%
  slice_head(n = 20)

up_go2 <- up_go %>%
  arrange(desc(Count)) %>%
  mutate(Description = factor(Description, levels = rev(Description)))



down_go <- as.data.frame(GO_results_down) %>%
  dplyr::filter(ONTOLOGY == "BP") %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 20)

down_go2 <- down_go %>%
  dplyr::arrange(Count) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description)))



plotEnrichAdv(
  up_go2, down_go2,
  plot_type = "two",
  term_metric = "Count",
  stats_metric = "p.adjust"
)







## =============================================================================
## Pathway enrichment + term–gene network (pathfindR) and in vivo vs in vitro
## pathway-gene overlap table
## =============================================================================

suppressPackageStartupMessages({
  library(pathfindR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggraph)
  library(grid)
})

# ----------------------------
# 1) Parameters (same logic)
# ----------------------------
padj_cut <- 0.05
lfc_cut  <- 0.5
pathways_focus <- c("Apoptosis", "Necroptosis", "Ferroptosis", "Cellular senescence")

# ----------------------------
# 2) Prepare input for pathfindR from DESeq2 results (same logic)
# ----------------------------
stopifnot(all(c("log2FoldChange", "padj") %in% colnames(res_clean)))

pathway_input <- data.frame(
  Gene.symbol = rownames(res_clean),
  logFC       = res_clean$log2FoldChange,
  adj.P.Val   = res_clean$padj,
  stringsAsFactors = FALSE
)

RA_input <- subset(pathway_input, adj.P.Val < padj_cut & abs(logFC) > lfc_cut)

# ----------------------------
# 3) Run pathfindR (default = KEGG; same logic)
# ----------------------------
RA_output <- run_pathfindR(RA_input)

# Save pathfindR output
write.csv(RA_output, "HNBK_res_clean_pathways.csv", row.names = FALSE)

# # ----------------------------
# # 4) Quick visualizations (same logic, just ordered execution)
# # ----------------------------
# # Network for all enriched terms
# term_gene_graph(result_df = RA_output, use_description = TRUE)
# 
# # Ranking by Fold_Enrichment (as you did)
# RA_ranked <- RA_output[order(RA_output$Fold_Enrichment), ]
# enrichment_chart(RA_ranked, top_terms = 25) +
#   scale_color_gradient(low = "blue", high = "red")
# 
# # Clustered chart for all terms
# RA_clustered <- cluster_enriched_terms(RA_output)
# enrichment_chart(RA_clustered, plot_by_cluster = TRUE) +
#   scale_color_gradient(low = "blue", high = "red")
# 
# # Clustered chart for top 100 by Fold_Enrichment
# RA_top100 <- RA_output[order(RA_output$Fold_Enrichment), ][1:100, ]
# RA_clustered_top100 <- cluster_enriched_terms(RA_top100)
# enrichment_chart(RA_clustered_top100, plot_by_cluster = TRUE) +
#   scale_color_gradient(low = "blue", high = "red")

# ----------------------------
# 5) Focused term–gene graph for selected pathways (same logic)
# ----------------------------
selected_res_path1 <- RA_output[RA_output$Term_Description %in% pathways_focus, ]

if (nrow(selected_res_path1) > 0) {
  p <- term_gene_graph(
    result_df       = selected_res_path1,
    use_description = TRUE,
    node_colors     = c("#E5D7BF", "red", "blue")   # enriched term / up / down
  )
  
  # Match your labeling customization
  options(ggrepel.max.overlaps = Inf)
  
  p <- p + ggraph::geom_node_text(
    aes(label = name),
    repel = TRUE,
    max.overlaps = Inf,
    size = 5,
    point.padding = grid::unit(0.6, "lines"),
    box.padding   = grid::unit(0.45, "lines"),
    force         = 12,
    force_pull    = 0.6,
    min.segment.length = 0,
    segment.alpha = 0.5,
    segment.size  = 0.3
  )
  
  print(p)
}


#Tabel 1 
#Common pathways of #Apoptosis", "Necroptosis", "Ferroptosis", "Cellular senescence" were selected from the RA_output and saved as csv named bulk
#pathways of #Apoptosis", "Necroptosis", "Ferroptosis", "Cellular senescence" from Fig.3 were also selected and saved as csv named snRNAseq_NEU_path

# =============================================================================
# Part B: Compare pathway gene lists between in vivo (bulk) and in vitro (snRNAseq)
# =============================================================================

# ----------------------------
# 6) Inputs (same objects, just standardized names)
# ----------------------------
bulk_death <- as.data.frame(bulk)
snRNAseq_Neu_death <- as.data.frame(snRNAseq_NEU_path)

stopifnot(all(c("Term_Description", "Up_regulated", "Down_regulated") %in% colnames(bulk_death)))
stopifnot(all(c("Term_Description", "Up_regulated", "Down_regulated") %in% colnames(snRNAseq_Neu_death)))

# ----------------------------
# 7) Helper: make pathway -> gene vector list (same logic)
# ----------------------------
get_gene_list <- function(df, gene_col, pathways_keep) {
  df %>%
    mutate(
      Term_Description = trimws(as.character(Term_Description)),
      gene_str = as.character(.data[[gene_col]])
    ) %>%
    filter(Term_Description %in% pathways_keep) %>%
    separate_rows(gene_str, sep = "\\s*,\\s*") %>%
    filter(!is.na(gene_str), gene_str != "") %>%
    transmute(Term_Description, Gene.symbol = gene_str) %>%
    distinct() %>%
    group_by(Term_Description) %>%
    summarise(genes = list(sort(unique(Gene.symbol))), .groups = "drop") %>%
    { setNames(.$genes, .$Term_Description) }
}

# ----------------------------
# 8) Pathway gene lists per dataset (Up and Down)
# ----------------------------
bulk_up_list   <- get_gene_list(bulk_death, "Up_regulated", pathways_focus)
bulk_down_list <- get_gene_list(bulk_death, "Down_regulated", pathways_focus)

sn_up_list     <- get_gene_list(snRNAseq_Neu_death, "Up_regulated", pathways_focus)
sn_down_list   <- get_gene_list(snRNAseq_Neu_death, "Down_regulated", pathways_focus)

# Ensure all pathways exist (even if missing)
for (pw in pathways_focus) {
  if (is.null(bulk_up_list[[pw]]))   bulk_up_list[[pw]]   <- character(0)
  if (is.null(bulk_down_list[[pw]])) bulk_down_list[[pw]] <- character(0)
  if (is.null(sn_up_list[[pw]]))     sn_up_list[[pw]]     <- character(0)
  if (is.null(sn_down_list[[pw]]))   sn_down_list[[pw]]   <- character(0)
}

# ----------------------------
# 9) Common genes (intersection) and summary table
# ----------------------------
common_up_list <- setNames(
  lapply(pathways_focus, function(pw) intersect(bulk_up_list[[pw]], sn_up_list[[pw]])),
  pathways_focus
)

common_down_list <- setNames(
  lapply(pathways_focus, function(pw) intersect(bulk_down_list[[pw]], sn_down_list[[pw]])),
  pathways_focus
)

common_by_pathway_tbl <- tibble(
  Term_Description = pathways_focus,
  n_common_up   = sapply(pathways_focus, function(pw) length(common_up_list[[pw]])),
  common_up     = sapply(pathways_focus, function(pw) paste(common_up_list[[pw]], collapse = ", ")),
  n_common_down = sapply(pathways_focus, function(pw) length(common_down_list[[pw]])),
  common_down   = sapply(pathways_focus, function(pw) paste(common_down_list[[pw]], collapse = ", "))
)

print(common_by_pathway_tbl)








##requrie output from Fig 6 HNRNPK_MAJIQ_significant_events_p0.05_dpsi0.1.csv
sig_genes <- as.data.frame(HNRNPK_MAJIQ_significant_events_p0_05_dpsi0_1)
# Genes for pathway analysis
stopifnot("gene_name" %in% colnames(sig_events))
sig_genes <- unique(sig_events$gene_name)

## =============================================================================
##  GO enrichment (BP) for genes from significant splicing events
## =============================================================================
sig_genes <- unique(as.character(sig_genes))
sig_genes <- sig_genes[!is.na(sig_genes) & sig_genes != ""]

GO_sig <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)


## =============================================================================
## Custom GO barplot (Top 20 by Count; fill = -log10(p.adjust))
## =============================================================================
go_df <- as.data.frame(GO_sig) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

p_go <- ggplot(go_df, aes(x = Count, y = Description, fill = -log10(p.adjust))) +
  geom_col() +
  scale_fill_gradient(low = "#1E90FF", high = "#FF6A6A", name = "-log10(p.adjust)") +
  theme_classic() +
  labs(
    title = "Pathways of significantly spliced genes",
    x = "Count",
    y = NULL
  )

p_go



