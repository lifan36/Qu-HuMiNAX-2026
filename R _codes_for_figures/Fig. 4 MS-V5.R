
#MG DEG table: MG_4R301k_vs_4R301ctrl
MG_4R301k_vs_4R301ctrl <- as.data.frame(MG_4R301k_vs_4R301ctrl)


## =============================================================================
## Figure 4: MG seeded vs Ctrl
## Directional pathway enrichment (pathfindR) + AD risk and IFN DotPlots
## Reorganized for publication (same logic, no extra analyses)
## =============================================================================

## -----------------------------------------------------------------------------
## 0) Libraries
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(pathfindR)
  library(ggplot2)
  library(patchwork)
  library(Seurat)
  library(pheatmap)
})

## -----------------------------------------------------------------------------
## 1) Directional KEGG enrichment with pathfindR (Up vs Down)
##    Input: MG_4R301k_vs_4R301ctrl (DEG table with columns: gene, avg_log2FC, p_val_adj)
## -----------------------------------------------------------------------------


# Make pathfindR formatted table
MG_tbl <- MG_4R301k_vs_4R301ctrl %>%
  transmute(
    Gene.symbol = as.character(gene),
    logFC       = as.numeric(avg_log2FC),
    adj.P.Val   = as.numeric(p_val_adj)
  ) %>%
  filter(!is.na(Gene.symbol), !is.na(logFC), !is.na(adj.P.Val)) %>%
  arrange(adj.P.Val) %>%
  distinct(Gene.symbol, .keep_all = TRUE)

# Split into up vs down DEG sets
padj_cutoff   <- 0.05
logfc_abs_cut <- 0.1

RA_input_up <- MG_tbl %>% filter(adj.P.Val < padj_cutoff,  logFC >  logfc_abs_cut)
RA_input_dn <- MG_tbl %>% filter(adj.P.Val < padj_cutoff,  logFC < -logfc_abs_cut)

# Guardrails
if (nrow(RA_input_up) == 0 && nrow(RA_input_dn) == 0) {
  stop("No genes passed thresholds for either direction.")
}

# Run pathfindR separately (same seed, same defaults)
set.seed(301)
res_up <- if (nrow(RA_input_up) > 0) {
  out <- run_pathfindR(RA_input_up)   # defaults (KEGG, p<0.05)
  if (nrow(out) > 0) out$Direction <- "Up"
  out
} else {
  NULL
}

set.seed(301)
res_dn <- if (nrow(RA_input_dn) > 0) {
  out <- run_pathfindR(RA_input_dn)
  if (nrow(out) > 0) out$Direction <- "Down"
  out
} else {
  NULL
}

RA_bidir <- bind_rows(res_up, res_dn)

# Plot: one chart per direction
if (!is.null(res_up) && nrow(res_up) > 0) {
  res_up <- res_up[order(res_up$Fold_Enrichment, decreasing = TRUE), ]
  p_up <- enrichment_chart(res_up, top_terms = 20) +
    scale_color_gradient(low = "blue", high = "red") +
    ggtitle("Up-regulated DEG pathways") +
    aes(y = reorder(Term_Description, Fold_Enrichment))
  print(p_up)
}

if (!is.null(res_dn) && nrow(res_dn) > 0) {
  res_dn <- res_dn[order(res_dn$Fold_Enrichment, decreasing = TRUE), ]
  p_dn <- enrichment_chart(res_dn, top_terms = 20) +
    scale_color_gradient(low = "blue", high = "red") +
    ggtitle("Down-regulated DEG pathways") +
    aes(y = reorder(Term_Description, Fold_Enrichment))
  print(p_dn)
}



# Save combined directional results
out_dir <- "G:/Other computers/My Laptop/sequencing analysis/MS/fig 3 microglia"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  RA_bidir,
  file.path(out_dir, "MG_4R301k_vs_4R301ctrl_pathways_bidir.csv"),
  row.names = FALSE
)

# Combine up/down panels side by side and save
if (exists("p_up") && exists("p_dn")) {
  p_combined <- p_up | p_dn
  print(p_combined)
  
  ggsave(
    "MG_4R301k_vs_4R301ctrl_pathways_up_down.png",
    p_combined, width = 14, height = 6, dpi = 300
  )
}

## -----------------------------------------------------------------------------
## Fig 4c AD risk genes intersection and DotPlot
## -----------------------------------------------------------------------------

#AD risk genes were downloaded from this paper: PMID: 36907103 PMCID: PMC10024184 DOI: 10.1016/j.ebiom.2023.104511
#Supplementary Table S1, column Gene, exported as .csv and loaded as AD_risk


# -----------------------------
# 1) Build AD-risk intersect DEG table
# -----------------------------
ad_risk <- as.data.frame(AD_risk)
AD_risk_genes <- unique(as.character(ad_risk$Gene))

deg <- as.data.frame(MG_4R301k_vs_4R301ctrl)
deg$gene <- as.character(deg$gene)

# Keep only DEGs whose gene symbol is in AD-risk list
deg_ad <- deg[deg$gene %in% AD_risk_genes, , drop = FALSE]

# Direction label (as in your code)
deg_ad$Direction <- ifelse(deg_ad$avg_log2FC >= 0, "Up", "Down")

# Final ordering used in your original workflow
deg_ad <- deg_ad[order(deg_ad$avg_log2FC, decreasing = TRUE), , drop = FALSE]
deg_ad <- deg_ad[order(deg_ad$p_val_adj, -abs(deg_ad$avg_log2FC)), , drop = FALSE]


deg_ad <- deg_ad[order(deg_ad$p_val_adj, -deg_ad$avg_log2FC), , drop = FALSE]



deg_ad <- deg_ad[order(
  deg_ad$Direction != "Up",   # Up first
  deg_ad$p_val_adj,           # then smaller adj p
  ifelse(deg_ad$Direction == "Up",
         -deg_ad$avg_log2FC,  # Up: bigger positive first
         deg_ad$avg_log2FC)  # Down: more negative first
), , drop = FALSE]



# Write out
out_csv <- "MG_4R301k_vs_4R301ctrl_intersect_AD_risk_genes.csv"
write.csv(deg_ad, file = out_csv, row.names = FALSE)

# -----------------------------
# 2) Prepare Seurat object and subset MG
# -----------------------------
conditions_oi <- c("4R_P301S_Ctrl", "4R_P301S_K18")

obj <- subset(xeno_5genotype, subset = Condition %in% conditions_oi)
DefaultAssay(obj) <- "RNA"

mg_obj <- subset(obj, subset = control_superclass == "MG")

# -----------------------------
# 3) DotPlot for AD-risk intersect genes
# -----------------------------
mg2 <- subset(mg_obj, subset = Condition %in% conditions_oi)
mg2$Condition <- factor(mg2$Condition, levels = conditions_oi)
Idents(mg2) <- "Condition"

genes_use <- unique(as.character(deg_ad$gene))
genes_use <- genes_use[!is.na(genes_use) & nzchar(genes_use)]
genes_use <- genes_use[genes_use %in% rownames(mg2)]


# genes_use is already your AD-risk intersect list (filtered to rownames)
gene_labels <- setNames(genes_use, genes_use)

p_ad_top <- DotPlot(
  object    = mg2,                     # or mg_obj if you prefer; mg2 is already 2 conditions
  features  = genes_use,
  group.by  = "Condition",
  cols      = c("lightgrey", "red"),
  dot.scale = 10
) +
  scale_x_discrete(position = "top", labels = gene_labels) +
  scale_y_discrete(labels = c(
    "4R_P301S_Ctrl" = "Ctrl",
    "4R_P301S_K18"  = "Seeded"
  )) +
  labs(
    x = NULL, y = NULL,
    title = paste0("MG DotPlot | AD-risk genes ∩ MG DEGs (n=", length(genes_use), ")")
  ) +
  theme(
    axis.text.y         = element_text(size = 20),
    axis.text.x.top     = element_text(size = 20, angle = 90, hjust = 0, vjust = 0.5),
    axis.text.x.bottom  = element_blank(),
    axis.ticks.x.bottom = element_blank()
  )

p_ad_top


## -----------------------------------------------------------------------------
##Fig 4d IFN genes intersection and DotPlot. 
#IFN gene list is downloaded from GSEA REACTOME interferon signaling:https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_INTERFERON_SIGNALING.html


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})


## -----------------------------
## 1) IFN gene list (Reactome)
## -----------------------------
IFN_genes <- as.data.frame(IFN_reactome)$GENE_SYMBOLS
IFN_genes <- unique(as.character(IFN_genes))
IFN_genes <- IFN_genes[!is.na(IFN_genes) & nzchar(IFN_genes)]

## -----------------------------
## 2) Intersect IFN genes with MG DEGs (same logic as your original)
## -----------------------------
deg <- as.data.frame(MG_4R301k_vs_4R301ctrl)
stopifnot("gene" %in% colnames(deg))



deg_ifn <- deg[deg$gene %in% IFN_genes, , drop = FALSE]
deg_ifn$Direction <- ifelse(deg_ifn$avg_log2FC >= 0, "Up", "Down")

deg_ifn <- deg_ifn[order(deg_ifn$avg_log2FC, decreasing = TRUE), , drop = FALSE]

message("IFN ∩ MG_DEG genes: ", nrow(deg_ifn))

## -----------------------------
## 3) Prepare plotting object: two conditions only
## -----------------------------
DefaultAssay(mg_obj) <- "RNA"
conditions_oi <- c("4R_P301S_Ctrl", "4R_P301S_K18")

mg2 <- subset(mg_obj, subset = Condition %in% conditions_oi)
mg2$Condition <- factor(mg2$Condition, levels = c("4R_P301S_K18", "4R_P301S_Ctrl"))  # Seeded bottom, Ctrl top
Idents(mg2) <- "Condition"

## -----------------------------
## 4) Genes to plot (keep deg_ifn order; drop missing genes)
## -----------------------------
genes_use <- unique(as.character(deg_ifn$gene))
genes_use <- genes_use[!is.na(genes_use) & nzchar(genes_use)]
genes_use <- genes_use[genes_use %in% rownames(mg2)]

if (length(genes_use) < 1) stop("No IFN∩DEG genes found in mg2 RNA rownames.")

gene_labels <- setNames(genes_use, genes_use)

## -----------------------------
## 5) DotPlot styled to match your figure (top x-axis, large labels)
## -----------------------------
p_ifn_top <- DotPlot(
  object    = mg2,
  features  = genes_use,
  group.by  = "Condition",
  cols      = c("lightgrey", "hotpink"),
  dot.scale = 10
) +
  scale_x_discrete(position = "top", labels = gene_labels) +
  scale_y_discrete(labels = c(
    "4R_P301S_Ctrl" = "Ctrl",
    "4R_P301S_K18"  = "Seeded"
  )) +
  labs(
    x = NULL, y = NULL,
    title = paste0("IFN genes (n=", length(genes_use), ")")
  ) +
  theme(
    axis.text.y         = element_text(size = 20),
    axis.text.x.top     = element_text(size = 20, angle = 90, hjust = 0, vjust = 0.5),
    axis.text.x.bottom  = element_blank(),
    axis.ticks.x.bottom = element_blank()
  )

p_ifn_top

#Fig. 3e
#Top 500 DEGs of MG were input into GSEA webset and Human Hallmark get sets were selected for computing overlaps
#top 20 pathways were selected for plotting

# -----------------------------
# 0) Input
# -----------------------------
df <- as.data.frame(MG_up500_hallmark_top20)

# -----------------------------
# 1) Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(readr)   # parse_number()
})

# -----------------------------
# 2) Clean + format table
#    (column names match your object exactly)
# -----------------------------
df2 <- df %>%
  rename(
    pathway_raw = `Gene Set Name`,
    K           = `# Genes in Gene Set (K)`,
    k           = `# Genes in Overlap (k)`,
    FDR         = `FDR q-value`
  ) %>%
  mutate(
    K = parse_number(as.character(K)),
    k = parse_number(as.character(k)),
    FDR = as.numeric(FDR),
    
    # display names
    pathway = pathway_raw %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_to_title(),
    
    # x-axis value in your figure
    score = -log10(pmax(FDR, 1e-300)),
    
    # text shown on bars (k/K)
    kk_label = paste0(k, "/", K)
  )

# -----------------------------
# 3) Select top pathways (Top 20)
# -----------------------------
n_top <- 20
plot_top <- df2 %>%
  arrange(FDR, desc(score)) %>%
  slice_head(n = n_top) %>%
  mutate(pathway = factor(pathway, levels = rev(pathway)))

# -----------------------------
# 4) Plot (matches your style)
# -----------------------------
p_hallmark <- ggplot(plot_top, aes(x = score, y = pathway)) +
  geom_col(fill = "#F79C7B") +
  geom_text(aes(label = kk_label), hjust = 1, nudge_x = -0.5, size = 3.5) +
  labs(
    x = expression(-log[10]~"(FDR)"),
    y = NULL,
    title = "Top 500 upregualted MG DEGs-Hallmark pathway enrichment"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_classic(base_size = 12)

p_hallmark

