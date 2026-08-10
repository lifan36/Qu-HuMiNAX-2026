library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(grid)

############################################################
## Set desired_base levels for control_celltype
############################################################

ct_all <- unique(as.character(xeno_wt$control_celltype))
ct_all <- ct_all[!is.na(ct_all)]

EN_levels  <- sort(ct_all[grepl("^Neu_EN_", ct_all)])
IN_levels  <- sort(ct_all[grepl("^Neu_IN_", ct_all)])
NE_levels  <- intersect("Neu_Early", ct_all)
AST_levels <- intersect("AST", ct_all)
MG_levels  <- intersect("MG", ct_all)
NPC_levels <- intersect("NPC_RG", ct_all)

desired_base <- c(EN_levels, IN_levels, NE_levels, AST_levels, MG_levels, NPC_levels)
desired_base <- c(desired_base, setdiff(ct_all, desired_base))

xeno_wt$control_celltype <- factor(
  xeno_wt$control_celltype,
  levels = desired_base
)

############################################################
## Make 3R WT control object
############################################################

xeno_3R <- subset(
  xeno_wt,
  subset = Condition == "3R_WT_Ctrl"
)

xeno_3R$Sample_Name <- droplevels(factor(xeno_3R$Sample_Name))
xeno_3R$Condition <- droplevels(factor(xeno_3R$Condition))
xeno_3R$control_celltype <- droplevels(factor(xeno_3R$control_celltype))

Idents(xeno_3R) <- droplevels(xeno_3R$control_celltype)

DefaultAssay(xeno_3R) <- "RNA"

table(xeno_3R$Sample_Name)
table(xeno_3R$Condition)
unique(xeno_3R$control_celltype)

############################################################
## Build percentage labels for cell type UMAP S4A
############################################################

ct <- as.character(xeno_3R$control_celltype)
pct <- prop.table(table(ct)) * 100

base_present <- desired_base[desired_base %in% names(pct)]

lab_map <- setNames(
  sprintf("%s (%.1f%%)", base_present, as.numeric(pct[base_present])),
  base_present
)

vals <- unname(lab_map[ct])
names(vals) <- colnames(xeno_3R)

xeno_3R <- AddMetaData(
  xeno_3R,
  metadata = vals,
  col.name = "control_celltype_pct"
)

xeno_3R$control_celltype_pct <- factor(
  xeno_3R$control_celltype_pct,
  levels = unname(lab_map)
)

levels(xeno_3R$control_celltype_pct)

p_f_3R <- DimPlot(
  xeno_3R,
  group.by = "control_celltype_pct",
  label = TRUE,
  repel = TRUE,
  pt.size = 1.2
) +
  NoAxes()

print(p_f_3R)

############################################################
## Cell type percentage by sample S4B
## Keep NPC_RG and VLMC as separate categories
############################################################

xeno_3R$broad_celltype <- case_when(
  grepl("^Neu", as.character(xeno_3R$control_celltype)) ~ "Neu",
  as.character(xeno_3R$control_celltype) == "AST" ~ "AST",
  as.character(xeno_3R$control_celltype) == "MG" ~ "MG",
  as.character(xeno_3R$control_celltype) == "NPC_RG" ~ "NPC_RG",
  as.character(xeno_3R$control_celltype) == "VLMC" ~ "VLMC",
  TRUE ~ NA_character_
)

sample_vec <- as.character(xeno_3R$Sample_Name)
broad_vec  <- as.character(xeno_3R$broad_celltype)

keep <- !is.na(sample_vec) & !is.na(broad_vec)

sample_broad_counts <- table(
  sample_vec[keep],
  broad_vec[keep]
)

sample_broad_pct <- prop.table(sample_broad_counts, margin = 1) * 100

celltype_order <- c("Neu", "AST", "MG", "NPC_RG", "VLMC")

sample_broad_pct_df <- as.data.frame.matrix(sample_broad_pct)

missing_cols <- setdiff(celltype_order, colnames(sample_broad_pct_df))
sample_broad_pct_df[missing_cols] <- 0

sample_broad_pct_df <- sample_broad_pct_df[, celltype_order, drop = FALSE] %>%
  rownames_to_column("Sample_Name") %>%
  mutate(across(all_of(celltype_order), ~ round(.x, 2)))

print(sample_broad_pct_df)

write.csv(
  sample_broad_pct_df,
  file = "3R_WT_Ctrl_broad_celltype_percent_by_sample.csv",
  row.names = FALSE
)


############################################################
## MG marker analysis S4C
############################################################

mg <- subset(
  xeno_3R,
  subset = control_celltype == "MG"
)

DefaultAssay(mg) <- "RNA"

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

set.seed(42)

mg <- AddModuleScore(
  object = mg,
  features = list(paper_mg_maturity_genes_use),
  name = "Paper_MG_Maturity_Score"
)

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
      plot.title = element_text(hjust = 0.5, size = 20),
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
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

final_plots <- c(gene_plots, list(score_plot))

p_mg_maturation_3x3_3R <- wrap_plots(
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

p_mg_maturation_3x3_3R

############################################################
## AST maturation analysis S4D
############################################################

ast <- subset(
  xeno_3R,
  subset = control_celltype == "AST"
)

DefaultAssay(ast) <- "RNA"

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

ast_maturation_genes <- wang_mature_ast_genes

ast_maturation_genes_use <- intersect(
  ast_maturation_genes,
  rownames(ast)
)

setdiff(ast_maturation_genes, rownames(ast))

genes_present <- intersect(wang_mature_ast_genes, rownames(ast))
genes_missing <- setdiff(wang_mature_ast_genes, rownames(ast))

genes_missing
length(genes_present)

expr_mat <- FetchData(
  ast,
  vars = genes_present,
  layer = "data"
)

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

set.seed(42)

ast <- AddModuleScore(
  object = ast,
  features = list(ast_maturation_genes_use),
  name = "AST_Maturation_Score"
)

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

purple_high <- "#7B3F98"

small_legend <- guides(
  color = guide_colorbar(
    barheight = unit(2, "cm"),
    barwidth = unit(.5, "cm"),
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

final_plots <- c(gene_plots, list(score_plot))

p_ast_maturation_8panel_3R <- wrap_plots(
  final_plots,
  nrow = 2
) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 14)
    )
  ) &
  theme(
    plot.margin = margin(1, 1)
  )

p_ast_maturation_2by4panel_3R
