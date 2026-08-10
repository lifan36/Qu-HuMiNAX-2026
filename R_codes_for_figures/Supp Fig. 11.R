library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)

############################################################
## Use KI final label object
############################################################

obj <- xeno_wKI_final_label

table(xeno_wKI_final_label$Condition)
ncol(xeno_wKI_final_label)

DefaultAssay(obj) <- "RNA"

############################################################
## Use the same UMAP reduction as original KI plot
############################################################

reduction_use <- "umap_res0.3"

Reductions(obj)

if (!reduction_use %in% Reductions(obj)) {
  stop(paste0("Reduction ", reduction_use, " not found. Check Reductions(obj)."))
}

############################################################
## Set cell type order
############################################################

ct_all <- unique(as.character(obj$control_celltype))
ct_all <- ct_all[!is.na(ct_all)]

EN_levels  <- sort(ct_all[grepl("^Neu_EN_", ct_all)])
IN_levels  <- sort(ct_all[grepl("^Neu_IN_", ct_all)])
NE_levels  <- intersect("Neu_Early", ct_all)
AST_levels <- intersect("AST", ct_all)
MG_levels  <- intersect("MG", ct_all)
NPC_levels <- intersect("NPC_RG", ct_all)

desired_base <- c(
  EN_levels,
  IN_levels,
  NE_levels,
  AST_levels,
  MG_levels,
  NPC_levels
)

desired_base <- c(desired_base, setdiff(ct_all, desired_base))

obj$control_celltype <- factor(
  obj$control_celltype,
  levels = desired_base
)

############################################################
## Set condition order
## ISOKI comes before KI
############################################################

condition_order <- c(
  "4R_P301S_Ctrl",
  "4R_P301S_K18",
  "ISOKI_4R_P301S_K18",
  "KI_4R_P301S_K18"
)

obj$Condition <- factor(
  obj$Condition,
  levels = condition_order
)

############################################################
## Helper function
############################################################

make_pct_umap <- function(seu, plot_title = NULL, reduction_use = "umap_res0.3") {
  
  if ("Sample_Name" %in% colnames(seu@meta.data)) {
    seu$Sample_Name <- droplevels(factor(seu$Sample_Name))
  }
  
  seu$Condition <- droplevels(factor(seu$Condition))
  seu$control_celltype <- droplevels(factor(seu$control_celltype))
  
  Idents(seu) <- seu$control_celltype
  
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
  ## Label positions from same UMAP reduction
  ############################################################
  
  umap_df <- as.data.frame(Embeddings(seu, reduction = reduction_use))
  umap_df <- umap_df[, 1:2]
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  
  umap_df$control_celltype_pct <- seu$control_celltype_pct
  
  label_df <- umap_df %>%
    group_by(control_celltype_pct) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      .groups = "drop"
    )
  
  ############################################################
  ## Plot
  ############################################################
  
  p <- DimPlot(
    seu,
    reduction = reduction_use,
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
      inherit.aes = FALSE,
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
## All conditions UMAP
############################################################

xeno_KI_all <- subset(
  obj,
  subset = Condition %in% condition_order
)

p_KI_all <- make_pct_umap(
  xeno_KI_all,
  plot_title = "All conditions",
  reduction_use = reduction_use
)

############################################################
## ISOKI and KI only
############################################################

xeno_ISOKI_4R_P301S_K18 <- subset(
  obj,
  subset = Condition == "ISOKI_4R_P301S_K18"
)

xeno_KI_4R_P301S_K18 <- subset(
  obj,
  subset = Condition == "KI_4R_P301S_K18"
)

p_ISOKI_4R_P301S_K18 <- make_pct_umap(
  xeno_ISOKI_4R_P301S_K18,
  plot_title = "ISOKI 4R P301S + Tau seeds",
  reduction_use = reduction_use
)

p_KI_4R_P301S_K18 <- make_pct_umap(
  xeno_KI_4R_P301S_K18,
  plot_title = "KI 4R P301S + Tau seeds",
  reduction_use = reduction_use
)

############################################################
## Combine into one row
## Order: all conditions, ISOKI, KI
############################################################

p_KI_combined <- 
  p_KI_all |
  p_ISOKI_4R_P301S_K18 |
  p_KI_4R_P301S_K18

p_KI_combined <- p_KI_combined +
  plot_layout(nrow = 1, widths = c(1, 1, 1)) &
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "plain"
    ),
    legend.title = element_blank()
  )

print(p_KI_combined)

############################################################
## Save
############################################################

ggsave(
  filename = "KI_ISOKI_celltype_percentage_UMAP_panel_umap_res0.3.pdf",
  plot = p_KI_combined,
  width = 15,
  height = 5
)
