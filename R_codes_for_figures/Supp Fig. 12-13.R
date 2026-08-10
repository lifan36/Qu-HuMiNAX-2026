#Ex. Fig 11

## =============================================================================
## 0) Packages
## =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(VennDiagram)
  library(ggpubr)
  library(ggrepel)
  library(ggVennDiagram)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
})

## =============================================================================
## 1) Microglia panels a to d
## Inputs:
##   MG_PGRNKI_2_vs_301k
##   MG_PGRNKIiso_vs_301k
##   MG_4R301k_vs_4R301ctrl
## =============================================================================

## 1.0 Coerce inputs to data.frame (exactly as used)
MG_PGRNKI_2_vs_301k       <- as.data.frame(MG_PGRNKI_2_vs_301k)
MG_PGRNKIiso_vs_301k      <- as.data.frame(MG_PGRNKIiso_vs_301k)
MG_4R301k_vs_4R301ctrl    <- as.data.frame(MG_4R301k_vs_4R301ctrl)

## 1.1 Input checks

required_cols <- c("gene", "p_val_adj", "avg_log2FC")
for (nm in c("MG_PGRNKI_2_vs_301k","MG_PGRNKIiso_vs_301k","MG_4R301k_vs_4R301ctrl")) {
  missing_cols <- setdiff(required_cols, colnames(get(nm)))
  if (length(missing_cols) > 0) {
    stop(paste0(nm, " is missing columns: ", paste(missing_cols, collapse = ", ")))
  }
}

## 1.2 Define significant DEGs (FDR < 0.05)
deg1 <- MG_PGRNKI_2_vs_301k     %>% filter(p_val_adj < 0.05)
deg2 <- MG_4R301k_vs_4R301ctrl  %>% filter(p_val_adj < 0.05)
deg3 <- MG_PGRNKIiso_vs_301k    %>% filter(p_val_adj < 0.05)

gene_list <- list(
  PGRNKI2_vs_301k    = unique(deg1$gene),
  p301k_vs_ctrl      = unique(deg2$gene),
  PGRNKIiso_vs_301k  = unique(deg3$gene)
)

## 1.3 Panel a: 3 way Venn (all significant DEGs, direction ignored)
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

## 1.4 Panel b: same direction overlap between the two PGRN comparisons, then correlation
up1   <- deg1 %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
down1 <- deg1 %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()

up3   <- deg3 %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
down3 <- deg3 %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()

up_same   <- intersect(up1, up3)
down_same <- intersect(down1, down3)
same_dir_genes <- union(up_same, down_same)

MG_PGRNKIiso_vs_301k_sameDir <- MG_PGRNKIiso_vs_301k %>%
  filter(gene %in% same_dir_genes)

p301k_vs_ctrl <- MG_4R301k_vs_4R301ctrl

df2 <- merge(
  p301k_vs_ctrl[, c("gene", "avg_log2FC")],
  MG_PGRNKIiso_vs_301k_sameDir[, c("gene", "avg_log2FC")],
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

ggscatter(
  df2,
  x = "avg_log2FC.x",
  y = "avg_log2FC.y",
  color = grDevices::adjustcolor("hotpink", alpha.f = 0.3),
  add = "reg.line",
  add.params = list(
    color = "black",
    size  = 0.5,
    fill  = grDevices::adjustcolor("grey", alpha.f = 0.3)
  ),
  conf.int = TRUE,
  cor.coef = TRUE,
  cor.method = "pearson",
  xlab = "log2FC MG p301k vs control",
  ylab = "log2FC MG PGRNKIiso vs 301k (same direction genes)",
  repel = TRUE
) +
  # geom_text_repel(aes(label = label_to_display), size = 3, max.overlaps = 150) +
  geom_hline(yintercept =  0.1, linetype = "dashed", color = "grey10") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey10") +
  geom_vline(xintercept =  0.1, linetype = "dashed", color = "grey10") +
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey10")




## =============================================================================
## 2) Astrocyte panels e to h
## Inputs:
##   AST_PGRNKI_2_vs_301k
##   AST_PGRNKIiso_vs_301k
##   AST_4R301k_vs_4R301ctrl
## =============================================================================

## 2.0 Input checks and coercion
stopifnot(
  exists("AST_PGRNKI_2_vs_301k"),
  exists("AST_PGRNKIiso_vs_301k"),
  exists("AST_4R301k_vs_4R301ctrl")
)

AST_PGRNKI_2_vs_301k     <- as.data.frame(AST_PGRNKI_2_vs_301k)
AST_PGRNKIiso_vs_301k    <- as.data.frame(AST_PGRNKIiso_vs_301k)
AST_4R301k_vs_4R301ctrl  <- as.data.frame(AST_4R301k_vs_4R301ctrl)

for (nm in c("AST_PGRNKI_2_vs_301k","AST_PGRNKIiso_vs_301k","AST_4R301k_vs_4R301ctrl")) {
  missing_cols <- setdiff(required_cols, colnames(get(nm)))
  if (length(missing_cols) > 0) {
    stop(paste0(nm, " is missing columns: ", paste(missing_cols, collapse = ", ")))
  }
}

## 2.1 Define significant DEGs (FDR < 0.05)
deg1 <- AST_PGRNKI_2_vs_301k     %>% filter(p_val_adj < 0.05)
deg2 <- AST_4R301k_vs_4R301ctrl  %>% filter(p_val_adj < 0.05)
deg3 <- AST_PGRNKIiso_vs_301k    %>% filter(p_val_adj < 0.05)

gene_list <- list(
  PGRNKI2_vs_301k    = unique(deg1$gene),
  p301k_vs_ctrl      = unique(deg2$gene),
  PGRNKIiso_vs_301k  = unique(deg3$gene)
)

## 2.2 Panel e: 3 way Venn (all significant DEGs, direction ignored)
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

## 2.3 Panel f: same direction overlap between the two PGRN comparisons, then correlation
up1   <- deg1 %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
down1 <- deg1 %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()

up3   <- deg3 %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
down3 <- deg3 %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()

up_same   <- intersect(up1, up3)
down_same <- intersect(down1, down3)
same_dir_genes <- union(up_same, down_same)

AST_PGRNKIiso_vs_301k_sameDir <- AST_PGRNKIiso_vs_301k %>%
  filter(gene %in% same_dir_genes)

p301k_vs_ctrl <- AST_4R301k_vs_4R301ctrl

df2 <- merge(
  p301k_vs_ctrl[, c("gene", "avg_log2FC")],
  AST_PGRNKIiso_vs_301k_sameDir[, c("gene", "avg_log2FC")],
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

ggscatter(
  df2,
  x = "avg_log2FC.x",
  y = "avg_log2FC.y",
  color = grDevices::adjustcolor("darkorange", alpha.f = 0.3),
  add = "reg.line",
  add.params = list(
    color = "black",
    size  = 0.5,
    fill  = grDevices::adjustcolor("grey", alpha.f = 0.3)
  ),
  conf.int = TRUE,
  cor.coef = TRUE,
  cor.method = "pearson",
  xlab = "log2FC AST p301k vs control",
  ylab = "log2FC AST PGRNKIiso vs 301k (same direction genes)",
  repel = TRUE
) +
  # geom_text_repel(aes(label = label_to_display), size = 3, max.overlaps = 150) +
  geom_hline(yintercept =  0.1, linetype = "dashed", color = "grey10") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey10") +
  geom_vline(xintercept =  0.1, linetype = "dashed", color = "grey10") +
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey10")


#Microglia and astrocyte mirror pathways 
## =============================================================================
## Shared helpers (define once, reuse for MG + AST)
## =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

if (!exists("normalize_symbols")) {
  normalize_symbols <- function(x) {
    x <- as.character(x); x <- trimws(x)
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
    df$label <- sprintf("%s / %s",
                        df$Count,
                        as.numeric(sub("/.*", "", df$BgRatio)))
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

if (!exists("plot_mirror_ranked_by_K18_signed")) {
  plot_mirror_ranked_by_K18_signed <- function(ov, title,
                                               k18_sign = +1, other_sign = -1,
                                               k18_color = "#CD6090", other_color = "#4E79A7") {
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
}


## -----------------------------------------------------------------------------
## 1) Color schemes (matched to Extended Fig 10 pastel tones)

## -----------------------------------------------------------------------------
scheme_MG_up   <- list(k18 = "#f5c8f5", other = "#dce1be")  # right lavender, left light green
scheme_MG_down <- list(k18 = "#a5e1e1", other = "#f0be9b")  # left teal, right peach

scheme_AST_up   <- list(k18 = "#f5c8cd", other = "#bee1cd") # right pink, left light green
scheme_AST_down <- list(k18 = "#a5b9e1", other = "#f0a59b") # left periwinkle, right salmon


## =============================================================================
## MG
## Requires:
##   MG_4R301k_vs_4R301ctrl
##   MG_PGRNKIiso_vs_301k_sameDir
## =============================================================================
# stopifnot(exists("MG_4R301k_vs_4R301ctrl"), exists("MG_PGRNKIiso_vs_301k_sameDir"))
# stopifnot(all(req_cols %in% colnames(MG_4R301k_vs_4R301ctrl)))
# stopifnot(all(req_cols %in% colnames(MG_PGRNKIiso_vs_301k_sameDir)))

## 1.1 Build DEG sets
mg_k18_sets <- make_updown(MG_4R301k_vs_4R301ctrl,       fdr = 0.05)
mg_iso_sets <- make_updown(MG_PGRNKIiso_vs_301k_sameDir, fdr = 0.05)

## 1.2 GO BP enrichment (rank side equals K18)
mg_k18_up_bp   <- do_go_bp(mg_k18_sets$up)
mg_k18_down_bp <- do_go_bp(mg_k18_sets$down)
mg_iso_up_bp   <- do_go_bp(mg_iso_sets$up)
mg_iso_down_bp <- do_go_bp(mg_iso_sets$down)

MG_K18_UP_all <- go_to_df_all(mg_k18_up_bp)
MG_K18_DN_all <- go_to_df_all(mg_k18_down_bp)
MG_ISO_UP_all <- go_to_df_all(mg_iso_up_bp)
MG_ISO_DN_all <- go_to_df_all(mg_iso_down_bp)

## 1.3 Overlaps ranked by K18
ov_MG_UPK18_DNres <- get_rescued_overlap_rank_by_K18(MG_K18_UP_all, MG_ISO_DN_all, top_n = 15)
ov_MG_DNK18_UPres <- get_rescued_overlap_rank_by_K18(MG_K18_DN_all, MG_ISO_UP_all, top_n = 15)

## 1.4 Mirror plots (colors match MG panels c and d)
p_MG_rescue1 <- plot_mirror_ranked_by_K18_signed(
  ov_MG_UPK18_DNres,
  "MG mirror pathways: Up in p301k vs ctrl AND Down in PGRN OE vs 301k",
  k18_sign    = +1,
  other_sign  = -1,
  k18_color   = scheme_MG_up$k18,
  other_color = scheme_MG_up$other
)

p_MG_rescue2 <- plot_mirror_ranked_by_K18_signed(
  ov_MG_DNK18_UPres,
  "MG mirror pathways: Down in p301k vs ctrl AND Up in PGRN OE vs 301k",
  k18_sign    = -1,
  other_sign  = +1,
  k18_color   = scheme_MG_down$k18,
  other_color = scheme_MG_down$other
)

p_MG_rescue1
p_MG_rescue2


## =============================================================================
## AST
## Requires:
##   AST_4R301k_vs_4R301ctrl
##   AST_PGRNKIiso_vs_301k_sameDir
## =============================================================================

## 2.1 Build DEG sets
ast_k18_sets <- make_updown(AST_4R301k_vs_4R301ctrl,       fdr = 0.05)
ast_iso_sets <- make_updown(AST_PGRNKIiso_vs_301k_sameDir, fdr = 0.05)

## 2.2 GO BP enrichment (rank side equals K18)
ast_k18_up_bp   <- do_go_bp(ast_k18_sets$up)
ast_k18_down_bp <- do_go_bp(ast_k18_sets$down)
ast_iso_up_bp   <- do_go_bp(ast_iso_sets$up)
ast_iso_down_bp <- do_go_bp(ast_iso_sets$down)

AST_K18_UP_all <- go_to_df_all(ast_k18_up_bp)
AST_K18_DN_all <- go_to_df_all(ast_k18_down_bp)
AST_ISO_UP_all <- go_to_df_all(ast_iso_up_bp)
AST_ISO_DN_all <- go_to_df_all(ast_iso_down_bp)

## 2.3 Overlaps ranked by K18
ov_AST_UPK18_DNres_ISO <- get_rescued_overlap_rank_by_K18(AST_K18_UP_all, AST_ISO_DN_all, top_n = 15)
ov_AST_DNK18_UPres_ISO <- get_rescued_overlap_rank_by_K18(AST_K18_DN_all, AST_ISO_UP_all, top_n = 15)

## 2.4 Mirror plots (colors match AST panels g and h)
p_AST_ISO_rescue1 <- plot_mirror_ranked_by_K18_signed(
  ov_AST_UPK18_DNres_ISO,
  "AST mirror pathways (ISO sameDir): Up in p301k vs ctrl AND Down in PGRN OE vs 301k",
  k18_sign    = +1,
  other_sign  = -1,
  k18_color   = scheme_AST_up$k18,
  other_color = scheme_AST_up$other
)

p_AST_ISO_rescue2 <- plot_mirror_ranked_by_K18_signed(
  ov_AST_DNK18_UPres_ISO,
  "AST mirror pathways (ISO sameDir): Down in p301k vs ctrl AND Up in PGRN OE vs 301k",
  k18_sign    = -1,
  other_sign  = +1,
  k18_color   = scheme_AST_down$k18,
  other_color = scheme_AST_down$other
)

p_AST_ISO_rescue1
p_AST_ISO_rescue2
