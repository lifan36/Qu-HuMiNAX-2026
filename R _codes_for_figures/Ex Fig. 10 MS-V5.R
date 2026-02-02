suppressPackageStartupMessages({
  library(dplyr)
  library(ggVennDiagram)
  library(ggplot2)
  library(patchwork)
})

# ---- ensure data.frames ----
MG_PGRNKI_2_vs_301k    <- as.data.frame(MG_PGRNKI_2_vs_301k)
MG_PGRNKIiso_vs_301k   <- as.data.frame(MG_PGRNKIiso_vs_301k)

AST_PGRNKI_2_vs_301k   <- as.data.frame(AST_PGRNKI_2_vs_301k)
AST_PGRNKIiso_vs_301k  <- as.data.frame(AST_PGRNKIiso_vs_301k)

NEU_PGRNKI_2_vs_301k   <- as.data.frame(NEU_PGRNKI_2_vs_301k)
NEU_PGRNKIiso_vs_301k  <- as.data.frame(NEU_PGRNKIiso_vs_301k)

# ============================================================
# A) MICROGLIA: p1 | p2
# ============================================================
deg1 <- MG_PGRNKI_2_vs_301k   %>% filter(p_val_adj < 0.05)
deg3 <- MG_PGRNKIiso_vs_301k  %>% filter(p_val_adj < 0.05)

up1   <- deg1 %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
down1 <- deg1 %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()

up3   <- deg3 %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
down3 <- deg3 %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()

up_list <- list(
  PGRNKI_2_vs_301k   = up1,
  PGRNKIiso_vs_301k  = up3
)

down_list <- list(
  PGRNKI_2_vs_301k   = down1,
  PGRNKIiso_vs_301k  = down3
)

p1 <- ggVennDiagram(up_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "hotpink") +
  labs(title = "Overlap of UP DEGs") +
  theme(plot.title = element_text(size = 14, face = "bold"))

p2 <- ggVennDiagram(down_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "darkseagreen1") +
  labs(title = "Overlap of DOWN DEGs") +
  theme(plot.title = element_text(size = 14, face = "bold"))

p_mg <- p1 | p2
p_mg


# ============================================================
# B) ASTROCYTE (AST): p1 | p2
# ============================================================
deg1 <- AST_PGRNKI_2_vs_301k   %>% filter(p_val_adj < 0.05)
deg3 <- AST_PGRNKIiso_vs_301k  %>% filter(p_val_adj < 0.05)

up1   <- deg1 %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
down1 <- deg1 %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()

up3   <- deg3 %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
down3 <- deg3 %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()

up_list <- list(
  PGRNKI_2_vs_301k   = up1,
  PGRNKIiso_vs_301k  = up3
)

down_list <- list(
  PGRNKI_2_vs_301k   = down1,
  PGRNKIiso_vs_301k  = down3
)

p1 <- ggVennDiagram(up_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "orchid1") +
  labs(title = "Overlap of UP DEGs") +
  theme(plot.title = element_text(size = 14, face = "bold"))

p2 <- ggVennDiagram(down_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "turquoise") +
  labs(title = "Overlap of DOWN DEGs") +
  theme(plot.title = element_text(size = 14, face = "bold"))

p_ast <- p1 | p2
p_ast





# ============================================================
# C) NEURON (NEU): p1 | p2
# ============================================================
deg1 <- NEU_PGRNKI_2_vs_301k   %>% filter(p_val_adj < 0.05)
deg3 <- NEU_PGRNKIiso_vs_301k  %>% filter(p_val_adj < 0.05)

up1   <- deg1 %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
down1 <- deg1 %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()

up3   <- deg3 %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
down3 <- deg3 %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()

up_list <- list(
  PGRNKI_2_vs_301k   = up1,
  PGRNKIiso_vs_301k  = up3
)

down_list <- list(
  PGRNKI_2_vs_301k   = down1,
  PGRNKIiso_vs_301k  = down3
)

p1 <- ggVennDiagram(up_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#FF6347") +
  labs(title = "Overlap of UP DEGs") +
  theme(plot.title = element_text(size = 14, face = "bold"))

p2 <- ggVennDiagram(down_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Overlap of DOWN DEGs") +
  theme(plot.title = element_text(size = 14, face = "bold"))

p_neu <- p1 | p2
p_neu




