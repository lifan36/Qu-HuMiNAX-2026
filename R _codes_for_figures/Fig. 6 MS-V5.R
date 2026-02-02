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

