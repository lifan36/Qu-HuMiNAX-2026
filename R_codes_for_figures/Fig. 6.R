#common rescued DEGs by PGRN OE by both cell lines 

NEU_PGRNKI_2_vs_301k <- as.data.frame(NEU_PGRNKI_2_vs_301k)
NEU_PGRNKIiso_vs_301k <-as.data.frame(NEU_PGRNKIiso_vs_301k)
NEU_p301k_vs_ct <- as.data.frame(NEU_p301k_vs_ct)

head(NEU_p301k_vs_ct)


library(dplyr)

# -----------------------------
# Settings
# -----------------------------
padj_cutoff <- 0.05
logfc_cutoff <- 0.1

# -----------------------------
# Make sure each DEG table has a clean gene column
# -----------------------------
NEU_PGRNKI_2_vs_301k <- NEU_PGRNKI_2_vs_301k %>%
  as.data.frame() %>%
  mutate(gene = ifelse(is.na(gene), rownames(.), gene))

NEU_PGRNKIiso_vs_301k <- NEU_PGRNKIiso_vs_301k %>%
  as.data.frame() %>%
  mutate(gene = ifelse(is.na(gene), rownames(.), gene))

NEU_p301k_vs_ct <- NEU_p301k_vs_ct %>%
  as.data.frame() %>%
  mutate(gene = ifelse(is.na(gene), rownames(.), gene))

# -----------------------------
# Optional: define significant DEGs
# If your input dataframes are already filtered DEGs, this step is still safe
# -----------------------------
ki2_deg <- NEU_PGRNKI_2_vs_301k %>%
  filter(p_val_adj < padj_cutoff, abs(avg_log2FC) > logfc_cutoff) %>%
  mutate(direction_KI2 = sign(avg_log2FC))

iso_deg <- NEU_PGRNKIiso_vs_301k %>%
  filter(p_val_adj < padj_cutoff, abs(avg_log2FC) > logfc_cutoff) %>%
  mutate(direction_iso = sign(avg_log2FC))

p301k_deg <- NEU_p301k_vs_ct %>%
  filter(p_val_adj < padj_cutoff, abs(avg_log2FC) > logfc_cutoff) %>%
  mutate(direction_p301k = sign(avg_log2FC))

# -----------------------------
# Step 1:
# Find overlapping DEGs between PGRN KI iso and PGRN KI 2
# that share the same avg_log2FC direction
# -----------------------------
common_deg_genes <- iso_deg %>%
  select(
    gene,
    avg_log2FC_iso = avg_log2FC,
    p_val_adj_iso = p_val_adj,
    direction_iso
  ) %>%
  inner_join(
    ki2_deg %>%
      select(
        gene,
        avg_log2FC_KI2 = avg_log2FC,
        p_val_adj_KI2 = p_val_adj,
        direction_KI2
      ),
    by = "gene"
  ) %>%
  filter(direction_iso == direction_KI2)

# -----------------------------
# Step 2:
# Subset this common gene list within NEU_PGRNKIiso_vs_301k
# Call this common_DEGs
# -----------------------------
common_DEGs <- NEU_PGRNKIiso_vs_301k %>%
  filter(gene %in% common_deg_genes$gene) %>%
  left_join(
    common_deg_genes %>%
      select(gene, avg_log2FC_KI2, p_val_adj_KI2, direction_KI2),
    by = "gene"
  ) %>%
  rename(
    avg_log2FC_iso = avg_log2FC,
    p_val_adj_iso = p_val_adj
  ) %>%
  mutate(direction_iso = sign(avg_log2FC_iso))

# -----------------------------
# Step 3:
# Within common_DEGs, find rescued DEGs:
# genes that move in the opposite direction in NEU_p301k_vs_ct
# -----------------------------
rescued_DEGs <- common_DEGs %>%
  inner_join(
    p301k_deg %>%
      select(
        gene,
        avg_log2FC_p301k_vs_ct = avg_log2FC,
        p_val_adj_p301k_vs_ct = p_val_adj,
        direction_p301k
      ),
    by = "gene"
  ) %>%
  filter(direction_iso == -direction_p301k) %>%
  mutate(
    rescue_pattern = case_when(
      avg_log2FC_p301k_vs_ct > 0 & avg_log2FC_iso < 0 ~ "Disease_up_PGRN_down",
      avg_log2FC_p301k_vs_ct < 0 & avg_log2FC_iso > 0 ~ "Disease_down_PGRN_up",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(p_val_adj_iso)

# -----------------------------
# Print summary
# -----------------------------
cat("Common same-direction DEGs between PGRN KI iso and PGRN KI 2:", 
    nrow(common_DEGs), "\n")

cat("Rescued DEGs with opposite direction to p301k vs ct:", 
    nrow(rescued_DEGs), "\n")

print(head(rescued_DEGs, 20))

# -----------------------------
# Export results
# -----------------------------

write.csv(
  rescued_DEGs,
  "NEU_rescued_DEGs_common_PGRN_opposite_to_p301k_vs_ct.csv",
  row.names = FALSE
)


rescued_DEGs_NEU <-as.data.frame(NEU_rescued_DEGs_common_PGRN_opposite_to_p301k_vs_ct) 

head(rescued_DEGs_NEU)







head(rescued_DEGs_NEU)








head(MG_PGRNKIiso_vs_301k)

MG_PGRNKIiso_vs_301k<- as.data.frame(MG_PGRNKIiso_vs_301k)
MG_p301k_vs_ct <- as.data.frame(MG_p301k_vs_ct)

MG_PGRNKI_2_vs_301k <-as.data.frame(MG_PGRNKI_2_vs_301k)

MG_PGRNKIiso_vs_301k<- as.data.frame(MG_PGRNKIiso_vs_301k)


# ============================================================
# Microglia rescued DEG list
# Common PGRN KI iso and PGRN KI 2 DEGs
# Opposite direction to P301S K18 disease effect
# ============================================================

library(dplyr)
library(stringr)
library(tidyr)
library(readr)

# -----------------------------
# Convert DEG tables to dataframes
# -----------------------------
MG_p301k_vs_ct <- as.data.frame(MG_p301k_vs_ct)
MG_PGRNKI_2_vs_301k <- as.data.frame(MG_PGRNKI_2_vs_301k)
MG_PGRNKIiso_vs_301k <- as.data.frame(MG_PGRNKIiso_vs_301k)

# -----------------------------
# Settings
# -----------------------------
padj_cutoff <- 0.05
logfc_cutoff <- 0.1

# -----------------------------
# Make sure each DEG table has a clean gene column
# -----------------------------
MG_PGRNKI_2_vs_301k <- MG_PGRNKI_2_vs_301k %>%
  as.data.frame() %>%
  mutate(gene = ifelse(is.na(gene), rownames(.), gene))

MG_PGRNKIiso_vs_301k <- MG_PGRNKIiso_vs_301k %>%
  as.data.frame() %>%
  mutate(gene = ifelse(is.na(gene), rownames(.), gene))

MG_p301k_vs_ct <- MG_p301k_vs_ct %>%
  as.data.frame() %>%
  mutate(gene = ifelse(is.na(gene), rownames(.), gene))

# -----------------------------
# Define significant DEGs
# -----------------------------
ki2_deg <- MG_PGRNKI_2_vs_301k %>%
  filter(p_val_adj < padj_cutoff, abs(avg_log2FC) > logfc_cutoff) %>%
  mutate(direction_KI2 = sign(avg_log2FC))

iso_deg <- MG_PGRNKIiso_vs_301k %>%
  filter(p_val_adj < padj_cutoff, abs(avg_log2FC) > logfc_cutoff) %>%
  mutate(direction_iso = sign(avg_log2FC))

p301k_deg <- MG_p301k_vs_ct %>%
  filter(p_val_adj < padj_cutoff, abs(avg_log2FC) > logfc_cutoff) %>%
  mutate(direction_p301k = sign(avg_log2FC))

# -----------------------------
# Step 1:
# Find overlapping DEGs between PGRN KI iso and PGRN KI 2
# that share the same avg_log2FC direction
# -----------------------------
MG_common_deg_genes <- iso_deg %>%
  select(
    gene,
    avg_log2FC_iso = avg_log2FC,
    p_val_adj_iso = p_val_adj,
    direction_iso
  ) %>%
  inner_join(
    ki2_deg %>%
      select(
        gene,
        avg_log2FC_KI2 = avg_log2FC,
        p_val_adj_KI2 = p_val_adj,
        direction_KI2
      ),
    by = "gene"
  ) %>%
  filter(direction_iso == direction_KI2)

# -----------------------------
# Step 2:
# Subset this common gene list within MG_PGRNKIiso_vs_301k
# Call this MG_common_DEGs
# -----------------------------
MG_common_DEGs <- MG_PGRNKIiso_vs_301k %>%
  filter(gene %in% MG_common_deg_genes$gene) %>%
  left_join(
    MG_common_deg_genes %>%
      select(gene, avg_log2FC_KI2, p_val_adj_KI2, direction_KI2),
    by = "gene"
  ) %>%
  rename(
    avg_log2FC_iso = avg_log2FC,
    p_val_adj_iso = p_val_adj
  ) %>%
  mutate(direction_iso = sign(avg_log2FC_iso))

# -----------------------------
# Step 3:
# Within MG_common_DEGs, find rescued DEGs:
# genes that move in the opposite direction in MG_p301k_vs_ct
# -----------------------------
MG_rescued_DEGs <- MG_common_DEGs %>%
  inner_join(
    p301k_deg %>%
      select(
        gene,
        avg_log2FC_p301k_vs_ct = avg_log2FC,
        p_val_adj_p301k_vs_ct = p_val_adj,
        direction_p301k
      ),
    by = "gene"
  ) %>%
  filter(direction_iso == -direction_p301k) %>%
  mutate(
    rescue_pattern = case_when(
      avg_log2FC_p301k_vs_ct > 0 & avg_log2FC_iso < 0 ~ "Disease_up_PGRN_down",
      avg_log2FC_p301k_vs_ct < 0 & avg_log2FC_iso > 0 ~ "Disease_down_PGRN_up",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(p_val_adj_iso)

# -----------------------------
# Print summary
# -----------------------------
cat("MG common same-direction DEGs between PGRN KI iso and PGRN KI 2:",
    nrow(MG_common_DEGs), "\n")

cat("MG rescued DEGs with opposite direction to p301k vs ct:",
    nrow(MG_rescued_DEGs), "\n")

print(head(MG_rescued_DEGs, 20))

# -----------------------------
# Export MG rescued DEG result
# -----------------------------
write.csv(
  MG_rescued_DEGs,
  "MG_rescued_DEGs_common_PGRN_opposite_to_p301k_vs_ct.csv",
  row.names = FALSE
)

# Optional object name for downstream use
rescued_DEGs_MG <- as.data.frame(MG_rescued_DEGs_common_PGRN_opposite_to_p301k_vs_ct)

head(rescued_DEGs_MG)



#HNRNPK-AS1 TF analysis
promoter_seq <- "
CCGTTAGGTCCGATGGCAAGTATGGCACCAACCTTCGTCGAGGGCTCTGTTGACGGCGGCGCACAGGGACCAGTACCCACTGTATGTTTTCCCTCTGTACCTCACCACACACTTGTGCTCGTCCTGCTCCGACCAAAAGGAGAAGTTGAGCGTGTCTGTCACGAACACCCAGTTGACCGCGGCCTCGTCGGCCGCCCTGGGGTTCAGCTCATGAAGGGCTTTCCACCCCTCCACGCGCAGCTCTGGCCCCGCCGCCTTGGCCAGCAGCAGCTCTGCCACCCTCCGTACGCCTCCGCTGTCAATAAACACATCCCGACTGTTTTCTGCAATGAATTTAGAGGATTCCCTGGGATTTAGGAGCCCGTCCATTCTGCCCTTTGGGACCCGGCGGGCCAAACTCTTCCCGCCCTAGGCGGGTCAAGGTCCACTCTGCGCCTTGCGCTACTACGAACCGCGCGATCACAGGGCCTAAACAAGGTCCCGCCCTGCCAGGGAGTGGCACCGAGAACTTAGTGCAGCCAAGGTAAAGGTTCTCTCCCGCGGCGTTGCCCTGATTAGGAAACCAGAGAGACCCCGGGGCAGCTCGCTGGGCAGAGGCCGCAAGCACTTCCGGGAAAGCGGTCCCGTTCACGCCCTCTCGGGCGCGCTCTCAGGCTCCGCTTCCGGTTCCAAAGGTCTGGAAACTGCTAAGGGTCCGGCGTTCACCGGCCGTGAGTTGCCCACCCTGGGAGGTGCATGGCTTTGCCATGTTTTGCCTTTTCTCACCCCGCTTCGGACAGGTGGGACTCGAAGACGGCGACTCGCCGCCTTGCCCTGGGTCGCGATGCCCCCACCCCGCTCCCACCCCTCTCCTTTCCCGCAGCCTGAGCAAGACGCCTTCGGCAACGCCCCGCATCCGGGCTCCTTGGCCTCGAAGAGCATTATGGCCGTAGATCTGGGTGCTGAGGACTGAGCCACCCCCAGACTGCGACATGGGCGGCGGTGCCTCCTTCCCCAAGCCCCAGGGAGTGTTTTTTTGTTTGTTTTGTTTTTTGTTTTTGTTTAAGTGAATTCCTGGCCCATTTGGCCGAACGAACGT
"

promoter_seq_clean <- promoter_seq |>
  toupper() |>
  gsub("[^ATCGN]", "", x = _)

writeLines(
  c(">HNRNPK_AS1_promoter", promoter_seq_clean),
  "HNRNPK_AS1_promoter.fa"
)


#the promoter sequence was saved as a FASTA file named HNRNPK_AS1_promoter.fa. Curated vertebrate transcription factor motifs were downloaded from JASPAR CORE in MEME format, specifically the nonredundant vertebrate position frequency matrix file. This motif file and the promoter FASTA file were then uploaded to the FIMO web server.FIMO was run using a p value threshold of 1e-4 and scanning both DNA strands. The output file, fimo.tsv, listed predicted motif occurrences in the HNRNPK-AS1 promoter, including the matched TF motif, motif position, strand, score, and statistical significance. These predicted motifs were then used to nominate candidate TFs that may regulate HNRNPK-AS1 expression.

fimo <- read.delim(
  "G:/Other computers/My Laptop/sequencing analysis/MS/R code/fimo.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


fimo <- read_tsv(fimo_file, comment = "#")

head(fimo)

fimo <- as.data.frame(fimo)



#predicted TF overlay with resuced DEG list

rescued_DEGs<-NEU_rescued_DEGs_common_PGRN_opposite_to_p301k_vs_ct
rescued_DEGs<- as.data.frame(rescued_DEGs)


library(dplyr)
library(stringr)

fimo_tfs_split <- fimo %>%
  as.data.frame() %>%
  filter(!is.na(motif_alt_id)) %>%
  mutate(motif_alt_id_clean = str_trim(motif_alt_id)) %>%
  tidyr::separate_rows(motif_alt_id_clean, sep = "::") %>%
  mutate(motif_alt_id_upper = toupper(motif_alt_id_clean)) %>%
  distinct(motif_alt_id_clean, motif_alt_id_upper)

rescued_DEGs_tf_overlap_split <- rescued_DEGs %>%
  mutate(gene_upper = toupper(gene)) %>%
  inner_join(
    fimo_tfs_split,
    by = c("gene_upper" = "motif_alt_id_upper")
  ) %>%
  arrange(p_val_adj_iso)

print(rescued_DEGs_tf_overlap_split)

write.csv(
  rescued_DEGs_tf_overlap_split,
  "NEU_rescued_DEGs_overlap_FIMO_motif_alt_id_split_complex_motifs.csv",
  row.names = FALSE
)


#DEGs input

MG_rescued_DEGs_common_PGRN_opposite_to_p301k_vs_ct <- as.data.frame(MG_rescued_DEGs_common_PGRN_opposite_to_p301k_vs_ct)
NEU_rescued_DEGs_common_PGRN_opposite_to_p301k_vs_ct<-as.data.frame(NEU_rescued_DEGs_common_PGRN_opposite_to_p301k_vs_ct)
NEU_rescued_DEGs_overlap_FIMO_motif_alt_id_split_complex_motifs <- as.data.frame(NEU_rescued_DEGs_overlap_FIMO_motif_alt_id_split_complex_motifs)


unique(xeno_wKI_final_label$control_celltype)
obj <- xeno_wKI_final_label

conditions_use <- c(
  "4R_P301S_Ctrl",
  "4R_P301S_K18",  
  "ISOKI_4R_P301S_K18",
  "KI_4R_P301S_K18"
)


#somethings urls doesn't load, can be downloaded from website and manually loadinto R
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))


# ============================================================
# Reliable MG ligand -> neuronal TF pair table
#
# Goal:
#   Identify PGRN-rescued MG ligands predicted by NicheNet to regulate
#   neuronal HNRNPK-AS1-associated TFs.
#
# Ligand filters:
#   1. PGRN rescued in MG by both ISO and KI_2
#   2. Known NicheNet ligand
#   3. Expressed in MG
#   4. Linked to at least one neuronal TF by nonzero NicheNet regulatory potential
#
# TF filters:
#   1. From NEU_rescued_DEGs_overlap_FIMO_motif_alt_id_split_complex_motifs$gene
#   2. Rescued neuronal DEG
#   3. Expressed in neurons
#   4. Present in NicheNet target matrix
#   5. Linked to at least one eligible MG ligand
#
# Output:
#   CSV of all reliable ligand -> neuronal TF pairs
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(Matrix)
  library(tibble)
})

# ============================================================
# 0. Settings
# ============================================================

obj <- xeno_wKI_final_label

celltype_col  <- "control_celltype"
condition_col <- "Condition"

sender_celltype <- "MG"
receiver_pattern <- "^Neu_"

assay_use <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
layer_use <- "data"

# Biological expression filter
pct_expr_cutoff <- 0.05

outdir <- "NicheNet_reliable_MG_ligand_to_neuron_TF_pairs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (!exists("ligand_target_matrix")) {
  stop("ligand_target_matrix is not loaded.")
}

# ============================================================
# 1. Helper functions
# ============================================================

pick_gene_col <- function(df) {
  preferred <- c("gene", "genes", "gene_symbol", "gene...1", "gene...7", "gene...16")
  hit <- preferred[preferred %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  
  gene_like <- grep("^gene($|\\.\\.\\.)", colnames(df), value = TRUE)
  if (length(gene_like) > 0) return(gene_like[1])
  
  stop("No gene-like column found.")
}

clean_deg_table <- function(df, gene_col = NULL) {
  if (is.null(gene_col)) {
    gene_col <- pick_gene_col(df)
  }
  
  if (!"avg_log2FC_iso" %in% colnames(df)) df$avg_log2FC_iso <- NA_real_
  if (!"avg_log2FC_KI2" %in% colnames(df)) df$avg_log2FC_KI2 <- NA_real_
  if (!"rescue_pattern" %in% colnames(df)) df$rescue_pattern <- NA_character_
  if (!"rescue_lines" %in% colnames(df)) df$rescue_lines <- NA_character_
  if (!"p_val_adj_p301k_vs_ct" %in% colnames(df)) df$p_val_adj_p301k_vs_ct <- NA_real_
  if (!"p_val_adj_iso" %in% colnames(df)) df$p_val_adj_iso <- NA_real_
  if (!"p_val_adj_KI2" %in% colnames(df)) df$p_val_adj_KI2 <- NA_real_
  
  df %>%
    mutate(
      gene_symbol = toupper(str_trim(as.character(.data[[gene_col]]))),
      avg_log2FC_p301k_vs_ct = as.numeric(avg_log2FC_p301k_vs_ct),
      avg_log2FC_iso = as.numeric(avg_log2FC_iso),
      avg_log2FC_KI2 = as.numeric(avg_log2FC_KI2),
      p_val_adj_p301k_vs_ct = as.numeric(p_val_adj_p301k_vs_ct),
      p_val_adj_iso = as.numeric(p_val_adj_iso),
      p_val_adj_KI2 = as.numeric(p_val_adj_KI2),
      rescue_pattern = as.character(rescue_pattern),
      rescue_lines = as.character(rescue_lines)
    ) %>%
    filter(
      !is.na(gene_symbol),
      gene_symbol != "",
      !is.na(avg_log2FC_p301k_vs_ct)
    ) %>%
    arrange(desc(abs(avg_log2FC_p301k_vs_ct))) %>%
    distinct(gene_symbol, .keep_all = TRUE)
}

get_data_layer <- function(seu, assay_use, layer_use = "data") {
  tryCatch(
    GetAssayData(seu, assay = assay_use, layer = layer_use),
    error = function(e) GetAssayData(seu, assay = assay_use, slot = layer_use)
  )
}

get_expressed_gene_stats <- function(seu, cells_use, assay, layer = "data", genes = NULL) {
  mat <- get_data_layer(seu, assay_use = assay, layer_use = layer)
  
  original_rownames <- rownames(mat)
  rownames(mat) <- toupper(original_rownames)
  
  cells_use <- intersect(cells_use, colnames(mat))
  if (length(cells_use) == 0) {
    stop("No matching cells found for expression calculation.")
  }
  
  if (is.null(genes)) {
    genes_use <- rownames(mat)
  } else {
    genes_use <- intersect(toupper(genes), rownames(mat))
  }
  
  pct_expr <- Matrix::rowMeans(mat[genes_use, cells_use, drop = FALSE] > 0)
  avg_expr <- Matrix::rowMeans(mat[genes_use, cells_use, drop = FALSE])
  
  tibble(
    gene_symbol = names(pct_expr),
    pct_expr = as.numeric(pct_expr),
    avg_expr = as.numeric(avg_expr)
  )
}

# ============================================================
# 2. Define sender MG cells and receiver neuron cells
# ============================================================

meta <- obj@meta.data

sender_cells <- rownames(meta)[
  as.character(meta[[celltype_col]]) == sender_celltype
]

receiver_cells <- rownames(meta)[
  grepl(receiver_pattern, as.character(meta[[celltype_col]]))
]

cat("MG sender cells:", length(sender_cells), "\n")
cat("Neuron receiver cells:", length(receiver_cells), "\n")

if (length(sender_cells) == 0) stop("No MG sender cells found.")
if (length(receiver_cells) == 0) stop("No Neu_ receiver cells found.")

# ============================================================
# 3. Clean input DEG tables
# ============================================================

mg_tbl <- clean_deg_table(
  MG_rescued_DEGs_common_PGRN_opposite_to_p301k_vs_ct
)

tf_tbl <- clean_deg_table(
  NEU_rescued_DEGs_overlap_FIMO_motif_alt_id_split_complex_motifs,
  gene_col = "gene"
)

# ============================================================
# 4. Standardize NicheNet ligand-target matrix orientation
# Expected final orientation:
#   rows = target genes
#   columns = ligands
# ============================================================

lt <- ligand_target_matrix
rownames(lt) <- toupper(rownames(lt))
colnames(lt) <- toupper(colnames(lt))

candidate_mg_genes <- unique(mg_tbl$gene_symbol)
candidate_tf_genes <- unique(tf_tbl$gene_symbol)

score_ligands_in_cols <- length(intersect(candidate_mg_genes, colnames(lt))) *
  length(intersect(candidate_tf_genes, rownames(lt)))

score_ligands_in_rows <- length(intersect(candidate_mg_genes, rownames(lt))) *
  length(intersect(candidate_tf_genes, colnames(lt)))

if (score_ligands_in_cols >= score_ligands_in_rows) {
  ligand_target_matrix_use <- lt
  cat("NicheNet matrix orientation: targets rows, ligands columns\n")
} else {
  ligand_target_matrix_use <- t(lt)
  cat("NicheNet matrix orientation: transposed to targets rows, ligands columns\n")
}

all_targets <- rownames(ligand_target_matrix_use)
all_ligands <- colnames(ligand_target_matrix_use)

# ============================================================
# 5. Expression stats in MG and neurons
# ============================================================

mg_expr_stats <- get_expressed_gene_stats(
  seu = obj,
  cells_use = sender_cells,
  assay = assay_use,
  layer = layer_use,
  genes = candidate_mg_genes
) %>%
  rename(
    ligand_pct_expr_MG = pct_expr,
    ligand_avg_expr_MG = avg_expr
  )

neuron_expr_stats <- get_expressed_gene_stats(
  seu = obj,
  cells_use = receiver_cells,
  assay = assay_use,
  layer = layer_use,
  genes = candidate_tf_genes
) %>%
  rename(
    tf_pct_expr_neuron = pct_expr,
    tf_avg_expr_neuron = avg_expr
  )

# ============================================================
# 6. Define reliable MG ligand universe
# ============================================================

eligible_mg_ligands <- mg_tbl %>%
  mutate(
    ligand = gene_symbol,
    ligand_rescue_mean_log2FC = rowMeans(
      cbind(avg_log2FC_iso, avg_log2FC_KI2),
      na.rm = TRUE
    ),
    ligand_rescue_direction = case_when(
      avg_log2FC_p301k_vs_ct > 0 & avg_log2FC_iso < 0 & avg_log2FC_KI2 < 0 ~ "Disease_up_PGRN_down",
      avg_log2FC_p301k_vs_ct < 0 & avg_log2FC_iso > 0 & avg_log2FC_KI2 > 0 ~ "Disease_down_PGRN_up",
      TRUE ~ "Other"
    )
  ) %>%
  filter(
    ligand_rescue_direction != "Other",
    ligand %in% all_ligands
  ) %>%
  left_join(
    mg_expr_stats,
    by = c("ligand" = "gene_symbol")
  ) %>%
  filter(
    !is.na(ligand_pct_expr_MG),
    ligand_pct_expr_MG >= pct_expr_cutoff
  ) %>%
  transmute(
    ligand,
    ligand_rescue_direction,
    ligand_rescue_pattern = rescue_pattern,
    ligand_rescue_lines = rescue_lines,
    ligand_disease_log2FC = avg_log2FC_p301k_vs_ct,
    ligand_iso_log2FC = avg_log2FC_iso,
    ligand_KI2_log2FC = avg_log2FC_KI2,
    ligand_rescue_mean_log2FC,
    ligand_padj_disease = p_val_adj_p301k_vs_ct,
    ligand_padj_iso = p_val_adj_iso,
    ligand_padj_KI2 = p_val_adj_KI2,
    ligand_pct_expr_MG,
    ligand_avg_expr_MG
  ) %>%
  distinct(ligand, .keep_all = TRUE)

cat("\nEligible MG ligands after both-line rescue, NicheNet ligand, and MG expression filters:\n")
print(eligible_mg_ligands %>% count(ligand_rescue_direction))

# ============================================================
# 7. Define reliable neuronal TF universe
# ============================================================

eligible_neuron_tfs <- tf_tbl %>%
  mutate(
    target_tf = gene_symbol,
    tf_rescue_mean_log2FC = rowMeans(
      cbind(avg_log2FC_iso, avg_log2FC_KI2),
      na.rm = TRUE
    ),
    tf_rescue_direction = case_when(
      avg_log2FC_p301k_vs_ct > 0 & avg_log2FC_iso < 0 & avg_log2FC_KI2 < 0 ~ "Disease_up_PGRN_down",
      avg_log2FC_p301k_vs_ct < 0 & avg_log2FC_iso > 0 & avg_log2FC_KI2 > 0 ~ "Disease_down_PGRN_up",
      TRUE ~ "Other"
    )
  ) %>%
  filter(
    tf_rescue_direction != "Other",
    target_tf %in% all_targets
  ) %>%
  left_join(
    neuron_expr_stats,
    by = c("target_tf" = "gene_symbol")
  ) %>%
  filter(
    !is.na(tf_pct_expr_neuron),
    tf_pct_expr_neuron >= pct_expr_cutoff
  ) %>%
  transmute(
    target_tf,
    tf_rescue_direction,
    tf_rescue_pattern = rescue_pattern,
    tf_rescue_lines = rescue_lines,
    tf_disease_log2FC = avg_log2FC_p301k_vs_ct,
    tf_iso_log2FC = avg_log2FC_iso,
    tf_KI2_log2FC = avg_log2FC_KI2,
    tf_rescue_mean_log2FC,
    tf_padj_disease = p_val_adj_p301k_vs_ct,
    tf_padj_iso = p_val_adj_iso,
    tf_padj_KI2 = p_val_adj_KI2,
    tf_pct_expr_neuron,
    tf_avg_expr_neuron
  ) %>%
  distinct(target_tf, .keep_all = TRUE)

cat("\nEligible neuronal TFs after FIMO/rescued DEG, NicheNet target, and neuron expression filters:\n")
print(eligible_neuron_tfs %>% count(tf_rescue_direction))

if (nrow(eligible_mg_ligands) == 0) {
  stop("No eligible MG ligands passed filters.")
}

if (nrow(eligible_neuron_tfs) == 0) {
  stop("No eligible neuronal TFs passed filters.")
}



# ============================================================
# 8. Extract ligand -> TF NicheNet regulatory potential weights
# Use a stronger relative cutoff instead of weight > 0
# ============================================================

weight_quantile_cutoff <- 0.95   # keep top 10% strongest nonzero weights
# Try 0.75 for broader, 0.90 for focused, 0.95 for very stringent

ligands_use <- eligible_mg_ligands$ligand
tfs_use <- eligible_neuron_tfs$target_tf

ligand_target_score_mat <- ligand_target_matrix_use[
  tfs_use,
  ligands_use,
  drop = FALSE
]

# First make all possible ligand -> TF pairs
ligand_tf_pair_tbl_raw <- as.data.frame(as.table(as.matrix(ligand_target_score_mat))) %>%
  transmute(
    target_tf = as.character(Var1),
    ligand = as.character(Var2),
    nichenet_ligand_target_weight = as.numeric(Freq)
  )

# Keep nonzero weights first
ligand_tf_pair_tbl_nonzero <- ligand_tf_pair_tbl_raw %>%
  filter(nichenet_ligand_target_weight > 0)

# Define data-driven cutoff
weight_cutoff_value <- quantile(
  ligand_tf_pair_tbl_nonzero$nichenet_ligand_target_weight,
  probs = weight_quantile_cutoff,
  na.rm = TRUE
)

cat("\nNicheNet weight cutoff quantile:", weight_quantile_cutoff, "\n")
cat("NicheNet weight cutoff value:", weight_cutoff_value, "\n")

cat("\nWeight distribution among nonzero links:\n")
print(
  quantile(
    ligand_tf_pair_tbl_nonzero$nichenet_ligand_target_weight,
    probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1),
    na.rm = TRUE
  )
)

# Apply stronger weight cutoff
ligand_tf_pair_tbl_all <- ligand_tf_pair_tbl_nonzero %>%
  filter(nichenet_ligand_target_weight >= weight_cutoff_value) %>%
  left_join(eligible_mg_ligands, by = "ligand") %>%
  left_join(eligible_neuron_tfs, by = "target_tf") %>%
  mutate(
    same_rescue_direction = ligand_rescue_direction == tf_rescue_direction
  ) %>%
  arrange(
    desc(nichenet_ligand_target_weight),
    ligand,
    target_tf
  )

cat("\nLigand -> TF links after weight cutoff:\n")
cat("Total pairs:", nrow(ligand_tf_pair_tbl_all), "\n")
cat("Unique ligands:", n_distinct(ligand_tf_pair_tbl_all$ligand), "\n")
cat("Unique TFs:", n_distinct(ligand_tf_pair_tbl_all$target_tf), "\n")

cat("\nLigands retained:\n")
print(sort(unique(ligand_tf_pair_tbl_all$ligand)))

cat("\nTFs retained:\n")
print(sort(unique(ligand_tf_pair_tbl_all$target_tf)))

if (nrow(ligand_tf_pair_tbl_all) == 0) {
  stop("No ligand -> TF pairs passed the NicheNet weight cutoff. Lower weight_quantile_cutoff.")
}



# ============================================================
# Ligand -> TF network plot
# Clean version:
#   - Uses ligand_tf_pair_tbl_all directly
#   - Keeps only plotted TFs with actual connections
#   - No disconnected TF labels
#   - Edge thickness = NicheNet regulatory potential weight
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(grid)
})

# ============================================================
# 0. User settings
# ============================================================

outdir <- "NicheNet_reliable_MG_ligand_to_neuron_TF_pairs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Use the filtered table from your current threshold result
plot_source_tbl <- ligand_tf_pair_tbl_all

# If too messy, use 3 to 5
top_links_per_ligand <- 5

# IMPORTANT:
# FALSE means only show TFs that actually have plotted edges
keep_all_final_tfs <- FALSE

# ============================================================
# 1. Build plotting edge table
# ============================================================

pair_plot_tbl <- plot_source_tbl %>%
  filter(
    !is.na(ligand),
    !is.na(target_tf),
    !is.na(nichenet_ligand_target_weight),
    nichenet_ligand_target_weight > 0
  ) %>%
  group_by(ligand) %>%
  slice_max(
    order_by = nichenet_ligand_target_weight,
    n = top_links_per_ligand,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    ligand,
    desc(nichenet_ligand_target_weight),
    target_tf
  )

if (nrow(pair_plot_tbl) == 0) {
  stop("No ligand -> TF pairs available for plotting.")
}

cat("\nPlotted edges:", nrow(pair_plot_tbl), "\n")
cat("Plotted ligands:", n_distinct(pair_plot_tbl$ligand), "\n")
cat("Plotted TFs:", n_distinct(pair_plot_tbl$target_tf), "\n")

# ============================================================
# 2. Define ligand order
# Prefer rescue direction, then strongest weight
# ============================================================

ligand_order_tbl <- pair_plot_tbl %>%
  group_by(ligand) %>%
  summarise(
    ligand_rescue_direction = first(ligand_rescue_direction),
    max_weight = max(nichenet_ligand_target_weight, na.rm = TRUE),
    mean_weight = mean(nichenet_ligand_target_weight, na.rm = TRUE),
    n_tfs = n_distinct(target_tf),
    .groups = "drop"
  ) %>%
  mutate(
    ligand_rescue_direction = factor(
      ligand_rescue_direction,
      levels = c("Disease_up_PGRN_down", "Disease_down_PGRN_up")
    )
  ) %>%
  arrange(
    ligand_rescue_direction,
    desc(max_weight),
    desc(n_tfs),
    ligand
  )

ligand_order <- ligand_order_tbl$ligand

# ============================================================
# 3. Define TF order
# Only use TFs that have plotted edges
# ============================================================

tf_order_tbl <- pair_plot_tbl %>%
  group_by(target_tf) %>%
  summarise(
    tf_rescue_direction = first(tf_rescue_direction),
    max_weight = max(nichenet_ligand_target_weight, na.rm = TRUE),
    mean_weight = mean(nichenet_ligand_target_weight, na.rm = TRUE),
    n_ligands = n_distinct(ligand),
    .groups = "drop"
  ) %>%
  mutate(
    tf_rescue_direction = factor(
      tf_rescue_direction,
      levels = c("Disease_up_PGRN_down", "Disease_down_PGRN_up")
    )
  ) %>%
  arrange(
    tf_rescue_direction,
    desc(max_weight),
    desc(n_ligands),
    target_tf
  )

tf_order <- tf_order_tbl$target_tf

# ============================================================
# 4. Build node positions
# ============================================================

ligand_nodes <- tibble(
  name = ligand_order,
  node_type = "Ligand",
  x = 0.28,
  y = rev(seq_along(ligand_order))
)

tf_nodes <- tibble(
  name = tf_order,
  node_type = "TF",
  x = 0.72,
  y = rev(seq_along(tf_order))
)

plot_height <- max(nrow(ligand_nodes), nrow(tf_nodes))

# Stretch shorter side to occupy similar vertical space
if (nrow(ligand_nodes) < plot_height) {
  ligand_nodes$y <- rev(seq(1, plot_height, length.out = nrow(ligand_nodes)))
}

if (nrow(tf_nodes) < plot_height) {
  tf_nodes$y <- rev(seq(1, plot_height, length.out = nrow(tf_nodes)))
}

# ============================================================
# 5. Join node coordinates to edges
# ============================================================

pair_plot_tbl2 <- pair_plot_tbl %>%
  left_join(
    ligand_nodes %>%
      select(ligand = name, x_ligand = x, y_ligand = y),
    by = "ligand"
  ) %>%
  left_join(
    tf_nodes %>%
      select(target_tf = name, x_tf = x, y_tf = y),
    by = "target_tf"
  ) %>%
  filter(
    !is.na(x_ligand),
    !is.na(y_ligand),
    !is.na(x_tf),
    !is.na(y_tf)
  )

if (nrow(pair_plot_tbl2) == 0) {
  stop("No plottable edges remained after joining node positions.")
}



# ============================================================
# MG ligand dot plot + ligand->TF network + neuron TF dot plot
# Uses pair_plot_tbl2 as the final ligand->TF edge table
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(cowplot)
  library(grid)
})

# ============================================================
# 0. Settings
# ============================================================

obj <- xeno_wKI_final_label

conditions_use <- c(
  "4R_P301S_Ctrl",
  "4R_P301S_K18",
  "ISOKI_4R_P301S_K18",
  "KI_4R_P301S_K18"
)

# outdir <- "NicheNet_MG_ligand_TF_main_figure"
# dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

dir_colors <- c(
  "Disease_up_PGRN_down" = "firebrick2",
  "Disease_down_PGRN_up" = "royalblue3"
)

assay_use <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)

# ============================================================
# 1. Build ordered ligand and TF tables from pair_plot_tbl2
#    Order: disease up / PGRN down first, then disease down / PGRN up
# ============================================================

ligand_info <- pair_plot_tbl2 %>%
  group_by(ligand) %>%
  summarise(
    direction = first(ligand_rescue_direction),
    max_weight = max(nichenet_ligand_target_weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    direction = factor(
      direction,
      levels = c("Disease_up_PGRN_down", "Disease_down_PGRN_up")
    )
  ) %>%
  arrange(direction, desc(max_weight), ligand)

tf_info <- pair_plot_tbl2 %>%
  group_by(target_tf) %>%
  summarise(
    direction = first(tf_rescue_direction),
    max_weight = max(nichenet_ligand_target_weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    direction = factor(
      direction,
      levels = c("Disease_up_PGRN_down", "Disease_down_PGRN_up")
    )
  ) %>%
  arrange(direction, desc(max_weight), target_tf)

ligand_order <- ligand_info$ligand
tf_order <- tf_info$target_tf

# ============================================================
# 2. Prepare MG and neuron Seurat objects
# ============================================================

obj$Condition <- factor(as.character(obj$Condition), levels = conditions_use)

mg_obj <- subset(
  obj,
  subset = control_celltype == "MG" & Condition %in% conditions_use
)

neu_obj <- subset(
  obj,
  subset = grepl("^Neu_", control_celltype) & Condition %in% conditions_use
)

DefaultAssay(mg_obj) <- assay_use
DefaultAssay(neu_obj) <- assay_use

# Keep only genes that are present
ligand_order_use <- ligand_order[ligand_order %in% rownames(mg_obj)]
tf_order_use <- tf_order[tf_order %in% rownames(neu_obj)]

ligand_info <- ligand_info %>% filter(ligand %in% ligand_order_use)
tf_info <- tf_info %>% filter(target_tf %in% tf_order_use)

# ============================================================
# 3. Helper functions
# ============================================================

make_dotplot_df <- function(seu, features_use, conditions_use) {
  dp <- DotPlot(
    seu,
    features = features_use,
    group.by = "Condition"
  )$data
  
  dp %>%
    mutate(
      id = factor(as.character(id), levels = conditions_use),
      features.plot = factor(as.character(features.plot), levels = rev(features_use))
    )
}

make_expression_dotplot <- function(dp_df, title_text, show_legend = TRUE) {
  ggplot(dp_df, aes(x = id, y = features.plot)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_size(range = c(1.5, 8), name = "% cells") +
    scale_color_gradient2(
      low = "blue",
      mid = "gray90",
      high = "firebrick3",
      midpoint = 0,
      name = "Avg expr\nscaled"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = title_text
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(face = "plain"),
      legend.position = if (show_legend) "right" else "none"
    )
}

make_direction_strip <- function(gene_df, gene_col, title_text = NULL) {
  strip_df <- gene_df %>%
    transmute(
      gene = .data[[gene_col]],
      direction = as.character(direction)
    ) %>%
    mutate(
      gene = factor(gene, levels = rev(unique(gene)))
    )
  
  ggplot(strip_df, aes(x = 1, y = gene, fill = direction)) +
    geom_tile(color = NA) +
    scale_fill_manual(values = dir_colors, drop = FALSE) +
    labs(title = title_text) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 0)
    )
}

# ============================================================
# 4. Make MG ligand dot plot
# ============================================================

mg_dp_df <- make_dotplot_df(
  seu = mg_obj,
  features_use = ligand_order_use,
  conditions_use = conditions_use
)

p_mg_dot <- make_expression_dotplot(
  dp_df = mg_dp_df,
  title_text = "MG ligands",
  show_legend = T
)+
  scale_y_discrete(position = "right") +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(3.5, "cm"),
      barheight = unit(0.35, "cm")
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1
    )
  ) +
  theme(
    axis.text.y.right = element_text(size = 11),
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    legend.position = "left",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )



print(p_mg_dot)





# ============================================================
# 5. Make neuron TF dot plot
# ============================================================

neu_dp_df <- make_dotplot_df(
  seu = neu_obj,
  features_use = tf_order_use,
  conditions_use = conditions_use
)

p_neu_dot <- make_expression_dotplot(
  dp_df = neu_dp_df,
  title_text = "Neuron TFs",
  show_legend = TRUE
) +
  scale_y_discrete(position = "right") +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(3.5, "cm"),
      barheight = unit(0.35, "cm")
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1
    )
  ) +
  theme(
    axis.text.y.right = element_text(size = 11),
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    legend.position = "left",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


print(p_neu_dot)

# ============================================================
# 6. Remake ligand -> TF network plot
# TFs on the left, ligands on the right
# Arrows point from MG ligands to neuronal TFs
# Edge thickness = NicheNet regulatory potential weight
# ============================================================

# ------------------------------------------------------------
# 6.1 Prepare edge table
# ------------------------------------------------------------

plot_edges <- pair_plot_tbl2 %>%
  mutate(
    ligand = as.character(ligand),
    target_tf = as.character(target_tf)
  ) %>%
  filter(
    ligand %in% ligand_order_use,
    target_tf %in% tf_order_use
  )

if (nrow(plot_edges) == 0) {
  stop("No ligand-TF links found for the selected ligand and TF orders.")
}

# ------------------------------------------------------------
# 6.2 Build node positions
# Left side = neuronal TFs
# Right side = MG ligands
# ------------------------------------------------------------

tf_nodes <- tibble(
  target_tf = as.character(tf_order_use),
  direction = as.character(tf_info$direction[match(tf_order_use, tf_info$target_tf)]),
  tf_x = 0.24,
  tf_y = rev(seq_along(tf_order_use))
)

ligand_nodes <- tibble(
  ligand = as.character(ligand_order_use),
  direction = as.character(ligand_info$direction[match(ligand_order_use, ligand_info$ligand)]),
  lig_x = 0.76,
  lig_y = rev(seq_along(ligand_order_use))
)

plot_height <- max(nrow(tf_nodes), nrow(ligand_nodes))

# Stretch shorter side so both columns fill similar vertical space
if (nrow(tf_nodes) < plot_height) {
  tf_nodes$tf_y <- rev(seq(1, plot_height, length.out = nrow(tf_nodes)))
}

if (nrow(ligand_nodes) < plot_height) {
  ligand_nodes$lig_y <- rev(seq(1, plot_height, length.out = nrow(ligand_nodes)))
}

# ------------------------------------------------------------
# 6.3 Join node coordinates to edges
# ------------------------------------------------------------

plot_edges2 <- plot_edges %>%
  left_join(
    ligand_nodes %>% select(ligand, lig_x, lig_y),
    by = "ligand"
  ) %>%
  left_join(
    tf_nodes %>% select(target_tf, tf_x, tf_y),
    by = "target_tf"
  ) %>%
  filter(
    !is.na(lig_x),
    !is.na(lig_y),
    !is.na(tf_x),
    !is.na(tf_y)
  )

cat("\nNetwork plot check:\n")
cat("Edges:", nrow(plot_edges2), "\n")
cat("Ligands:", n_distinct(plot_edges2$ligand), "\n")
cat("TFs:", n_distinct(plot_edges2$target_tf), "\n")

if (nrow(plot_edges2) == 0) {
  stop("No network edges remain after joining ligand and TF coordinates.")
}

# ------------------------------------------------------------
# 6.4 Plot
# Arrows point from ligand to TF:
#   start = ligand side
#   end   = TF side
# ------------------------------------------------------------

y_top <- plot_height + 2.6
y_sub <- plot_height + 1.5

p_net <- ggplot() +
  annotate(
    "rect",
    xmin = 0.02, xmax = 0.46,
    ymin = 0.3, ymax = plot_height + 0.7,
    fill = "grey92",
    color = NA
  ) +
  annotate(
    "rect",
    xmin = 0.54, xmax = 0.98,
    ymin = 0.3, ymax = plot_height + 0.7,
    fill = "grey92",
    color = NA
  ) +
  
  # arrows from MG ligands to neuronal TFs
  geom_segment(
    data = plot_edges2,
    aes(
      x = lig_x - 0.035,
      y = lig_y,
      xend = tf_x + 0.035,
      yend = tf_y,
      linewidth = nichenet_ligand_target_weight
    ),
    color = "black",
    alpha = .8,
    lineend = "round",
    arrow = grid::arrow(
      length = unit(0.08, "inches"),
      type = "closed"
    )
  ) +
  
  # TF nodes, left
  geom_point(
    data = tf_nodes,
    aes(x = tf_x, y = tf_y, fill = direction),
    shape = 22,
    size = 4.2,
    color = "black",
    stroke = 0.2
  ) +
  
  # ligand nodes, right
  geom_point(
    data = ligand_nodes,
    aes(x = lig_x, y = lig_y, fill = direction),
    shape = 22,
    size = 4.2,
    color = "black",
    stroke = 0.2
  ) +
  
  # TF labels, left
  geom_text(
    data = tf_nodes,
    aes(x = tf_x - 0.02, y = tf_y, label = target_tf),
    hjust = 1,
    size = 4.3
  ) +
  
  # ligand labels, right
  geom_text(
    data = ligand_nodes,
    aes(x = lig_x + 0.02, y = lig_y, label = ligand),
    hjust = 0,
    size = 4.3
  ) +
  
  # headers
  annotate(
    "text",
    x = 0.24, y = y_top,
    label = "TFs",
    fontface = "bold",
    size = 7
  ) +
  annotate(
    "text",
    x = 0.24, y = y_sub,
    label = "in neurons",
    size = 5.8
  ) +
  annotate(
    "text",
    x = 0.76, y = y_top,
    label = "Ligands",
    fontface = "bold",
    size = 7
  ) +
  annotate(
    "text",
    x = 0.76, y = y_sub,
    label = "in microglia",
    size = 5.8
  ) +
  
  scale_fill_manual(values = dir_colors, guide = "none") +
  scale_linewidth_continuous(
    range = c(0.3, 1),
    guide = "none"
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, plot_height + 3.2),
    clip = "off"
  ) +
  theme_void(base_size = 14) +
  ggtitle("Ligand → TF predicted links") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(10, 70, 10, 70)
  )

print(p_net)






fimo <- read.delim(
  "G:/Other computers/My Laptop/sequencing analysis/MS/R code/fimo.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


fimo <- read_tsv(fimo_file, comment = "#")

head(fimo)

fimo <- as.data.frame(fimo)



#Fig 6A use unbiasedly selected neuron TFs and MG ligands to annotate HNRNPK-AS1 promoter

# ============================================================
# HNRNPK-AS1 promoter FIMO motif map
# Highlight TFs predicted downstream of plotted MG ligands
#
# Logic:
#   1. Read FIMO output
#   2. Use promoter_seq to define promoter length
#   3. Define all rescued neuronal TFs from your FIMO-overlap TF table
#   4. Define highlighted TFs as TFs connected to plotted MG ligands
#   5. Plot all rescued TF motif hits in black
#   6. Highlight ligand-linked TF motif hits in red
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

# # ============================================================
# # 0. User settings
# # ============================================================
# 
# outdir <- "HNRNPK_AS1_FIMO_ligand_linked_TF_map"
# dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
# 
# # Use your FIMO file path
# fimo_file <- "G:/Other computers/My Laptop/sequencing analysis/MS/R code/fimo.tsv"
# 
# # Use either pair_plot_tbl2 if you want the exact plotted network TFs,
# # or ligand_tf_pair_tbl_all if you want all TFs linked after NicheNet filtering.
use_exact_plotted_network <- TRUE

# ============================================================
# 1. Clean promoter sequence and get length
# ============================================================

promoter_seq_clean <- promoter_seq |>
  toupper() |>
  gsub("[^ATCGN]", "", x = _)

promoter_length <- nchar(promoter_seq_clean)

cat("Promoter length:", promoter_length, "bp\n")

# ============================================================
# 2. Read FIMO
# ============================================================


fimo <- as.data.frame(fimo)

# Standardize p/q value column names
if ("p-value" %in% colnames(fimo)) {
  fimo <- fimo %>% rename(p.value = `p-value`)
}
if ("q-value" %in% colnames(fimo)) {
  fimo <- fimo %>% rename(q.value = `q-value`)
}

if (!all(c("motif_alt_id", "start", "stop", "strand", "p.value") %in% colnames(fimo))) {
  stop("FIMO table is missing one or more required columns: motif_alt_id, start, stop, strand, p.value")
}

# ============================================================
# 3. Get all rescued neuronal TFs from your FIMO overlap table
# ============================================================

tf_col <- if ("gene" %in% colnames(NEU_rescued_DEGs_overlap_FIMO_motif_alt_id_split_complex_motifs)) {
  "gene"
} else if ("genes" %in% colnames(NEU_rescued_DEGs_overlap_FIMO_motif_alt_id_split_complex_motifs)) {
  "genes"
} else {
  stop("Could not find gene or genes column in NEU_rescued_DEGs_overlap_FIMO_motif_alt_id_split_complex_motifs")
}

rescued_tf_list <- NEU_rescued_DEGs_overlap_FIMO_motif_alt_id_split_complex_motifs %>%
  mutate(tf = toupper(str_trim(as.character(.data[[tf_col]])))) %>%
  filter(!is.na(tf), tf != "") %>%
  pull(tf) %>%
  unique()

cat("Number of rescued neuronal TFs:", length(rescued_tf_list), "\n")

# ============================================================
# 4. Define TFs linked to plotted MG ligands
# ============================================================

if (use_exact_plotted_network) {
  
  # This highlights only TFs in the final network plot
  ligand_linked_tfs <- pair_plot_tbl2 %>%
    mutate(target_tf = toupper(str_trim(as.character(target_tf)))) %>%
    filter(!is.na(target_tf), target_tf != "") %>%
    pull(target_tf) %>%
    unique()
  
  plotted_ligands <- pair_plot_tbl2 %>%
    mutate(ligand = toupper(str_trim(as.character(ligand)))) %>%
    filter(!is.na(ligand), ligand != "") %>%
    pull(ligand) %>%
    unique()
  
} else {
  
  # This highlights all TFs linked to all filtered MG ligands
  ligand_linked_tfs <- ligand_tf_pair_tbl_all %>%
    mutate(target_tf = toupper(str_trim(as.character(target_tf)))) %>%
    filter(!is.na(target_tf), target_tf != "") %>%
    pull(target_tf) %>%
    unique()
  
  plotted_ligands <- ligand_tf_pair_tbl_all %>%
    mutate(ligand = toupper(str_trim(as.character(ligand)))) %>%
    filter(!is.na(ligand), ligand != "") %>%
    pull(ligand) %>%
    unique()
}

# Keep only TFs that are also in your rescued TF/FIMO-overlap list
ligand_linked_tfs <- intersect(ligand_linked_tfs, rescued_tf_list)

cat("Number of plotted or filtered MG ligands:", length(plotted_ligands), "\n")
cat("Number of ligand-linked TFs to highlight:", length(ligand_linked_tfs), "\n")
print(ligand_linked_tfs)

# ============================================================
# 5. Prepare FIMO table
# ============================================================

fimo2 <- fimo %>%
  mutate(
    motif_alt_id = as.character(motif_alt_id),
    motif_alt_upper = toupper(str_trim(motif_alt_id))
  ) %>%
  separate_rows(motif_alt_upper, sep = "::|,|;|/") %>%
  rename(motif_tf = motif_alt_upper) %>%
  mutate(
    motif_tf = str_trim(motif_tf),
    motif_center = (start + stop) / 2,
    neglog10p = -log10(p.value)
  )

# ============================================================
# 6. All rescued TF motif hits
# These are black dots
# ============================================================

all_rescued_hits <- fimo2 %>%
  filter(motif_tf %in% rescued_tf_list) %>%
  distinct(motif_tf, start, stop, strand, .keep_all = TRUE) %>%
  arrange(start, p.value)

cat("All rescued TF motif hits:", nrow(all_rescued_hits), "\n")
cat("Unique rescued TFs with FIMO hits:", n_distinct(all_rescued_hits$motif_tf), "\n")

# ============================================================
# 7. Ligand-linked TF motif hits
# These are highlighted red dots
# ============================================================

ligand_linked_hits_all <- all_rescued_hits %>%
  filter(motif_tf %in% ligand_linked_tfs) %>%
  arrange(start, p.value)

cat("Ligand-linked TF motif hits:", nrow(ligand_linked_hits_all), "\n")
cat("Unique ligand-linked TFs with FIMO hits:", n_distinct(ligand_linked_hits_all$motif_tf), "\n")

# Pick one best hit per highlighted TF for labeling
ligand_linked_hits_labeled <- ligand_linked_hits_all %>%
  group_by(motif_tf) %>%
  slice_min(order_by = p.value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(start) %>%
  mutate(
    label_y = rep(c(1.08, 1.12, 1.16, 1.20), length.out = n())
  )


# ============================================================
# 9. Accurate vertical FIMO motif map
# Red dots = labeled best motif hit per ligand-linked TF only
# Black dots = all rescued TF motif hits
# Promoter position 1 is shown at the top
# ============================================================

# Use the original promoter coordinate for plotting.
# The y-axis is reversed below so position 1 appears at the top.
all_rescued_hits_v <- all_rescued_hits %>%
  mutate(
    motif_y = motif_center
  )

# Pick ONE best motif hit per highlighted TF.
# These are the only red dots, so every red dot has a label.
ligand_linked_hits_labeled_v <- ligand_linked_hits_all %>%
  group_by(motif_tf) %>%
  slice_min(order_by = p.value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(motif_center) %>%
  mutate(
    motif_y = motif_center,
    
    # stagger labels slightly to the right to reduce overlap
    label_x = rep(c(1.05, 1.1, 1.15), length.out = n())
  )




# Check that every red dot will be labeled
cat("Number of red dots:", nrow(ligand_linked_hits_labeled_v), "\n")
cat("Number of labels:", nrow(ligand_linked_hits_labeled_v), "\n")
print(ligand_linked_hits_labeled_v %>% select(motif_tf, start, stop, motif_center, p.value))

p_motif_map_ligand_linked_vertical <- ggplot() +
  
  # promoter backbone
  geom_segment(
    aes(x = 1, xend = 1, y = 1, yend = promoter_length),
    linewidth = 0.8,
    color = "gray50"
  ) +
  
  # all rescued TF motif hits as black dots
  geom_point(
    data = all_rescued_hits_v,
    aes(x = 1, y = motif_y),
    color = "black",
    size = 2,
    alpha = 0.65
  ) +
  
  # only labeled best ligand-linked TF motif hits as red dots
  geom_point(
    data = ligand_linked_hits_labeled_v,
    aes(x = 1, y = motif_y),
    color = "firebrick3",
    size = 3.2,
    alpha = 0.98
  ) +
  
  # horizontal connector from each red dot to its own label
  geom_segment(
    data = ligand_linked_hits_labeled_v,
    aes(
      x = 1.02,
      xend = label_x - 0.02,
      y = motif_y,
      yend = motif_y
    ),
    color = "firebrick3",
    linewidth = 0.35,
    alpha = 0.85
  ) +
  
  # labels aligned with the exact red dot position
  geom_text(
    data = ligand_linked_hits_labeled_v,
    aes(
      x = label_x,
      y = motif_y,
      label = motif_tf
    ),
    color = "black",
    size = 4.2,
    hjust = 0,
    vjust = 0.5
  ) +
  
  # Reverse y-axis so promoter position 1 is at the top
  scale_y_reverse(
    limits = c(promoter_length, 1),
    breaks = seq(0, promoter_length, by = 200),
    expand = c(0.01, 0.01)
  ) +
  
  scale_x_continuous(
    limits = c(0.94, 1.65),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  
  labs(
    x = NULL,
    y = "Position in HNRNPK-AS1 promoter sequence (bp)"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_motif_map_ligand_linked_vertical)

# ============================================================
# 9. Horizontal FIMO motif map
# Clean version:
#   black dots = all rescued TF motif hits
#   red dots   = one best motif hit per highlighted TF
#   every red dot has a label
# ============================================================

# Rebuild the labeled table so it is explicit and consistent
ligand_linked_hits_labeled <- ligand_linked_hits_all %>%
  group_by(motif_tf) %>%
  slice_min(order_by = p.value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(motif_center) %>%
  mutate(
    label_y = rep(c(1.08, 1.12, 1.16, 1.20), length.out = n())
  )

cat("Number of red dots:", nrow(ligand_linked_hits_labeled), "\n")
cat("Number of labels:", nrow(ligand_linked_hits_labeled), "\n")

p_motif_map_ligand_linked <- ggplot() +
  
  # promoter backbone
  geom_segment(
    aes(x = 1, xend = promoter_length, y = 1, yend = 1),
    linewidth = 0.8,
    color = "gray50"
  ) +
  
  # all rescued TF motif hits
  geom_point(
    data = all_rescued_hits,
    aes(x = motif_center, y = 1),
    color = "black",
    size = 2,
    alpha = 0.65
  ) +
  
  # ONLY the labeled best ligand-linked TF hits
  geom_point(
    data = ligand_linked_hits_labeled,
    aes(x = motif_center, y = 1),
    color = "firebrick3",
    size = 3.0,
    alpha = 0.95
  ) +
  
  # vertical lines from each selected red dot to its label
  geom_segment(
    data = ligand_linked_hits_labeled,
    aes(
      x = motif_center,
      xend = motif_center,
      y = 1.01,
      yend = label_y - 0.012
    ),
    color = "firebrick3",
    linewidth = 0.35,
    alpha = 0.8
  ) +
  
  # labels centered above their own red dots
  geom_text(
    data = ligand_linked_hits_labeled,
    aes(
      x = motif_center,
      y = label_y,
      label = motif_tf
    ),
    color = "black",
    size = 4.4,
    hjust = 0.5,
    vjust = 0
  ) +
  
  scale_x_continuous(
    limits = c(1, promoter_length),
    breaks = seq(0, promoter_length, by = 200),
    expand = c(0.01, 0.01)
  ) +
  
  scale_y_continuous(
    limits = c(0.96, 1.24),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Position in HNRNPK-AS1 promoter sequence (bp)",
    y = NULL,
    title = "Ligand-linked rescued TF motifs in the HNRNPK-AS1 promoter"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_motif_map_ligand_linked)

