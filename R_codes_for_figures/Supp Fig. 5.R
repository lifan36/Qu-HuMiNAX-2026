############################################################
# Process the GSE168408 / Herring et al. human PFC reference
#
# Output:
#   GSE168408_reference_EN_IN_MG_AST.rds
#   GSE168408_reference_cell_counts.csv
#   GSE168408_reference_metadata.csv
############################################################

rm(list = ls())
gc()

############################################################
# 1. Packages
############################################################

cran_packages <- c(
  "Seurat",
  "Matrix",
  "dplyr",
  "curl"
)

cran_missing <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(cran_missing) > 0) {
  install.packages(cran_missing)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "zellkonverter",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "DelayedArray"
)

bioc_missing <- bioc_packages[
  !vapply(
    bioc_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(bioc_missing) > 0) {
  BiocManager::install(
    bioc_missing,
    ask = FALSE,
    update = FALSE
  )
}

library(Seurat)
library(Matrix)
library(dplyr)
library(zellkonverter)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(DelayedArray)

set.seed(1234)

############################################################
# 2. Directories
############################################################

reference_dir <- paste0(
  "C:/Users/wenhq/Documents/",
  "compare cell mature/reference 2 human"
)

output_dir <- file.path(
  reference_dir,
  "processed_reference"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# 3. Download the authors' processed object
#
# This contains:
#   - published cell annotations
#   - corrected developmental-stage metadata
#   - raw counts in X
#   - downsampled CPM in ds_norm_cts
############################################################

h5ad_file <- file.path(
  reference_dir,
  "RNA-all_full-counts-and-downsampled-CPM.h5ad"
)

h5ad_url <- paste0(
  "https://storage.googleapis.com/neuro-dev/",
  "Processed_data/",
  "RNA-all_full-counts-and-downsampled-CPM.h5ad"
)

if (!file.exists(h5ad_file)) {
  
  message("Downloading the official processed H5AD file...")
  
  curl::curl_download(
    url = h5ad_url,
    destfile = h5ad_file,
    quiet = FALSE,
    mode = "wb"
  )
}

stopifnot(file.exists(h5ad_file))

############################################################
# 4. Read the H5AD as a SingleCellExperiment
#
# use_hdf5 = TRUE reduces initial memory use.
############################################################

message("Reading H5AD file...")

reference_sce <- zellkonverter::readH5AD(
  file = h5ad_file,
  use_hdf5 = TRUE
)

cat("\nAvailable assays:\n")
print(SummarizedExperiment::assayNames(reference_sce))

reference_metadata <- as.data.frame(
  SummarizedExperiment::colData(reference_sce)
)

rownames(reference_metadata) <- colnames(reference_sce)

cat("\nAvailable metadata columns:\n")
print(colnames(reference_metadata))



unique(reference_metadata$cell_type)
############################################################
# Create annotation mapping robustly
############################################################

published_annotation_vector <- reference_metadata[[annotation_column]]

broad_celltype_vector <- reference_metadata[["broad_celltype"]]

# Convert possible list/DataFrame columns into simple vectors
published_annotation_vector <- as.character(
  unlist(
    published_annotation_vector,
    use.names = FALSE
  )
)

broad_celltype_vector <- as.character(
  unlist(
    broad_celltype_vector,
    use.names = FALSE
  )
)

# Confirm that both vectors correspond to the same nuclei
stopifnot(
  length(published_annotation_vector) ==
    length(broad_celltype_vector)
)

annotation_mapping <- data.frame(
  published_annotation = published_annotation_vector,
  broad_celltype = broad_celltype_vector,
  stringsAsFactors = FALSE
) %>%
  dplyr::group_by(
    published_annotation,
    broad_celltype
  ) %>%
  dplyr::summarise(
    n_nuclei = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    broad_celltype,
    dplyr::desc(n_nuclei)
  )

print(
  annotation_mapping,
  n = Inf
)












############################################################
# Correct GSE168408 broad-cell-type mapping
############################################################

# Extract the published annotations as a simple character vector
published_labels <- as.character(
  unlist(
    reference_metadata[[annotation_column]],
    use.names = FALSE
  )
)

# Exact mapping based on the published major cell classes
reference_metadata$broad_celltype <- dplyr::case_when(
  
  # Excitatory / projection neurons
  published_labels %in% c(
    "L2-3_CUX2",
    "L4_RORB",
    "L5-6_TLE4",
    "L5-6_THEMIS",
    "PN_dev"
  ) ~ "EN",
  
  # Inhibitory interneurons
  published_labels %in% c(
    "SST",
    "VIP",
    "ID2",
    "PV",
    "MGE_dev",
    "CGE_dev",
    "PV_SCUBE3",
    "LAMP5_NOS1"
  ) ~ "IN",
  
  # Astrocytes
  published_labels == "Astro" ~ "AST",
  
  # Microglia
  published_labels == "Micro" ~ "MG",
  
  # Exclude oligodendrocytes, OPCs, vascular cells,
  # poor-quality nuclei, and any other populations
  TRUE ~ NA_character_
)

############################################################
# Check the annotation mapping
############################################################

annotation_mapping <- data.frame(
  published_annotation = published_labels,
  broad_celltype = reference_metadata$broad_celltype,
  stringsAsFactors = FALSE
) %>%
  dplyr::group_by(
    published_annotation,
    broad_celltype
  ) %>%
  dplyr::summarise(
    n_nuclei = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    broad_celltype,
    dplyr::desc(n_nuclei)
  )

print(
  annotation_mapping,
  n = Inf
)

cat("\nBroad cell-type counts before subsetting:\n")

print(
  table(
    reference_metadata$broad_celltype,
    useNA = "ifany"
  )
)

############################################################
# Add corrected annotations to the SCE
############################################################

SummarizedExperiment::colData(
  reference_sce
)$published_annotation <- published_labels

SummarizedExperiment::colData(
  reference_sce
)$broad_celltype <-
  reference_metadata$broad_celltype

############################################################
# Retain EN, IN, AST, and MG
############################################################

keep_nuclei <- SummarizedExperiment::colData(
  reference_sce
)$broad_celltype %in% c(
  "EN",
  "IN",
  "AST",
  "MG"
)

reference_sce_4types <- reference_sce[
  ,
  keep_nuclei
]

############################################################
# Confirm retained cell counts
############################################################

cat("\nRetained EN, IN, AST, and MG nuclei:\n")

print(
  table(
    SummarizedExperiment::colData(
      reference_sce_4types
    )$broad_celltype
  )
)

cat("\nPublished subtypes retained within each broad class:\n")

print(
  table(
    SummarizedExperiment::colData(
      reference_sce_4types
    )$published_annotation,
    SummarizedExperiment::colData(
      reference_sce_4types
    )$broad_celltype
  )
)















############################################################
# 6. Standardize donor/sample and developmental age metadata
############################################################

first_existing_column <- function(
    metadata,
    possible_columns
) {
  
  found <- possible_columns[
    possible_columns %in% colnames(metadata)
  ]
  
  if (length(found) == 0) {
    return(NA_character_)
  }
  
  found[1]
}

sample_column <- first_existing_column(
  reference_metadata,
  c(
    "batch",
    "sample",
    "sample_id",
    "donor",
    "donor_id",
    "orig.ident"
  )
)

age_column <- first_existing_column(
  reference_metadata,
  c(
    "age",
    "age_id",
    "age_label",
    "developmental_age",
    "numerical_age"
  )
)

stage_column <- first_existing_column(
  reference_metadata,
  c(
    "stage_id",
    "stage",
    "developmental_stage",
    "age_group"
  )
)

if (!is.na(sample_column)) {
  
  reference_metadata$reference_sample <-
    as.character(
      reference_metadata[[sample_column]]
    )
  
} else {
  
  reference_metadata$reference_sample <-
    "GSE168408"
}

if (!is.na(age_column)) {
  
  reference_metadata$reference_age <-
    as.character(
      reference_metadata[[age_column]]
    )
  
} else if (!is.na(stage_column)) {
  
  reference_metadata$reference_age <-
    as.character(
      reference_metadata[[stage_column]]
    )
  
} else {
  
  stop(
    "Could not locate an age or developmental-stage column."
  )
}

if (!is.na(stage_column)) {
  
  reference_metadata$reference_stage <-
    as.character(
      reference_metadata[[stage_column]]
    )
  
} else {
  
  reference_metadata$reference_stage <-
    reference_metadata$reference_age
}

############################################################
# 7. Create a chronological age variable
#
# Prenatal gestational ages are negative.
# Postnatal days are converted to years.
# Ages labeled yr are retained as years.
############################################################

parse_age_to_years <- function(x) {
  
  x_clean <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  age_numeric <- rep(
    NA_real_,
    length(x_clean)
  )
  
  # Extract the numeric portion
  age_value <- suppressWarnings(
    as.numeric(
      gsub(
        pattern = "[^0-9.]",
        replacement = "",
        x = x_clean
      )
    )
  )
  
  ##########################################################
  # Prenatal gestational age
  #
  # ga22 means gestational week 22.
  # Convert relative to an approximate 40-week birth:
  # ga22 = -18 weeks before birth
  ##########################################################
  
  prenatal_index <- grepl(
    pattern = "^(ga|gw|pcw)",
    x = x_clean
  )
  
  age_numeric[prenatal_index] <- (
    age_value[prenatal_index] - 40
  ) / 52.1429
  
  ##########################################################
  # Postnatal days
  ##########################################################
  
  day_index <- grepl(
    pattern = "^[0-9.]+d$",
    x = x_clean
  )
  
  age_numeric[day_index] <-
    age_value[day_index] / 365.25
  
  ##########################################################
  # Postnatal months
  ##########################################################
  
  month_index <- grepl(
    pattern = "^[0-9.]+(m|mo|mos|month|months)$",
    x = x_clean
  )
  
  age_numeric[month_index] <-
    age_value[month_index] / 12
  
  ##########################################################
  # Postnatal years
  ##########################################################
  
  year_index <- grepl(
    pattern = "^[0-9.]+(y|yr|yrs|year|years)$",
    x = x_clean
  )
  
  age_numeric[year_index] <-
    age_value[year_index]
  
  age_numeric
}

reference_metadata$reference_age_numeric <-
  parse_age_to_years(
    reference_metadata$reference_age
  )

############################################################
# Check the resulting chronological order
############################################################

age_check <- reference_metadata %>%
  dplyr::select(
    reference_age,
    reference_age_numeric
  ) %>%
  dplyr::distinct() %>%
  dplyr::arrange(
    reference_age_numeric
  )

print(
  age_check
)


############################################################
# 8. Add standardized metadata back to the SCE
############################################################

SummarizedExperiment::colData(
  reference_sce
)$broad_celltype <-
  reference_metadata$broad_celltype

SummarizedExperiment::colData(
  reference_sce
)$reference_sample <-
  reference_metadata$reference_sample

SummarizedExperiment::colData(
  reference_sce
)$reference_age <-
  reference_metadata$reference_age

SummarizedExperiment::colData(
  reference_sce
)$reference_stage <-
  reference_metadata$reference_stage

SummarizedExperiment::colData(
  reference_sce
)$reference_age_numeric <-
  reference_metadata$reference_age_numeric

SummarizedExperiment::colData(
  reference_sce
)$published_annotation <-
  as.character(
    reference_metadata[[annotation_column]]
  )

############################################################
# 9. Retain EN, IN, MG, and AST only
############################################################

keep_nuclei <- !is.na(
  SummarizedExperiment::colData(
    reference_sce
  )$broad_celltype
)

reference_sce_4types <- reference_sce[
  ,
  keep_nuclei
]

cat("\nNumber of retained nuclei:\n")

print(
  table(
    SummarizedExperiment::colData(
      reference_sce_4types
    )$broad_celltype
  )
)

############################################################
# 10. Identify the raw-count assay
############################################################

available_assays <- SummarizedExperiment::assayNames(
  reference_sce_4types
)

preferred_count_assays <- c(
  "X",
  "counts",
  "raw_counts",
  "raw"
)

count_assay <- preferred_count_assays[
  preferred_count_assays %in%
    available_assays
][1]

if (
  length(count_assay) == 0 ||
  is.na(count_assay)
) {
  
  count_assay <- available_assays[1]
  
  warning(
    "Could not identify an assay explicitly named X/counts. ",
    "Using assay: ",
    count_assay
  )
}

message(
  "Raw-count assay selected: ",
  count_assay
)

############################################################
# 11. Convert the selected counts to a sparse matrix
############################################################

reference_counts <- SummarizedExperiment::assay(
  reference_sce_4types,
  count_assay
)

if (!inherits(reference_counts, "dgCMatrix")) {
  
  message(
    "Converting the HDF5-backed matrix to dgCMatrix. ",
    "This step may use substantial RAM."
  )
  
  reference_counts <- as(
    reference_counts,
    "dgCMatrix"
  )
}

reference_counts <- Matrix::drop0(
  reference_counts
)

rownames(reference_counts) <- make.unique(
  rownames(reference_counts)
)

reference_metadata_4types <- as.data.frame(
  SummarizedExperiment::colData(
    reference_sce_4types
  )
)

rownames(reference_metadata_4types) <-
  colnames(reference_sce_4types)

############################################################
# 12. Construct the Seurat reference object
#
# The official H5AD is already QC filtered, so no additional
# nFeature or mitochondrial filtering is applied here.
############################################################

reference <- CreateSeuratObject(
  counts = reference_counts,
  meta.data = reference_metadata_4types,
  project = "GSE168408_Herring_PFC",
  assay = "RNA",
  min.cells = 3,
  min.features = 0
)

reference$dataset <- "GSE168408"
reference$reference_name <- "Herring_human_PFC"

reference$broad_celltype <- factor(
  reference$broad_celltype,
  levels = c(
    "EN",
    "IN",
    "AST",
    "MG"
  )
)

############################################################
# 13. Order the exact age labels chronologically
############################################################

age_order <- reference@meta.data %>%
  select(
    reference_age,
    reference_age_numeric
  ) %>%
  distinct() %>%
  arrange(
    reference_age_numeric,
    reference_age
  ) %>%
  pull(reference_age)

reference$reference_age <- factor(
  reference$reference_age,
  levels = unique(age_order),
  ordered = TRUE
)

############################################################
# 14. Standard log normalization
#
# Raw counts remain in the counts layer.
# Log-normalized expression is placed in the data layer.
############################################################

reference <- NormalizeData(
  object = reference,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

############################################################
# 15. Basic validation
############################################################

cat("\nFinal cell-type counts:\n")
print(
  table(
    reference$broad_celltype
  )
)

cat("\nFinal age-by-cell-type counts:\n")
print(
  table(
    reference$reference_age,
    reference$broad_celltype
  )
)

cat("\nDevelopmental stages:\n")
print(
  table(
    reference$reference_stage,
    reference$broad_celltype
  )
)

cat("\nSamples retained:\n")
print(
  table(
    reference$reference_sample,
    reference$broad_celltype
  )
)

############################################################
# 16. Save cell-count summary
############################################################

reference_cell_count_metadata <- data.frame(
  reference_age = as.character(
    unlist(
      reference@meta.data$reference_age,
      use.names = FALSE
    )
  ),
  
  reference_age_numeric = as.numeric(
    unlist(
      reference@meta.data$reference_age_numeric,
      use.names = FALSE
    )
  ),
  
  reference_stage = as.character(
    unlist(
      reference@meta.data$reference_stage,
      use.names = FALSE
    )
  ),
  
  reference_sample = as.character(
    unlist(
      reference@meta.data$reference_sample,
      use.names = FALSE
    )
  ),
  
  broad_celltype = as.character(
    unlist(
      reference@meta.data$broad_celltype,
      use.names = FALSE
    )
  ),
  
  stringsAsFactors = FALSE
)

# Confirm all columns have one value per nucleus
stopifnot(
  nrow(reference_cell_count_metadata) ==
    ncol(reference)
)

reference_cell_counts <- reference_cell_count_metadata %>%
  dplyr::group_by(
    reference_age,
    reference_age_numeric,
    reference_stage,
    reference_sample,
    broad_celltype
  ) %>%
  dplyr::summarise(
    n_nuclei = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    reference_age_numeric,
    broad_celltype,
    reference_sample
  )

print(
  reference_cell_counts,
  n = Inf
)
# 
# write.csv(
#   reference_cell_counts,
#   file = file.path(
#     output_dir,
#     "GSE168408_reference_cell_counts.csv"
#   ),
#   row.names = FALSE
# )

############################################################
# 17. Save standardized metadata
############################################################

metadata_to_save <- reference@meta.data %>%
  mutate(
    barcode = rownames(reference@meta.data)
  ) %>%
  select(
    barcode,
    dataset,
    reference_name,
    reference_sample,
    reference_age,
    reference_age_numeric,
    reference_stage,
    broad_celltype,
    published_annotation,
    everything()
  )

# write.csv(
#   metadata_to_save,
#   file = file.path(
#     output_dir,
#     "GSE168408_reference_metadata.csv"
#   ),
#   row.names = FALSE
# )

############################################################
# 18. Save the processed reference object
############################################################

saveRDS(
  reference,
  file = file.path(
    output_dir,
    "GSE168408_reference_EN_IN_MG_AST.rds"
  ),
  compress = FALSE
)

############################################################
# 19. Final object summary
############################################################

reference

cat(
  "\nSaved processed reference to:\n",
  file.path(
    output_dir,
    "GSE168408_reference_EN_IN_MG_AST.rds"
  ),
  "\n"
)


reference <-GSE168408_reference_EN_IN_MG_AST 







############################################################
# GSE168408 broad developmental-stage expression signatures
#
# Input:
#   reference
#
# Output:
#   reference_expression
#
# Columns:
#   gene
#   Prenatal_EN, Prenatal_IN, Prenatal_AST, Prenatal_MG
#   Neonatal_EN, ...
#   Adult_MG
#
# This matches the HuMiNAX workflow:
# average linear normalized expression across nuclei.
############################################################

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tibble)
library(readr)

############################################################
# 1. Use RNA assay
############################################################

DefaultAssay(reference) <- "RNA"

############################################################
# 2. Join Seurat v5 layers if necessary
############################################################

if (inherits(reference[["RNA"]], "Assay5")) {
  
  data_layers <- SeuratObject::Layers(
    reference[["RNA"]],
    search = "^data"
  )
  
  if (length(data_layers) > 1) {
    
    message(
      "Multiple RNA data layers detected. Joining layers."
    )
    
    reference <- SeuratObject::JoinLayers(
      reference,
      assay = "RNA"
    )
  }
  
  data_layers <- SeuratObject::Layers(
    reference[["RNA"]],
    search = "^data$"
  )
  
  if (length(data_layers) == 0) {
    
    message(
      "No normalized RNA data layer found. Running NormalizeData."
    )
    
    reference <- Seurat::NormalizeData(
      reference,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }
}

############################################################
# 3. Validate required metadata
############################################################

required_metadata <- c(
  "reference_age",
  "broad_celltype"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(reference@meta.data)
)

if (length(missing_metadata) > 0) {
  
  stop(
    "Missing required metadata columns: ",
    paste(
      missing_metadata,
      collapse = ", "
    )
  )
}

############################################################
# 4. Create broader developmental-age brackets
############################################################

reference$reference_age_character <- as.character(
  reference$reference_age
)

reference$age_bracket <- dplyr::case_when(
  
  reference$reference_age_character %in%
    c(
      "ga22",
      "ga24",
      "ga34"
    ) ~ "Prenatal",
  
  reference$reference_age_character %in%
    c(
      "2d",
      "34d"
    ) ~ "Neonatal",
  
  reference$reference_age_character %in%
    c(
      "86d",
      "118d",
      "179d"
    ) ~ "Early_infancy",
  
  reference$reference_age_character %in%
    c(
      "301d",
      "422d",
      "627d"
    ) ~ "Late_infancy",
  
  reference$reference_age_character %in%
    c(
      "2yr",
      "3yr",
      "4yr"
    ) ~ "Early_childhood",
  
  reference$reference_age_character %in%
    c(
      "6yr",
      "8yr",
      "10yr"
    ) ~ "Childhood",
  
  reference$reference_age_character %in%
    c(
      "12yr",
      "14yr",
      "16yr",
      "17yr"
    ) ~ "Adolescence",
  
  reference$reference_age_character %in%
    c(
      "20yr",
      "25yr",
      "40yr"
    ) ~ "Adult",
  
  TRUE ~ NA_character_
)

############################################################
# 5. Set chronological bracket order
############################################################

age_bracket_order <- c(
  "Prenatal",
  "Neonatal",
  "Early_infancy",
  "Late_infancy",
  "Early_childhood",
  "Childhood",
  "Adolescence",
  "Adult"
)

reference$age_bracket <- factor(
  reference$age_bracket,
  levels = age_bracket_order,
  ordered = TRUE
)

reference$broad_celltype <- factor(
  as.character(reference$broad_celltype),
  levels = c(
    "EN",
    "IN",
    "AST",
    "MG"
  )
)

############################################################
# 6. Confirm exact age-to-bracket mapping
############################################################

age_mapping <- data.frame(
  reference_age = as.character(
    reference$reference_age
  ),
  age_bracket = as.character(
    reference$age_bracket
  ),
  stringsAsFactors = FALSE
) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    age_bracket = factor(
      age_bracket,
      levels = age_bracket_order,
      ordered = TRUE
    )
  ) %>%
  dplyr::arrange(
    age_bracket
  )

print(
  age_mapping
)

unassigned_ages <- unique(
  as.character(
    reference$reference_age[
      is.na(reference$age_bracket)
    ]
  )
)

if (length(unassigned_ages) > 0) {
  
  stop(
    "The following ages were not assigned to a bracket: ",
    paste(
      unassigned_ages,
      collapse = ", "
    )
  )
}

############################################################
# 7. Keep only EN, IN, AST, and MG
############################################################

cells_keep <- colnames(reference)[
  !is.na(reference$age_bracket) &
    as.character(reference$broad_celltype) %in%
    c(
      "EN",
      "IN",
      "AST",
      "MG"
    )
]

reference_signature <- subset(
  reference,
  cells = cells_keep
)

############################################################
# 8. Create bracket × cell-type signature groups
############################################################

reference_signature$signature_group <- paste(
  as.character(
    reference_signature$age_bracket
  ),
  as.character(
    reference_signature$broad_celltype
  ),
  sep = "_"
)

desired_groups <- as.vector(
  unlist(
    lapply(
      age_bracket_order,
      function(age_bracket) {
        
        paste(
          age_bracket,
          c(
            "EN",
            "IN",
            "AST",
            "MG"
          ),
          sep = "_"
        )
      }
    )
  )
)

reference_signature$signature_group <- factor(
  reference_signature$signature_group,
  levels = desired_groups
)

############################################################
# 9. Check cell numbers per broad bracket
#
# Explicit group_by/summarise avoids the masked count()
# problem encountered earlier.
############################################################

cell_number_metadata <- data.frame(
  reference_age = as.character(
    reference_signature$reference_age
  ),
  age_bracket = as.character(
    reference_signature$age_bracket
  ),
  broad_celltype = as.character(
    reference_signature$broad_celltype
  ),
  signature_group = as.character(
    reference_signature$signature_group
  ),
  stringsAsFactors = FALSE
)

cell_number_table <- cell_number_metadata %>%
  dplyr::group_by(
    age_bracket,
    broad_celltype,
    signature_group
  ) %>%
  dplyr::summarise(
    n_nuclei = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    age_bracket = factor(
      age_bracket,
      levels = age_bracket_order,
      ordered = TRUE
    ),
    broad_celltype = factor(
      broad_celltype,
      levels = c(
        "EN",
        "IN",
        "AST",
        "MG"
      )
    )
  ) %>%
  dplyr::arrange(
    age_bracket,
    broad_celltype
  )

print(
  cell_number_table,
  n = Inf
)

############################################################
# 10. Also check exact ages contributing to each bracket
############################################################

exact_age_cell_number_table <- cell_number_metadata %>%
  dplyr::group_by(
    age_bracket,
    reference_age,
    broad_celltype
  ) %>%
  dplyr::summarise(
    n_nuclei = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    age_bracket = factor(
      age_bracket,
      levels = age_bracket_order,
      ordered = TRUE
    )
  ) %>%
  dplyr::arrange(
    age_bracket,
    reference_age,
    broad_celltype
  )

print(
  exact_age_cell_number_table,
  n = Inf
)

############################################################
# 11. Confirm all expected groups exist
############################################################

observed_groups <- unique(
  as.character(
    reference_signature$signature_group
  )
)

observed_groups <- observed_groups[
  !is.na(observed_groups)
]

missing_groups <- setdiff(
  desired_groups,
  observed_groups
)

if (length(missing_groups) > 0) {
  
  warning(
    "The following bracket/cell-type groups have no nuclei: ",
    paste(
      missing_groups,
      collapse = ", "
    )
  )
}

############################################################
# 12. Calculate average normalized expression
#
# AverageExpression with layer = "data":
#   1. exponentiates log-normalized values back to linear space
#   2. averages linear normalized expression across nuclei
#
# This matches the xeno_expression workflow.
############################################################

average_expression <- tryCatch(
  
  Seurat::AverageExpression(
    object = reference_signature,
    assays = "RNA",
    layer = "data",
    group.by = "signature_group",
    return.seurat = FALSE,
    verbose = FALSE
  ),
  
  error = function(e) {
    
    message(
      "Retrying AverageExpression using slot = 'data'."
    )
    
    Seurat::AverageExpression(
      object = reference_signature,
      assays = "RNA",
      slot = "data",
      group.by = "signature_group",
      return.seurat = FALSE,
      verbose = FALSE
    )
  }
)

average_matrix <- average_expression$RNA

############################################################
# 13. Standardize returned column names
############################################################

colnames(average_matrix) <- gsub(
  "-",
  "_",
  colnames(average_matrix),
  fixed = TRUE
)

cat(
  "\nAverageExpression returned these columns:\n"
)

print(
  colnames(average_matrix)
)

############################################################
# 14. Keep groups in chronological order
############################################################

available_desired_groups <- desired_groups[
  desired_groups %in%
    colnames(average_matrix)
]

missing_expression_groups <- setdiff(
  desired_groups,
  colnames(average_matrix)
)

if (length(missing_expression_groups) > 0) {
  
  warning(
    "These groups were absent from AverageExpression output: ",
    paste(
      missing_expression_groups,
      collapse = ", "
    )
  )
}

average_matrix <- average_matrix[
  ,
  available_desired_groups,
  drop = FALSE
]

############################################################
# 15. Create reference expression table
############################################################

reference_expression <- as.data.frame(
  as.matrix(average_matrix)
) %>%
  tibble::rownames_to_column(
    "gene"
  )

############################################################
# 16. Remove invalid and duplicated gene symbols
############################################################

reference_expression <- reference_expression %>%
  dplyr::filter(
    !is.na(gene),
    gene != ""
  ) %>%
  dplyr::group_by(
    gene
  ) %>%
  dplyr::summarise(
    dplyr::across(
      where(is.numeric),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

############################################################
# 17. Inspect final result
############################################################

cat(
  "\nFinal reference expression dimensions:\n"
)

print(
  dim(reference_expression)
)

cat(
  "\nFinal reference expression columns:\n"
)

print(
  colnames(reference_expression)
)

print(
  head(reference_expression)
)

############################################################
# 18. Save output
############################################################

output_dir <- file.path(
  getwd(),
  "GSE168408_reference_expression"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

readr::write_csv(
  reference_expression,
  file.path(
    output_dir,
    "GSE168408_broad_age_celltype_average_expression.csv"
  )
)

readr::write_csv(
  cell_number_table,
  file.path(
    output_dir,
    "GSE168408_broad_age_celltype_nuclei_counts.csv"
  )
)

readr::write_csv(
  exact_age_cell_number_table,
  file.path(
    output_dir,
    "GSE168408_exact_age_celltype_nuclei_counts.csv"
  )
)

saveRDS(
  reference_expression,
  file.path(
    output_dir,
    "GSE168408_broad_age_celltype_average_expression.rds"
  )
)

cat(
  "\nFiles saved to:\n",
  normalizePath(output_dir),
  "\n"
)








############################################################
# HuMiNAX 3R and 4R broad-cell-type expression signatures
#
# Input:
#   xeno_wt
#
# Output:
#   xeno_expression
#
# Columns:
#   gene
#   3R_EN, 3R_IN, 3R_AST, 3R_MG
#   4R_EN, 4R_IN, 4R_AST, 4R_MG
############################################################

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(readr)

############################################################
# 1. Use the RNA assay
############################################################

DefaultAssay(xeno_wt) <- "RNA"

# Seurat v5 objects can contain separate expression layers.
# Join them if multiple normalized data layers are present.
if (inherits(xeno_wt[["RNA"]], "Assay5")) {
  
  data_layers <- SeuratObject::Layers(
    xeno_wt[["RNA"]],
    search = "^data"
  )
  
  if (length(data_layers) > 1) {
    
    message(
      "Multiple RNA data layers detected. Joining layers."
    )
    
    xeno_wt <- SeuratObject::JoinLayers(
      xeno_wt,
      assay = "RNA"
    )
  }
  
  data_layers <- SeuratObject::Layers(
    xeno_wt[["RNA"]],
    search = "^data$"
  )
  
  # Normalize only if the RNA data layer is absent
  if (length(data_layers) == 0) {
    
    message(
      "No normalized RNA data layer found. Running NormalizeData."
    )
    
    xeno_wt <- NormalizeData(
      xeno_wt,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }
}

############################################################
# 2. Collapse the detailed cell types
############################################################

xeno_wt$broad_celltype <- case_when(
  
  grepl(
    "^Neu_EN_",
    as.character(xeno_wt$control_celltype)
  ) ~ "EN",
  
  grepl(
    "^Neu_IN_",
    as.character(xeno_wt$control_celltype)
  ) ~ "IN",
  
  as.character(xeno_wt$control_celltype) == "AST" ~ "AST",
  
  as.character(xeno_wt$control_celltype) == "MG" ~ "MG",
  
  TRUE ~ NA_character_
)

############################################################
# 3. Keep only the two control conditions and four cell types
############################################################

cells_keep <- colnames(xeno_wt)[
  as.character(xeno_wt$Condition) %in%
    c("3R_WT_Ctrl", "4R_WT_Ctrl") &
    !is.na(xeno_wt$broad_celltype)
]

xeno_signature <- subset(
  xeno_wt,
  cells = cells_keep
)

############################################################
# 4. Create the eight signature groups
############################################################

xeno_signature$isoform_group <- case_when(
  
  as.character(xeno_signature$Condition) ==
    "3R_WT_Ctrl" ~ "R3",
  
  as.character(xeno_signature$Condition) ==
    "4R_WT_Ctrl" ~ "R4",
  
  TRUE ~ NA_character_
)

xeno_signature$signature_group <- paste(
  xeno_signature$isoform_group,
  xeno_signature$broad_celltype,
  sep = "_"
)

desired_groups <- c(
  "R3_EN",
  "R3_IN",
  "R3_AST",
  "R3_MG",
  "R4_EN",
  "R4_IN",
  "R4_AST",
  "R4_MG"
)

xeno_signature$signature_group <- factor(
  xeno_signature$signature_group,
  levels = desired_groups
)

############################################################
# 5. Check the number of nuclei in each signature
############################################################

cell_number_table <- xeno_signature@meta.data %>%
  count(
    Condition,
    broad_celltype,
    signature_group,
    name = "n_nuclei"
  ) %>%
  arrange(
    signature_group
  )

print(cell_number_table)

missing_groups <- setdiff(
  desired_groups,
  unique(
    as.character(
      xeno_signature$signature_group
    )
  )
)

if (length(missing_groups) > 0) {
  
  stop(
    "The following groups have no nuclei: ",
    paste(missing_groups, collapse = ", ")
  )
}

############################################################
# 6. Calculate average normalized expression
#
# AverageExpression with layer = "data" converts the
# log-normalized values back to linear normalized expression
# before averaging across cells.
############################################################

average_expression <- tryCatch(
  
  Seurat::AverageExpression(
    object = xeno_signature,
    assays = "RNA",
    layer = "data",
    group.by = "signature_group",
    return.seurat = FALSE,
    verbose = FALSE
  ),
  
  error = function(e) {
    
    message(
      "Retrying AverageExpression using slot = 'data'."
    )
    
    Seurat::AverageExpression(
      object = xeno_signature,
      assays = "RNA",
      slot = "data",
      group.by = "signature_group",
      return.seurat = FALSE,
      verbose = FALSE
    )
  }
)

average_matrix <- average_expression$RNA

# Ensure the columns are ordered consistently
average_matrix <- average_matrix[
  ,
  desired_groups,
  drop = FALSE
]



colnames(average_matrix)

colnames(average_matrix) <- gsub(
  "-",
  "_",
  colnames(average_matrix),
  fixed = TRUE
)

print(colnames(average_matrix))

desired_groups <- c(
  "R3_EN",
  "R3_IN",
  "R3_AST",
  "R3_MG",
  "R4_EN",
  "R4_IN",
  "R4_AST",
  "R4_MG"
)

average_matrix <- average_matrix[
  ,
  desired_groups,
  drop = FALSE
]

xeno_expression <- as.data.frame(
  as.matrix(average_matrix)
) %>%
  tibble::rownames_to_column("gene")

colnames(xeno_expression) <- c(
  "gene",
  "3R_EN",
  "3R_IN",
  "3R_AST",
  "3R_MG",
  "4R_EN",
  "4R_IN",
  "4R_AST",
  "4R_MG"
)

head(xeno_expression)
dim(xeno_expression)

############################################################
# 7. Create the final expression table
############################################################

xeno_expression <- as.data.frame(
  as.matrix(average_matrix)
) %>%
  rownames_to_column("gene")

colnames(xeno_expression) <- c(
  "gene",
  "3R_EN",
  "3R_IN",
  "3R_AST",
  "3R_MG",
  "4R_EN",
  "4R_IN",
  "4R_AST",
  "4R_MG"
)


############################################################
# 8. Save 
############################################################

outdir <- "HuMiNAX_3R_4R_expression_signatures"

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

write_csv(
  xeno_expression,
  file.path(
    outdir,
    "HuMiNAX_3R_4R_EN_IN_AST_MG_expression.csv"
  )
)









############################################################
# GSE190815 hCS and transplanted hCS broad-cell signatures
#
# Output columns may include:
#   hCS_EN, hCS_IN, hCS_AST
#   t-hCS_EN, t-hCS_IN, t-hCS_AST
#
# Missing cell types are skipped automatically.
############################################################

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tibble)
library(purrr)
library(readr)

############################################################
# 1. File locations
############################################################

data_dir <- paste0(
  "C:/Users/wenhq/Documents/",
  "compare cell mature/hcs thcs"
)



hcs_file <- file.path(
  data_dir,
  "GSE190815_hCS_processed_SeuratObject.rds"
)

thcs_file <- file.path(
  data_dir,
  "GSE190815_t-hCS_processed_SeuratObject.rds"
)

stopifnot(
  file.exists(hcs_file),
  file.exists(thcs_file)
)



outdir <- file.path(
  data_dir,
  "GSE190815_average_signatures"
)

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

stopifnot(
  file.exists(hcs_file),
  file.exists(thcs_file)
)

############################################################
# 2. Load a gzipped RDS object
############################################################
read_rds_auto <- function(file_path) {
  
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  cat(
    "\nLoading:\n",
    file_path,
    "\nFile size:",
    round(file.info(file_path)$size / 1024^2, 2),
    "MB\n"
  )
  
  # Read the first few bytes without altering the file
  con <- file(file_path, open = "rb")
  header <- readBin(con, what = "raw", n = 8)
  close(con)
  
  header_hex <- paste(
    sprintf("%02x", as.integer(header)),
    collapse = " "
  )
  
  cat("File header:", header_hex, "\n")
  
  # Gzip magic bytes: 1f 8b
  is_gzip <- length(header) >= 2 &&
    as.integer(header[1]) == 0x1f &&
    as.integer(header[2]) == 0x8b
  
  if (is_gzip) {
    
    message("Detected a gzip-compressed RDS file.")
    
    con <- gzfile(file_path, open = "rb")
    on.exit(close(con), add = TRUE)
    
    object <- readRDS(con)
    
  } else {
    
    message(
      "File is not gzip-compressed. ",
      "Trying to read it directly as an RDS file."
    )
    
    object <- readRDS(file_path)
  }
  
  object
}




############################################################
# 3. Convert published labels to broad cell types
############################################################

map_broad_celltype <- function(labels) {
  
  labels_original <- as.character(labels)
  
  labels_clean <- tolower(
    trimws(labels_original)
  )
  
  broad <- rep(
    NA_character_,
    length(labels_clean)
  )
  
  # Astrocytes
  broad[
    grepl(
      "astro|astroglia",
      labels_clean
    )
  ] <- "AST"
  
  # Microglia, if present
  broad[
    grepl(
      "microglia|microglial|^mg$",
      labels_clean
    )
  ] <- "MG"
  
  # Inhibitory neurons
  broad[
    grepl(
      paste0(
        "gaba|inhibitory|interneuron|",
        "gaban|gabaergic"
      ),
      labels_clean
    )
  ] <- "IN"
  
  # Excitatory neurons
  broad[
    grepl(
      paste0(
        "glun|glutamatergic|excitatory|",
        "projection neuron"
      ),
      labels_clean
    )
  ] <- "EN"
  
  broad
}

############################################################
# 4. Identify the most likely annotation metadata column
############################################################

find_annotation_column <- function(object) {
  
  metadata <- object@meta.data
  
  candidate_results <- lapply(
    colnames(metadata),
    function(column_name) {
      
      values <- metadata[[column_name]]
      
      # Skip complex/list metadata columns
      if (
        !is.character(values) &&
        !is.factor(values)
      ) {
        return(NULL)
      }
      
      broad <- map_broad_celltype(values)
      
      data.frame(
        column = column_name,
        n_mapped = sum(
          !is.na(broad)
        ),
        n_celltypes = length(
          unique(
            broad[!is.na(broad)]
          )
        ),
        n_unique_labels = length(
          unique(
            as.character(values)
          )
        ),
        stringsAsFactors = FALSE
      )
    }
  )
  
  candidate_results <- bind_rows(
    candidate_results
  )
  
  candidate_results <- candidate_results %>%
    filter(
      n_mapped > 0
    ) %>%
    arrange(
      desc(n_celltypes),
      desc(n_mapped),
      n_unique_labels
    )
  
  if (nrow(candidate_results) == 0) {
    
    stop(
      paste0(
        "No metadata column with recognizable ",
        "EN, IN, AST, or MG labels was found.\n",
        "Run colnames(object@meta.data) and inspect ",
        "the annotations."
      )
    )
  }
  
  cat(
    "\nCandidate annotation columns:\n"
  )
  
  print(candidate_results)
  
  selected_column <- candidate_results$column[1]
  
  message(
    "\nAutomatically selected annotation column: ",
    selected_column
  )
  
  selected_column
}

############################################################
# 5. Ensure RNA log-normalized expression is available
############################################################

prepare_rna_assay <- function(object) {
  
  if (!"RNA" %in% Assays(object)) {
    
    stop(
      "The object does not contain an RNA assay. ",
      "Available assays: ",
      paste(
        Assays(object),
        collapse = ", "
      )
    )
  }
  
  DefaultAssay(object) <- "RNA"
  
  # Join layers only for Seurat v5 Assay5 objects
  if (
    inherits(
      object[["RNA"]],
      "Assay5"
    )
  ) {
    
    count_layers <- Layers(
      object[["RNA"]],
      search = "^counts"
    )
    
    data_layers <- Layers(
      object[["RNA"]],
      search = "^data"
    )
    
    if (
      length(count_layers) > 1 ||
      length(data_layers) > 1
    ) {
      
      message(
        "Joining multiple RNA layers."
      )
      
      object <- JoinLayers(
        object,
        assay = "RNA"
      )
    }
  }
  
  # Determine whether the normalized data slot/layer exists
  normalized_available <- tryCatch(
    {
      
      x <- GetAssayData(
        object,
        assay = "RNA",
        layer = "data"
      )
      
      nrow(x) > 0 &&
        ncol(x) > 0
      
    },
    error = function(e) {
      
      tryCatch(
        {
          
          x <- GetAssayData(
            object,
            assay = "RNA",
            slot = "data"
          )
          
          nrow(x) > 0 &&
            ncol(x) > 0
          
        },
        error = function(e2) {
          FALSE
        }
      )
    }
  )
  
  if (!normalized_available) {
    
    message(
      "RNA normalized data are absent. ",
      "Running LogNormalize on RNA counts."
    )
    
    object <- NormalizeData(
      object,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }
  
  object
}

############################################################
# 6. Process one object
############################################################

make_organoid_expression <- function(
    file_path,
    model_name,
    annotation_column = NULL
) {
  
  message(
    "\nLoading:\n",
    file_path
  )
  
  object <- readRDS(
    file_path
  )
  
  
  object <- prepare_rna_assay(
    object
  )
  
  cat(
    "\n============================\n",
    model_name,
    " metadata columns\n",
    "============================\n",
    sep = ""
  )
  
  print(
    colnames(object@meta.data)
  )
  
  annotation_column <- "cluster_label"
  
  original_labels <- as.character(
    object$cluster_label
  )
  
  object$broad_celltype <- dplyr::case_when(
    
    grepl(
      "^GluN_",
      original_labels
    ) ~ "EN",
    
    grepl(
      "^RELN_",
      original_labels
    ) ~ "EN",
    
    grepl(
      "^IN_",
      original_labels
    ) ~ "IN",
    
    grepl(
      "^Astroglia_",
      original_labels
    ) ~ "AST",
    
    TRUE ~ NA_character_
  )
  
  cat(
    "\n============================\n",
    model_name,
    " original annotations\n",
    "============================\n",
    sep = ""
  )
  
  annotation_summary <- data.frame(
    original_annotation = original_labels,
    broad_celltype = object$broad_celltype
  ) %>%
    count(
      original_annotation,
      broad_celltype,
      name = "n_nuclei"
    ) %>%
    arrange(
      broad_celltype,
      desc(n_nuclei)
    )
  
  print(annotation_summary)
  
  keep_cells <- colnames(object)[
    !is.na(
      object$broad_celltype
    )
  ]
  
  object_subset <- subset(
    object,
    cells = keep_cells
  )
  
  cell_counts <- object_subset@meta.data %>%
    count(
      broad_celltype,
      name = "n_nuclei"
    ) %>%
    arrange(
      broad_celltype
    )
  
  cat(
    "\n============================\n",
    model_name,
    " retained broad cell types\n",
    "============================\n",
    sep = ""
  )
  
  print(cell_counts)
  
  available_celltypes <- intersect(
    c(
      "EN",
      "IN",
      "AST",
      "MG"
    ),
    unique(
      as.character(
        object_subset$broad_celltype
      )
    )
  )
  
  if (length(available_celltypes) == 0) {
    
    stop(
      "No EN, IN, AST, or MG cells were identified in ",
      model_name,
      "."
    )
  }
  
  ##########################################################
  # Average all nuclei within each available cell type
  ##########################################################
  
  average_result <- tryCatch(
    {
      
      AverageExpression(
        object = object_subset,
        assays = "RNA",
        layer = "data",
        group.by = "broad_celltype",
        return.seurat = FALSE,
        verbose = FALSE
      )
      
    },
    error = function(e) {
      
      AverageExpression(
        object = object_subset,
        assays = "RNA",
        slot = "data",
        group.by = "broad_celltype",
        return.seurat = FALSE,
        verbose = FALSE
      )
    }
  )
  
  average_matrix <- average_result$RNA
  
  # Standardize any Seurat-generated column formatting
  clean_column_names <- colnames(
    average_matrix
  )
  
  clean_column_names <- sub(
    "^g",
    "",
    clean_column_names
  )
  
  clean_column_names <- gsub(
    "\\.",
    "_",
    clean_column_names
  )
  
  clean_column_names <- gsub(
    "-",
    "_",
    clean_column_names,
    fixed = TRUE
  )
  
  colnames(average_matrix) <-
    clean_column_names
  
  available_celltypes <- intersect(
    c(
      "EN",
      "IN",
      "AST",
      "MG"
    ),
    colnames(average_matrix)
  )
  
  average_matrix <- average_matrix[
    ,
    available_celltypes,
    drop = FALSE
  ]
  
  colnames(average_matrix) <- paste0(
    model_name,
    "_",
    colnames(average_matrix)
  )
  
  expression_table <- as.data.frame(
    as.matrix(average_matrix)
  ) %>%
    rownames_to_column(
      "gene"
    ) %>%
    filter(
      !is.na(gene),
      gene != ""
    ) %>%
    group_by(gene) %>%
    summarise(
      across(
        where(is.numeric),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  output <- list(
    expression = expression_table,
    cell_counts = cell_counts,
    annotation_summary = annotation_summary,
    annotation_column = annotation_column
  )
  
  rm(
    object,
    object_subset,
    average_result,
    average_matrix
  )
  
  gc()
  
  output
}

# unique(GSE190815_hCS_processed_SeuratObject$cluster_label)

############################################################
# 7. Process non-transplanted hCS
############################################################
hcs_result <- make_organoid_expression(
  file_path = hcs_file,
  model_name = "hCS"
)

############################################################
# 8. Process transplanted hCS
############################################################

thcs_result <- make_organoid_expression(
  file_path = thcs_file,
  model_name = "t-hCS"
)


############################################################
# 9. Merge hCS and t-hCS expression tables
############################################################

GSE190815_expression <- full_join(
  hcs_result$expression,
  thcs_result$expression,
  by = "gene"
)

############################################################
# 10. Put columns in a consistent order
############################################################

desired_columns <- c(
  "gene",
  "hCS_EN",
  "hCS_IN",
  "hCS_AST",
  "hCS_MG",
  "t-hCS_EN",
  "t-hCS_IN",
  "t-hCS_AST",
  "t-hCS_MG"
)

desired_columns <- intersect(
  desired_columns,
  colnames(GSE190815_expression)
)

GSE190815_expression <- GSE190815_expression %>%
  select(
    all_of(desired_columns)
  )

############################################################
# 11. Generate top 2,000 genes for every available signature
############################################################

signature_columns <- setdiff(
  colnames(GSE190815_expression),
  "gene"
)

GSE190815_top2000_long <- map_dfr(
  signature_columns,
  function(current_signature) {
    
    GSE190815_expression %>%
      select(
        gene,
        expression = all_of(
          current_signature
        )
      ) %>%
      filter(
        !is.na(expression),
        expression > 0
      ) %>%
      arrange(
        desc(expression)
      ) %>%
      slice_head(
        n = 2000
      ) %>%
      mutate(
        signature = current_signature,
        rank = row_number()
      ) %>%
      select(
        signature,
        rank,
        gene,
        expression
      )
  }
)

############################################################
# 12. Save final results
############################################################

write_csv(
  GSE190815_expression,
  file.path(
    outdir,
    "GSE190815_hCS_t-hCS_average_expression.csv"
  )
)









#Gage Cell microglia orgnaoid transplant paper
colnames(STS_integrated_CCA)


head(STS_integrated_CCA@meta.data )



sdrf_file <- paste0(
  "C:/Users/wenhq/Documents/",
  "compare cell mature/STS organoid/",
  "E-MTAB-11522.sdrf.txt"
)

sdrf <- read.delim(
  sdrf_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

dim(sdrf)
colnames(sdrf)
head(sdrf)


library(dplyr)
library(stringr)
library(tibble)

############################################################
# Build a clean table using only the needed SDRF columns
# This avoids the duplicated Protocol REF column names
############################################################

sdrf_selected <- tibble(
  source_name = sdrf[["Source Name"]],
  assay_name = sdrf[["Assay Name"]],
  submitted_file =
    sdrf[["Comment[SUBMITTED_FILE_NAME]"]],
  sampling_time_point = as.numeric(
    sdrf[["Factor Value[sampling time point]"]]
  ),
  sampling_time_unit =
    sdrf[["Unit[time unit]"]]
)

############################################################
# Convert assay names into Seurat orig.ident names
############################################################

sts_sample_map <- sdrf_selected %>%
  mutate(
    sample_prefix = assay_name,
    
    # 12_S2_L001 -> 12_S2
    # STS-G1_S4_L001 -> STS-G1_S4
    sample_prefix = str_remove(
      sample_prefix,
      "_L[0-9]+$"
    ),
    
    # 12_S2 -> 12
    # STS-G1_S4 -> STS-G1
    sample_prefix = str_remove(
      sample_prefix,
      "_S[0-9]+$"
    ),
    
    # STS-G1 -> STS_G1
    sample_prefix = str_replace_all(
      sample_prefix,
      "-",
      "_"
    ),
    
    # 12 -> STS_12
    # STS_G1 remains STS_G1
    seurat_orig_ident = if_else(
      str_detect(sample_prefix, "^STS_"),
      sample_prefix,
      paste0("STS_", sample_prefix)
    )
  ) %>%
  select(
    seurat_orig_ident,
    source_name,
    sampling_time_point,
    sampling_time_unit
  ) %>%
  distinct() %>%
  arrange(
    sampling_time_point,
    seurat_orig_ident
  )

print(sts_sample_map, n = Inf)



timepoint_lookup <- setNames(
  sts_sample_map$sampling_time_point,
  sts_sample_map$seurat_orig_ident
)

unit_lookup <- setNames(
  sts_sample_map$sampling_time_unit,
  sts_sample_map$seurat_orig_ident
)

source_lookup <- setNames(
  sts_sample_map$source_name,
  sts_sample_map$seurat_orig_ident
)

sample_ids <- as.character(
  STS_integrated_CCA$orig.ident
)

STS_integrated_CCA$sampling_time_point <- unname(
  timepoint_lookup[sample_ids]
)

STS_integrated_CCA$sampling_time_unit <- unname(
  unit_lookup[sample_ids]
)

STS_integrated_CCA$source_name <- unname(
  source_lookup[sample_ids]
)

STS_integrated_CCA$sampling_time_label <- paste0(
  STS_integrated_CCA$sampling_time_point,
  "_week"
)




table(
  STS_integrated_CCA$orig.ident,
  STS_integrated_CCA$sampling_time_label,
  useNA = "ifany"
)









############################################################
# Rusty Gage / Schafer transplanted hMG
# Average expression at 6, 12, and 24 weeks
############################################################

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(readr)

############################################################
# 1. Confirm the time-point annotation
############################################################

table(
  STS_integrated_CCA$orig.ident,
  STS_integrated_CCA$sampling_time_point,
  useNA = "ifany"
)

table(
  STS_integrated_CCA$sampling_time_point,
  useNA = "ifany"
)

############################################################
# 2. Create simple expression-signature labels
############################################################

STS_integrated_CCA$signature_group <- case_when(
  
  STS_integrated_CCA$sampling_time_point == 6 ~
    "Gage_6w_MG",
  
  STS_integrated_CCA$sampling_time_point == 12 ~
    "Gage_12w_MG",
  
  STS_integrated_CCA$sampling_time_point == 24 ~
    "Gage_24w_MG",
  
  TRUE ~ NA_character_
)

table(
  STS_integrated_CCA$signature_group,
  useNA = "ifany"
)

############################################################
# 3. Keep only cells with an assigned time point
############################################################

cells_keep <- colnames(STS_integrated_CCA)[
  !is.na(STS_integrated_CCA$signature_group)
]

STS_MG <- subset(
  STS_integrated_CCA,
  cells = cells_keep
)

############################################################
# 4. Use the original RNA assay
############################################################

DefaultAssay(STS_MG) <- "RNA"

# Join Seurat v5 layers only if needed
if (inherits(STS_MG[["RNA"]], "Assay5")) {
  
  data_layers <- Layers(
    STS_MG[["RNA"]],
    search = "^data"
  )
  
  count_layers <- Layers(
    STS_MG[["RNA"]],
    search = "^counts"
  )
  
  if (
    length(data_layers) > 1 ||
    length(count_layers) > 1
  ) {
    
    STS_MG <- JoinLayers(
      STS_MG,
      assay = "RNA"
    )
  }
}

############################################################
# 5. Normalize RNA expression
#
# Rerunning LogNormalize ensures that the three signatures
# are generated consistently from the RNA counts.
############################################################

STS_MG <- NormalizeData(
  STS_MG,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

############################################################
# 6. Average all cells within each time point
############################################################

gage_average <- AverageExpression(
  object = STS_MG,
  assays = "RNA",
  layer = "data",
  group.by = "signature_group",
  return.seurat = FALSE,
  verbose = FALSE
)$RNA

############################################################
# 7. Clean Seurat-generated column names
############################################################

colnames(gage_average) <- gsub(
  "-",
  "_",
  colnames(gage_average),
  fixed = TRUE
)

print(colnames(gage_average))

desired_columns <- c(
  "Gage_6w_MG",
  "Gage_12w_MG",
  "Gage_24w_MG"
)

missing_columns <- setdiff(
  desired_columns,
  colnames(gage_average)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing average-expression groups: ",
    paste(missing_columns, collapse = ", "),
    "\nAvailable columns: ",
    paste(colnames(gage_average), collapse = ", ")
  )
}

gage_average <- gage_average[
  ,
  desired_columns,
  drop = FALSE
]

############################################################
# 8. Create the standard expression table
############################################################

Gage_MG_expression <- as.data.frame(
  as.matrix(gage_average)
) %>%
  rownames_to_column("gene") %>%
  filter(
    !is.na(gene),
    gene != ""
  ) %>%
  group_by(gene) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

head(Gage_MG_expression)
dim(Gage_MG_expression)
colnames(Gage_MG_expression)









############################################################
#Process Park et al co-culture microglia with organoid
############################################################


library(Seurat)
library(Matrix)
library(dplyr)
library(tibble)
library(readr)

data_dir <- paste0(
  "C:/Users/wenhq/Documents/",
  "compare cell mature/coculture"
)

files <- c(
  coculture1 = file.path(
    data_dir,
    "GSM7774436_coculture1_hi_rc.txt.gz"
  ),
  coculture2 = file.path(
    data_dir,
    "GSM7774437_coculture2_hi_rc.txt.gz"
  ),
  coculture3 = file.path(
    data_dir,
    "GSM7774438_coculture3_hi_rc.txt.gz"
  )
)

stopifnot(all(file.exists(files)))








read_coculture_counts <- function(file, sample_name) {
  
  cat("\nReading:", basename(file), "\n")
  
  x <- read.delim(
    gzfile(file),
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # First column should contain gene names
  genes <- as.character(x[[1]])
  
  x <- x[, -1, drop = FALSE]
  
  counts <- Matrix(
    as.matrix(x),
    sparse = TRUE
  )
  
  rownames(counts) <- make.unique(genes)
  
  # Ensure cell barcodes are unique across samples
  colnames(counts) <- paste0(
    sample_name,
    "_",
    colnames(counts)
  )
  
  cat(
    "Genes:", nrow(counts),
    " Cells:", ncol(counts),
    "\n"
  )
  
  counts
}

count_list <- Map(
  read_coculture_counts,
  file = files,
  sample_name = names(files)
)







coculture_objects <- lapply(
  names(count_list),
  function(sample_name) {
    
    object <- CreateSeuratObject(
      counts = count_list[[sample_name]],
      project = sample_name,
      min.cells = 3,
      min.features = 200
    )
    
    object$sample <- sample_name
    
    object
  }
)

names(coculture_objects) <- names(count_list)

coculture <- merge(
  coculture_objects[[1]],
  y = coculture_objects[2:3],
  project = "Park_coculture"
)









coculture_objects <- lapply(
  names(count_list),
  function(sample_name) {
    
    object <- CreateSeuratObject(
      counts = count_list[[sample_name]],
      project = sample_name,
      min.cells = 3,
      min.features = 200
    )
    
    object$sample <- sample_name
    
    object
  }
)

names(coculture_objects) <- names(count_list)

coculture <- merge(
  coculture_objects[[1]],
  y = coculture_objects[2:3],
  project = "Park_coculture"
)




coculture[["percent.mt"]] <- PercentageFeatureSet(
  coculture,
  pattern = "^MT-"
)

VlnPlot(
  coculture,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "sample",
  ncol = 3
)










coculture <- subset(
  coculture,
  subset =
    nFeature_RNA >= 200 &
    nFeature_RNA <= 7000 &
    percent.mt < 20
)










DefaultAssay(coculture) <- "RNA"

coculture <- NormalizeData(
  coculture,
  verbose = FALSE
)

coculture <- FindVariableFeatures(
  coculture,
  nfeatures = 2000,
  verbose = FALSE
)

coculture <- ScaleData(
  coculture,
  verbose = FALSE
)

coculture <- RunPCA(
  coculture,
  npcs = 30,
  verbose = FALSE
)

coculture <- FindNeighbors(
  coculture,
  dims = 1:20,
  verbose = FALSE
)

coculture <- FindClusters(
  coculture,
  resolution = 0.4,
  verbose = FALSE
)

coculture <- RunUMAP(
  coculture,
  dims = 1:20,
  verbose = FALSE
)

DimPlot(
  coculture,
  group.by = "seurat_clusters",
  label = TRUE
)

DimPlot(
  coculture,
  group.by = "sample"
)
















microglia_markers <- c(
  "PTPRC",
  "SPI1",
  "AIF1",
  "CSF1R",
  "C1QA",
  "C1QB",
  "C1QC",
  "CX3CR1",
  "P2RY12",
  "P2RY13",
  "GPR34",
  "SALL1",
  "TMEM119",
  "APOE",
  "PLIN2"
)

neural_markers <- c(
  "SOX2",
  "VIM",
  "DCX",
  "STMN2",
  "TUBB3",
  "MAP2",
  "SLC17A7",
  "GAD1",
  "GAD2"
)

DotPlot(
  coculture,
  features = c(
    microglia_markers,
    neural_markers
  ),
  group.by = "seurat_clusters"
) +
  RotatedAxis()







microglia_clusters <- c("9")





Idents(coculture) <- "seurat_clusters"

iMicro <- subset(
  coculture,
  idents = microglia_clusters
)





iMicro$signature_group <- "Park_iMicro"

Park_iMicro_average <- AverageExpression(
  object = iMicro,
  assays = "RNA",
  layer = "data",
  group.by = "signature_group",
  return.seurat = FALSE,
  verbose = FALSE
)$RNA

colnames(Park_iMicro_average) <- "Park_iMicro"


Park_iMicro_expression <- as.data.frame(
  as.matrix(Park_iMicro_average)
) %>%
  rownames_to_column("gene") %>%
  filter(
    !is.na(gene),
    gene != ""
  ) %>%
  group_by(gene) %>%
  summarise(
    Park_iMicro = mean(
      Park_iMicro,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

head(Park_iMicro_expression)




write_csv(
  Park_iMicro_expression,
  file.path(
    data_dir,
    "Park_GSE242894_iMicro_average_expression.csv"
  )
)






#from Li et al Nature 2023


colnames(CHOOSE_CTRL_annot_srt@meta.data)

unique(CHOOSE_CTRL_annot_srt$celltype_cl_coarse2)

unique(CHOOSE_CTRL_annot_srt$celltype_cl)





############################################################
# Aggregate WT/control organoid EN, IN, and AST expression
############################################################

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tibble)

############################################################
# 1. Assign broad cell types
############################################################

CHOOSE_CTRL_annot_srt$broad_celltype <- case_when(
  
  CHOOSE_CTRL_annot_srt$celltype_cl %in% c(
    "L23",
    "L4",
    "L56",
    "L6_CThPN"
  ) ~ "EN",
  
  CHOOSE_CTRL_annot_srt$celltype_cl %in% c(
    "CGE_IN",
    "CGE_LGE_IN",
    "LGE_IN"
  ) ~ "IN",
  
  CHOOSE_CTRL_annot_srt$celltype_cl == "Astrocytes" ~
    "AST",
  
  TRUE ~ NA_character_
)

table(
  CHOOSE_CTRL_annot_srt$celltype_cl,
  CHOOSE_CTRL_annot_srt$broad_celltype,
  useNA = "ifany"
)

############################################################
# 2. Keep only EN, IN, and AST
############################################################

cells_keep <- colnames(CHOOSE_CTRL_annot_srt)[
  !is.na(CHOOSE_CTRL_annot_srt$broad_celltype)
]

CHOOSE_CTRL_broad <- subset(
  CHOOSE_CTRL_annot_srt,
  cells = cells_keep
)

table(CHOOSE_CTRL_broad$broad_celltype)

############################################################
# 3. Use RNA assay
############################################################

DefaultAssay(CHOOSE_CTRL_broad) <- "RNA"

if (inherits(CHOOSE_CTRL_broad[["RNA"]], "Assay5")) {
  
  CHOOSE_CTRL_broad <- JoinLayers(
    CHOOSE_CTRL_broad,
    assay = "RNA"
  )
}

############################################################
# 4. Normalize
############################################################

CHOOSE_CTRL_broad <- NormalizeData(
  CHOOSE_CTRL_broad,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

############################################################
# 5. Create WT/control organoid labels
############################################################

CHOOSE_CTRL_broad$signature_group <- case_when(
  
  CHOOSE_CTRL_broad$broad_celltype == "EN" ~
    "WT_organoid_EN",
  
  CHOOSE_CTRL_broad$broad_celltype == "IN" ~
    "WT_organoid_IN",
  
  CHOOSE_CTRL_broad$broad_celltype == "AST" ~
    "WT_organoid_AST"
)

table(CHOOSE_CTRL_broad$signature_group)

############################################################
# 6. Average all cells in each broad cell type
############################################################

CHOOSE_CTRL_average <- AverageExpression(
  object = CHOOSE_CTRL_broad,
  assays = "RNA",
  layer = "data",
  group.by = "signature_group",
  return.seurat = FALSE,
  verbose = FALSE
)$RNA

colnames(CHOOSE_CTRL_average) <- gsub(
  "-",
  "_",
  colnames(CHOOSE_CTRL_average),
  fixed = TRUE
)

desired_columns <- c(
  "WT_organoid_EN",
  "WT_organoid_IN",
  "WT_organoid_AST"
)

CHOOSE_CTRL_average <- CHOOSE_CTRL_average[
  ,
  desired_columns,
  drop = FALSE
]

############################################################
# 7. Final expression table
############################################################

CHOOSE_CTRL_expression <- as.data.frame(
  as.matrix(CHOOSE_CTRL_average)
) %>%
  rownames_to_column("gene") %>%
  filter(
    !is.na(gene),
    gene != ""
  ) %>%
  group_by(gene) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

head(CHOOSE_CTRL_expression)
colnames(CHOOSE_CTRL_expression)










############################################################
# GSE185472 in-vitro and transplanted organoid signatures
#
# Included:
#   10-week in-vitro organoids:
#     GSM6679388  HUES6
#     GSM6679389  iPSC822
#
#   5-month in-vitro organoids:
#     GSM6679390  HUES6
#     GSM6679391  iPSC822
#
#   Transplanted organoids:
#     GSM6679384  HUES6,   5 months
#     GSM6679385  iPSC822, 5 months
#     GSM5615952  iPSC822, 6 months
#     GSM5615953  HUES6,   6 months
#     GSM6679386  HUES6,   8 months
#     GSM6679387  iPSC822, 8 months
#
# Excluded completely:
#   GSM7286711  HUES6 8-month transplant CTRL scRNA-seq
#   GSM7286712  HUES6 8-month transplant TNFa scRNA-seq
#
# Cell-type definitions from the authors' annotation:
#   EN  = UL + DL
#   IN  = iN
#   AST = Ast
#
# Main output:
#   GSE185472_EN_IN_AST_expression
#
# CSV:
#   GSE185472_invitro_transplant_EN_IN_AST_expression.csv
#
# The two cell lines are pooled within each condition using
# nucleus-number weighting. Therefore, every retained nucleus
# contributes equally to the final condition/cell-type mean.
############################################################

library(Seurat)
library(SeuratObject)
library(Matrix)
library(dplyr)
library(tibble)
library(readr)
library(stringr)

############################################################
# 1. File locations
############################################################

base_directory <- "C:/Users/wenhq/Documents/GSE185472"

raw_directory <- file.path(
  base_directory,
  "GSE185472_RAW"
)

annotation_files <- c(
  ten_week = file.path(
    base_directory,
    "GSE185472_10wOrg_withAnnot.txt.gz"
  ),
  integrated = file.path(
    base_directory,
    "GSE185472_integratedorganoid_withAnnot.txt.gz"
  )
)

output_directory <- file.path(
  base_directory,
  "GSE185472_EN_IN_AST_processed"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  dir.exists(raw_directory),
  all(file.exists(annotation_files))
)

############################################################
# 2. Define the ten samples to include
#
# GSM7286711 and GSM7286712 are intentionally absent.
############################################################

sample_information <- tribble(
  ~sample_key,       ~GSM,         ~cell_line, ~condition,       ~annotation_source,
  "HUES6_10w_IV",    "GSM6679388", "HUES6",    "10w_invitro",    "ten_week",
  "iPSC822_10w_IV",  "GSM6679389", "iPSC822",  "10w_invitro",    "ten_week",
  "HUES6_5m_IV",     "GSM6679390", "HUES6",    "5m_invitro",     "integrated",
  "iPSC822_5m_IV",   "GSM6679391", "iPSC822",  "5m_invitro",     "integrated",
  "HUES6_5m_TX",     "GSM6679384", "HUES6",    "5m_transplant",  "integrated",
  "iPSC822_5m_TX",   "GSM6679385", "iPSC822",  "5m_transplant",  "integrated",
  "iPSC822_6m_TX",   "GSM5615952", "iPSC822",  "6m_transplant",  "integrated",
  "HUES6_6m_TX",     "GSM5615953", "HUES6",    "6m_transplant",  "integrated",
  "HUES6_8m_TX",     "GSM6679386", "HUES6",    "8m_transplant",  "integrated",
  "iPSC822_8m_TX",   "GSM6679387", "iPSC822",  "8m_transplant",  "integrated"
)

desired_sample_order <- sample_information$sample_key

desired_conditions <- c(
  "10w_invitro",
  "5m_invitro",
  "5m_transplant",
  "6m_transplant",
  "8m_transplant"
)

desired_celltypes <- c(
  "EN",
  "IN",
  "AST"
)

############################################################
# 3. Helper functions
############################################################

############################################################
# 3A. Read a compressed or uncompressed table
############################################################

read_table_no_header <- function(filename) {
  
  if (!file.exists(filename)) {
    stop("File does not exist: ", filename)
  }
  
  if (grepl("\\.gz$", filename, ignore.case = TRUE)) {
    
    connection <- gzfile(
      filename,
      open = "rt"
    )
    
    on.exit(
      close(connection),
      add = TRUE
    )
    
    output <- read.delim(
      connection,
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    )
    
  } else {
    
    output <- read.delim(
      filename,
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    )
  }
  
  output
}

############################################################
# 3B. Read barcodes
############################################################

read_barcode_file <- function(filename) {
  
  if (!file.exists(filename)) {
    stop("Barcode file does not exist: ", filename)
  }
  
  if (grepl("\\.gz$", filename, ignore.case = TRUE)) {
    
    connection <- gzfile(
      filename,
      open = "rt"
    )
    
    on.exit(
      close(connection),
      add = TRUE
    )
    
    barcodes <- readLines(connection)
    
  } else {
    
    barcodes <- readLines(filename)
  }
  
  barcodes
}

############################################################
# 3C. Read an MTX matrix
############################################################

read_matrix_market_file <- function(filename) {
  
  if (!file.exists(filename)) {
    stop("Matrix file does not exist: ", filename)
  }
  
  if (grepl("\\.gz$", filename, ignore.case = TRUE)) {
    
    connection <- gzfile(
      filename,
      open = "rb"
    )
    
    on.exit(
      close(connection),
      add = TRUE
    )
    
    counts <- Matrix::readMM(connection)
    
  } else {
    
    counts <- Matrix::readMM(filename)
  }
  
  methods::as(
    counts,
    "dgCMatrix"
  )
}

############################################################
# 3D. Extract the nucleotide barcode
############################################################

extract_barcode_core <- function(x) {
  
  stringr::str_extract(
    toupper(as.character(x)),
    "[ACGTN]{12,}"
  )
}

############################################################
# 3E. Require one matching file
############################################################

choose_single_file <- function(
    candidate_files,
    file_description
) {
  
  candidate_files <- unique(
    candidate_files[
      file.exists(candidate_files)
    ]
  )
  
  if (length(candidate_files) > 1) {
    
    filtered_candidates <- candidate_files[
      grepl(
        "filtered",
        candidate_files,
        ignore.case = TRUE
      )
    ]
    
    if (length(filtered_candidates) == 1) {
      candidate_files <- filtered_candidates
    }
  }
  
  if (length(candidate_files) == 0) {
    
    stop(
      "\nNo ",
      file_description,
      " was found."
    )
  }
  
  if (length(candidate_files) > 1) {
    
    stop(
      "\nMore than one ",
      file_description,
      " was found:\n",
      paste(
        candidate_files,
        collapse = "\n"
      )
    )
  }
  
  normalizePath(
    candidate_files,
    winslash = "/",
    mustWork = TRUE
  )
}

############################################################
# 3F. Locate matrix, barcode, and feature files
############################################################

find_matrix_files <- function(
    GSM,
    raw_directory
) {
  
  all_files <- list.files(
    path = raw_directory,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )
  
  tar_files <- all_files[
    grepl(
      GSM,
      basename(all_files),
      fixed = TRUE
    ) &
      grepl(
        "\\.tar\\.gz$",
        basename(all_files),
        ignore.case = TRUE
      )
  ]
  
  if (length(tar_files) > 1) {
    stop(
      "Multiple tar.gz files found for ",
      GSM
    )
  }
  
  if (length(tar_files) == 1) {
    
    extraction_directory <- file.path(
      raw_directory,
      "extracted",
      GSM
    )
    
    dir.create(
      extraction_directory,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    extracted_files <- list.files(
      extraction_directory,
      recursive = TRUE,
      full.names = TRUE,
      include.dirs = FALSE
    )
    
    if (length(extracted_files) == 0) {
      
      message(
        "Extracting ",
        basename(tar_files)
      )
      
      utils::untar(
        tarfile = tar_files,
        exdir = extraction_directory
      )
    }
    
    sample_files <- list.files(
      extraction_directory,
      recursive = TRUE,
      full.names = TRUE,
      include.dirs = FALSE
    )
    
  } else {
    
    sample_files <- all_files[
      grepl(
        GSM,
        basename(all_files),
        fixed = TRUE
      )
    ]
  }
  
  matrix_candidates <- sample_files[
    grepl(
      "matrix\\.mtx(\\.gz)?$",
      basename(sample_files),
      ignore.case = TRUE
    )
  ]
  
  barcode_candidates <- sample_files[
    grepl(
      "barcodes\\.tsv(\\.gz)?$",
      basename(sample_files),
      ignore.case = TRUE
    )
  ]
  
  feature_candidates <- sample_files[
    grepl(
      "(features|genes)\\.tsv(\\.gz)?$",
      basename(sample_files),
      ignore.case = TRUE
    )
  ]
  
  list(
    matrix = choose_single_file(
      matrix_candidates,
      paste0(GSM, " matrix file")
    ),
    barcodes = choose_single_file(
      barcode_candidates,
      paste0(GSM, " barcode file")
    ),
    features = choose_single_file(
      feature_candidates,
      paste0(GSM, " feature file")
    )
  )
}

############################################################
# 3G. Strip genome prefixes from gene names
############################################################

strip_genome_prefix <- function(x) {
  
  x <- sub(
    "^(hg19|GRCh37|GRCh38)[-_]",
    "",
    x,
    ignore.case = TRUE
  )
  
  x <- sub(
    "^(mm10|GRCm38|GRCm39)[-_]",
    "",
    x,
    ignore.case = TRUE
  )
  
  x
}

############################################################
# 3H. Load one sample and retain human gene counts
#
# This handles both:
#   - human-only in-vitro matrices
#   - combined hg19/mm10 transplant matrices
############################################################

load_human_count_matrix <- function(
    matrix_file,
    barcode_file,
    feature_file
) {
  
  message(
    "Reading matrix: ",
    basename(matrix_file)
  )
  
  counts <- read_matrix_market_file(
    matrix_file
  )
  
  barcodes <- read_barcode_file(
    barcode_file
  )
  
  features <- read_table_no_header(
    feature_file
  )
  
  if (ncol(features) < 2) {
    
    gene_ids <- as.character(
      features[[1]]
    )
    
    gene_symbols <- gene_ids
    
  } else {
    
    gene_ids <- as.character(
      features[[1]]
    )
    
    gene_symbols <- as.character(
      features[[2]]
    )
  }
  
  correct_orientation <- (
    nrow(counts) == nrow(features) &&
      ncol(counts) == length(barcodes)
  )
  
  reversed_orientation <- (
    nrow(counts) == length(barcodes) &&
      ncol(counts) == nrow(features)
  )
  
  if (reversed_orientation) {
    
    message(
      "Transposing the count matrix."
    )
    
    counts <- Matrix::t(
      counts
    )
    
    counts <- methods::as(
      counts,
      "dgCMatrix"
    )
  }
  
  if (!correct_orientation && !reversed_orientation) {
    
    stop(
      "\nMatrix dimensions do not match.",
      "\nMatrix dimensions: ",
      nrow(counts),
      " x ",
      ncol(counts),
      "\nFeature rows: ",
      nrow(features),
      "\nBarcodes: ",
      length(barcodes)
    )
  }
  
  ##########################################################
  # Keep Gene Expression features when a feature-type column
  # is present.
  ##########################################################
  
  feature_type_keep <- rep(
    TRUE,
    nrow(features)
  )
  
  if (ncol(features) >= 3) {
    
    feature_types <- as.character(
      features[[3]]
    )
    
    if (
      any(
        grepl(
          "Gene Expression",
          feature_types,
          ignore.case = TRUE
        )
      )
    ) {
      
      feature_type_keep <- grepl(
        "Gene Expression",
        feature_types,
        ignore.case = TRUE
      )
    }
  }
  
  ##########################################################
  # Identify species.
  ##########################################################
  
  human_prefix <- (
    grepl(
      "^(hg19|GRCh37|GRCh38)[-_]",
      gene_ids,
      ignore.case = TRUE
    ) |
      grepl(
        "^(hg19|GRCh37|GRCh38)[-_]",
        gene_symbols,
        ignore.case = TRUE
      )
  )
  
  mouse_prefix <- (
    grepl(
      "^(mm10|GRCm38|GRCm39)[-_]",
      gene_ids,
      ignore.case = TRUE
    ) |
      grepl(
        "^(mm10|GRCm38|GRCm39)[-_]",
        gene_symbols,
        ignore.case = TRUE
      )
  )
  
  plain_human_ensembl <- grepl(
    "^ENSG",
    gene_ids,
    ignore.case = TRUE
  )
  
  plain_mouse_ensembl <- grepl(
    "^ENSMUSG",
    gene_ids,
    ignore.case = TRUE
  )
  
  if (sum(human_prefix) >= 10000) {
    
    species_keep <- human_prefix
    
    message(
      "Combined-genome matrix detected; retaining ",
      sum(species_keep),
      " human features."
    )
    
  } else if (
    sum(plain_human_ensembl) >= 10000 &&
    sum(plain_mouse_ensembl) > 0
  ) {
    
    species_keep <- plain_human_ensembl
    
    message(
      "Mixed Ensembl matrix detected; retaining ",
      sum(species_keep),
      " human ENSG features."
    )
    
  } else if (sum(mouse_prefix) > 0) {
    
    species_keep <- !mouse_prefix
    
    message(
      "Removing ",
      sum(mouse_prefix),
      " explicitly mouse-prefixed features."
    )
    
  } else if (sum(plain_mouse_ensembl) > 0) {
    
    species_keep <- !plain_mouse_ensembl
    
    message(
      "Removing ",
      sum(plain_mouse_ensembl),
      " mouse Ensembl features."
    )
    
  } else {
    
    species_keep <- rep(
      TRUE,
      length(gene_ids)
    )
    
    message(
      "Human-only matrix assumed; all gene features retained."
    )
  }
  
  ##########################################################
  # Clean gene names and apply filters.
  ##########################################################
  
  gene_symbols_clean <- strip_genome_prefix(
    gene_symbols
  )
  
  gene_ids_clean <- strip_genome_prefix(
    gene_ids
  )
  
  use_gene_id <- (
    is.na(gene_symbols_clean) |
      trimws(gene_symbols_clean) == ""
  )
  
  gene_symbols_clean[
    use_gene_id
  ] <- gene_ids_clean[
    use_gene_id
  ]
  
  keep_features <- (
    feature_type_keep &
      species_keep &
      !is.na(gene_symbols_clean) &
      trimws(gene_symbols_clean) != ""
  )
  
  counts <- counts[
    keep_features,
    ,
    drop = FALSE
  ]
  
  gene_symbols_clean <- gene_symbols_clean[
    keep_features
  ]
  
  ##########################################################
  # Collapse duplicated symbols by summing raw counts.
  ##########################################################
  
  unique_gene_symbols <- unique(
    gene_symbols_clean
  )
  
  gene_mapping_matrix <- Matrix::sparseMatrix(
    i = match(
      gene_symbols_clean,
      unique_gene_symbols
    ),
    j = seq_along(
      gene_symbols_clean
    ),
    x = 1,
    dims = c(
      length(unique_gene_symbols),
      length(gene_symbols_clean)
    )
  )
  
  counts <- gene_mapping_matrix %*% counts
  
  counts <- methods::as(
    counts,
    "dgCMatrix"
  )
  
  rownames(counts) <- unique_gene_symbols
  colnames(counts) <- barcodes
  
  message(
    "Final human matrix: ",
    nrow(counts),
    " genes x ",
    ncol(counts),
    " nuclei"
  )
  
  counts
}

############################################################
# 3I. Read and standardize one annotation table
############################################################

############################################################
# Read and standardize a GSE185472 annotation table
############################################################

read_annotation_table <- function(
    filename,
    dataset_name,
    default_sample_id = NULL
) {
  
  ##########################################################
  # Read annotation file
  ##########################################################
  
  file_connection <- if (
    grepl(
      "\\.gz$",
      filename,
      ignore.case = TRUE
    )
  ) {
    gzfile(filename)
  } else {
    filename
  }
  
  annotation <- read.delim(
    file_connection,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  ##########################################################
  # Clean column names
  #
  # The first column in the 10-week file has a blank name.
  # Rename blank columns before trying to extract them.
  ##########################################################
  
  cleaned_column_names <- trimws(
    colnames(annotation)
  )
  
  blank_column_indices <- which(
    is.na(cleaned_column_names) |
      cleaned_column_names == ""
  )
  
  if (length(blank_column_indices) > 0) {
    
    cleaned_column_names[
      blank_column_indices
    ] <- paste0(
      "unnamed_column_",
      seq_along(blank_column_indices)
    )
  }
  
  colnames(annotation) <- make.unique(
    cleaned_column_names,
    sep = "_"
  )
  
  cat(
    "\nColumns detected in ",
    basename(filename),
    ":\n",
    sep = ""
  )
  
  print(
    colnames(annotation)
  )
  
  ##########################################################
  # Find one column using case-insensitive matching
  ##########################################################
  
  find_column <- function(candidates) {
    
    existing_names <- colnames(annotation)
    
    match_index <- match(
      tolower(candidates),
      tolower(existing_names)
    )
    
    match_index <- match_index[
      !is.na(match_index)
    ]
    
    if (length(match_index) == 0) {
      return(NA_character_)
    }
    
    existing_names[
      match_index[1]
    ]
  }
  
  ##########################################################
  # Identify barcode column
  ##########################################################
  
  barcode_column <- find_column(
    c(
      "barcode",
      "cell_barcode",
      "cellbarcode",
      "cell_id",
      "cellid",
      "cell",
      "unnamed_column_1",
      "x"
    )
  )
  
  if (is.na(barcode_column)) {
    
    warning(
      "Could not identify a named barcode column in ",
      basename(filename),
      ". Using the first column."
    )
    
    barcode_values <- annotation[[1]]
    
  } else {
    
    barcode_values <- annotation[[barcode_column]]
  }
  
  ##########################################################
  # Identify cluster column
  ##########################################################
  
  cluster_column <- find_column(
    c(
      "seurat_clusters",
      "seurat_cluster",
      "cluster",
      "clusters",
      "author_cluster"
    )
  )
  
  if (is.na(cluster_column)) {
    
    stop(
      "Could not identify the cluster column in ",
      basename(filename),
      ". Detected columns: ",
      paste(
        colnames(annotation),
        collapse = ", "
      )
    )
  }
  
  ##########################################################
  # Identify author cell-type column
  ##########################################################
  
  celltype_column <- find_column(
    c(
      "CellType",
      "cell_type",
      "celltype",
      "Cell.Type",
      "annotation",
      "author_celltype"
    )
  )
  
  if (is.na(celltype_column)) {
    
    stop(
      "Could not identify the cell-type column in ",
      basename(filename),
      ". Detected columns: ",
      paste(
        colnames(annotation),
        collapse = ", "
      )
    )
  }
  
  ##########################################################
  # Identify SampleID column
  ##########################################################
  
  sample_id_column <- find_column(
    c(
      "SampleID",
      "sample_id",
      "sampleID",
      "sample",
      "orig.ident",
      "orig_ident"
    )
  )
  
  ##########################################################
  # Create standardized columns
  ##########################################################
  
  annotation$barcode <- trimws(
    as.character(
      barcode_values
    )
  )
  
  annotation$seurat_clusters <- trimws(
    as.character(
      annotation[[cluster_column]]
    )
  )
  
  annotation$CellType <- trimws(
    as.character(
      annotation[[celltype_column]]
    )
  )
  
  if (!is.na(sample_id_column)) {
    
    annotation$SampleID <- trimws(
      as.character(
        annotation[[sample_id_column]]
      )
    )
    
  } else {
    
    if (is.null(default_sample_id)) {
      
      stop(
        "Could not identify the SampleID column in ",
        basename(filename),
        ", and default_sample_id was not provided."
      )
    }
    
    annotation$SampleID <- rep(
      default_sample_id,
      nrow(annotation)
    )
    
    message(
      "No SampleID column was present in ",
      basename(filename),
      ". Assigned SampleID = ",
      default_sample_id,
      "."
    )
  }
  
  ##########################################################
  # Remove rows without a usable barcode
  ##########################################################
  
  annotation <- annotation[
    !is.na(annotation$barcode) &
      annotation$barcode != "",
    ,
    drop = FALSE
  ]
  
  ##########################################################
  # Print annotation summary
  ##########################################################
  
  cat(
    "\nAnnotation summary for ",
    dataset_name,
    ":\n",
    sep = ""
  )
  
  cat(
    "Cells: ",
    nrow(annotation),
    "\n",
    sep = ""
  )
  
  cat(
    "\nSample IDs:\n"
  )
  
  print(
    table(
      annotation$SampleID,
      useNA = "ifany"
    )
  )
  
  cat(
    "\nAuthor cell types:\n"
  )
  
  print(
    table(
      annotation$CellType,
      useNA = "ifany"
    )
  )
  
  return(annotation)
}

############################################################
# 3J. Normalize SampleID strings for robust matching
############################################################

normalize_sample_id <- function(x) {
  
  x %>%
    as.character() %>%
    trimws() %>%
    tolower() %>%
    stringr::str_replace_all(
      "[^a-z0-9]",
      ""
    )
}

############################################################
# 4. Read both author annotation files
############################################################

# ############################################################
# # Read 10-week in vitro organoid annotation
# ############################################################
# 
# annotation_10w <- read_annotation_table(
#   filename = annotation_files[["ten_week"]],
#   dataset_name = "GSE185472 10-week organoid",
#   default_sample_id = "GSE185472_10w_organoid"
# )


############################################################
# Read integrated organoid/transplant annotation
############################################################

annotation_integrated <- read_annotation_table(
  filename = annotation_files[["integrated"]],
  dataset_name = "GSE185472 integrated organoid and transplant"
)

# 
# 
# 
# 
# 
# 
# 
# 
# annotation_10w <- read_annotation_table(
#   filename = annotation_files[["ten_week"]],
#   annotation_source = "ten_week"
# )
# 
# annotation_integrated <- read_annotation_table(
#   filename = annotation_files[["integrated"]],
#   annotation_source = "integrated"
# )

cat(
  "\n10-week SampleID values:\n"
)

print(
  sort(
    unique(
      as.character(
        annotation_10w$SampleID
      )
    )
  )
)

cat(
  "\nIntegrated SampleID values:\n"
)

print(
  sort(
    unique(
      as.character(
        annotation_integrated$SampleID
      )
    )
  )
)

annotation <- dplyr::bind_rows(
  annotation_10w,
  annotation_integrated
)


############################################################
# 5. Prepare annotation
#
# Included:
#   5-month in-vitro organoids
#   5-month transplanted organoids
#   6-month transplanted organoids
#   8-month transplanted organoids
#
# Excluded:
#   10-week in-vitro organoids
############################################################


############################################################
# 5A. Use only the integrated annotation
############################################################

required_annotation_columns <- c(
  "barcode",
  "SampleID",
  "seurat_clusters",
  "CellType"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  colnames(annotation_integrated)
)

if (length(missing_annotation_columns) > 0) {
  
  stop(
    "annotation_integrated is missing: ",
    paste(
      missing_annotation_columns,
      collapse = ", "
    )
  )
}

annotation <- annotation_integrated


############################################################
# 5B. Clean annotation metadata
############################################################

annotation <- annotation %>%
  dplyr::mutate(
    
    barcode = trimws(
      as.character(barcode)
    ),
    
    SampleID = trimws(
      as.character(SampleID)
    ),
    
    seurat_clusters = trimws(
      as.character(seurat_clusters)
    ),
    
    CellType = trimws(
      as.character(CellType)
    ),
    
    annotation_source = "integrated",
    
    sample_id_clean = tolower(
      trimws(
        as.character(SampleID)
      )
    )
  )


############################################################
# 5C. Identify cell line
############################################################

annotation <- annotation %>%
  dplyr::mutate(
    
    cell_line = dplyr::case_when(
      
      grepl(
        "^hues6",
        sample_id_clean
      ) ~ "HUES6",
      
      grepl(
        "^ip822",
        sample_id_clean
      ) ~ "iP822",
      
      TRUE ~ NA_character_
    )
  )


############################################################
# 5D. Identify experimental condition
############################################################

annotation <- annotation %>%
  dplyr::mutate(
    
    condition_detected = dplyr::case_when(
      
      grepl(
        "_5m_org$",
        sample_id_clean
      ) ~ "5m_in_vitro",
      
      grepl(
        "_5m_t$",
        sample_id_clean
      ) ~ "5m_transplant",
      
      grepl(
        "_6m_t$",
        sample_id_clean
      ) ~ "6m_transplant",
      
      grepl(
        "_8m_t$",
        sample_id_clean
      ) ~ "8m_transplant",
      
      TRUE ~ NA_character_
    )
  )


############################################################
# 5E. Create unique sample keys
############################################################

annotation <- annotation %>%
  dplyr::mutate(
    
    sample_key = dplyr::case_when(
      
      sample_id_clean == "hues6_5m_org" ~
        "HUES6_5m_in_vitro",
      
      sample_id_clean == "ip822_5m_org" ~
        "iP822_5m_in_vitro",
      
      sample_id_clean == "hues6_5m_t" ~
        "HUES6_5m_transplant",
      
      sample_id_clean == "ip822_5m_t" ~
        "iP822_5m_transplant",
      
      sample_id_clean == "hues6_6m_t" ~
        "HUES6_6m_transplant",
      
      sample_id_clean == "ip822_6m_t" ~
        "iP822_6m_transplant",
      
      sample_id_clean == "hues6_8m_t" ~
        "HUES6_8m_transplant",
      
      sample_id_clean == "ip822_8m_t" ~
        "iP822_8m_transplant",
      
      TRUE ~ NA_character_
    )
  )


############################################################
# 5F. Remove unusable annotation rows
############################################################

annotation <- annotation %>%
  dplyr::filter(
    !is.na(barcode),
    barcode != ""
  )


############################################################
# 5G. Check for unrecognized SampleID values
############################################################

unmatched_annotation <- annotation %>%
  dplyr::filter(
    is.na(cell_line) |
      is.na(condition_detected) |
      is.na(sample_key)
  ) %>%
  dplyr::distinct(
    SampleID,
    sample_id_clean,
    cell_line,
    condition_detected,
    sample_key
  )

if (nrow(unmatched_annotation) > 0) {
  
  cat(
    "\nUnrecognized SampleID values:\n"
  )
  
  print(
    unmatched_annotation
  )
  
  stop(
    "Some integrated annotation rows could not be assigned."
  )
}


############################################################
# 5H. Print annotation summary
############################################################

cat(
  "\n========================================\n",
  "Integrated annotation summary\n",
  "10-week organoids are excluded\n",
  "========================================\n"
)

print(
  annotation %>%
    dplyr::count(
      SampleID,
      cell_line,
      condition_detected,
      sample_key,
      name = "n_annotated_cells"
    ) %>%
    dplyr::arrange(
      condition_detected,
      cell_line
    )
)

cat(
  "\nAuthor cell types by condition:\n"
)

print(
  with(
    annotation,
    table(
      condition_detected,
      CellType,
      useNA = "ifany"
    )
  )
)


############################################################
# 6. Create sample metadata and attach GSM identifiers
#
# Eight samples total:
#   Two 5-month in-vitro organoids
#   Six transplanted organoids
############################################################


############################################################
# 6A. Define included samples
############################################################

sample_information <- tibble::tribble(
  
  ~sample_key,                ~GSM,         ~cell_line, ~condition_detected,
  
  "HUES6_5m_in_vitro",       "GSM6679390", "HUES6",    "5m_in_vitro",
  "iP822_5m_in_vitro",       "GSM6679391", "iP822",    "5m_in_vitro",
  
  "HUES6_5m_transplant",     "GSM6679384", "HUES6",    "5m_transplant",
  "iP822_5m_transplant",     "GSM6679385", "iP822",    "5m_transplant",
  
  "iP822_6m_transplant",     "GSM5615952", "iP822",    "6m_transplant",
  "HUES6_6m_transplant",     "GSM5615953", "HUES6",    "6m_transplant",
  
  "HUES6_8m_transplant",     "GSM6679386", "HUES6",    "8m_transplant",
  "iP822_8m_transplant",     "GSM6679387", "iP822",    "8m_transplant"
)


############################################################
# 6B. Validate sample table
############################################################

if (anyDuplicated(sample_information$sample_key) > 0) {
  stop("sample_information contains duplicated sample_key values.")
}

if (anyDuplicated(sample_information$GSM) > 0) {
  stop("sample_information contains duplicated GSM identifiers.")
}


############################################################
# 6C. Confirm annotation and sample metadata agree
############################################################

annotation_sample_keys <- sort(
  unique(
    annotation$sample_key
  )
)

metadata_sample_keys <- sort(
  sample_information$sample_key
)

missing_from_metadata <- setdiff(
  annotation_sample_keys,
  metadata_sample_keys
)

if (length(missing_from_metadata) > 0) {
  
  stop(
    "Annotation sample keys missing from sample_information: ",
    paste(
      missing_from_metadata,
      collapse = ", "
    )
  )
}

missing_from_annotation <- setdiff(
  metadata_sample_keys,
  annotation_sample_keys
)

if (length(missing_from_annotation) > 0) {
  
  stop(
    "sample_information entries missing from annotation: ",
    paste(
      missing_from_annotation,
      collapse = ", "
    )
  )
}


############################################################
# 6D. Confirm cell line and condition assignments
############################################################

metadata_check <- annotation %>%
  dplyr::distinct(
    sample_key,
    cell_line,
    condition_detected
  ) %>%
  dplyr::left_join(
    sample_information %>%
      dplyr::rename(
        expected_cell_line = cell_line,
        expected_condition = condition_detected
      ),
    by = "sample_key"
  ) %>%
  dplyr::mutate(
    metadata_matches = (
      cell_line == expected_cell_line &
        condition_detected == expected_condition
    )
  )

cat(
  "\n========================================\n",
  "Sample metadata validation\n",
  "========================================\n"
)

print(
  metadata_check
)

if (
  any(
    is.na(metadata_check$metadata_matches) |
    !metadata_check$metadata_matches
  )
) {
  
  stop(
    "Annotation assignments do not agree with sample_information."
  )
}


############################################################
# 6E. Attach GSM identifiers
############################################################

annotation <- annotation %>%
  dplyr::left_join(
    sample_information %>%
      dplyr::select(
        sample_key,
        GSM
      ),
    by = "sample_key"
  )


############################################################
# 6F. Confirm every annotation row received a GSM
############################################################

if (any(is.na(annotation$GSM))) {
  
  print(
    annotation %>%
      dplyr::filter(
        is.na(GSM)
      ) %>%
      dplyr::distinct(
        SampleID,
        sample_key,
        cell_line,
        condition_detected
      )
  )
  
  stop(
    "Some annotation rows did not receive a GSM identifier."
  )
}


############################################################
# 6G. Print final sample assignments
############################################################

cat(
  "\n========================================\n",
  "Final included samples\n",
  "========================================\n"
)

sample_assignment_check <- annotation %>%
  dplyr::distinct(
    SampleID,
    sample_id_clean,
    sample_key,
    GSM,
    cell_line,
    condition_detected
  ) %>%
  dplyr::arrange(
    factor(
      sample_key,
      levels = sample_information$sample_key
    )
  )

print(
  sample_assignment_check
)


############################################################
# 6H. Print cell numbers per sample
############################################################

cat(
  "\n========================================\n",
  "Annotated cells per sample\n",
  "========================================\n"
)

sample_cell_counts <- annotation %>%
  dplyr::count(
    sample_key,
    GSM,
    cell_line,
    condition_detected,
    name = "n_annotated_cells"
  ) %>%
  dplyr::arrange(
    factor(
      sample_key,
      levels = sample_information$sample_key
    )
  )

print(
  sample_cell_counts
)


############################################################
# 6I. Store sample order
############################################################

desired_sample_order <- sample_information$sample_key




############################################################
# 7. Load raw matrices and calculate sample-level
#    EN, IN, and AST average normalized expression
#
# Included conditions:
#   5-month in-vitro organoids
#   5-month transplanted organoids
#   6-month transplanted organoids
#   8-month transplanted organoids
#
# Excluded:
#   10-week organoids
############################################################


############################################################
# 7A. Confirm required objects and functions exist
############################################################

required_objects <- c(
  "annotation",
  "sample_information",
  "raw_directory"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1)
  )
]

if (length(missing_objects) > 0) {
  
  stop(
    "Missing required objects: ",
    paste(
      missing_objects,
      collapse = ", "
    )
  )
}

required_functions <- c(
  "extract_barcode_core",
  "find_matrix_files",
  "load_human_count_matrix"
)

missing_functions <- required_functions[
  !vapply(
    required_functions,
    exists,
    logical(1),
    mode = "function"
  )
]

if (length(missing_functions) > 0) {
  
  stop(
    "Missing required functions: ",
    paste(
      missing_functions,
      collapse = ", "
    )
  )
}


############################################################
# 7B. Confirm Step 6 produced the required columns
############################################################

required_sample_columns <- c(
  "sample_key",
  "GSM",
  "cell_line",
  "condition_detected"
)

missing_sample_columns <- setdiff(
  required_sample_columns,
  colnames(sample_information)
)

if (length(missing_sample_columns) > 0) {
  
  stop(
    "sample_information is missing: ",
    paste(
      missing_sample_columns,
      collapse = ", "
    )
  )
}

required_annotation_columns <- c(
  "barcode",
  "SampleID",
  "sample_key",
  "GSM",
  "cell_line",
  "condition_detected",
  "CellType"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  colnames(annotation)
)

if (length(missing_annotation_columns) > 0) {
  
  stop(
    "annotation is missing: ",
    paste(
      missing_annotation_columns,
      collapse = ", "
    )
  )
}


############################################################
# 7C. Define broad cell types from author annotations
#
# Author labels:
#   UL + DL = EN
#   iN      = IN
#   Ast     = AST
############################################################

desired_celltypes <- c(
  "EN",
  "IN",
  "AST"
)

annotation <- annotation %>%
  dplyr::mutate(
    
    barcode = trimws(
      as.character(barcode)
    ),
    
    barcode_core = extract_barcode_core(
      barcode
    ),
    
    author_celltype = trimws(
      as.character(CellType)
    ),
    
    author_celltype_clean = toupper(
      author_celltype
    ),
    
    broad_celltype = dplyr::case_when(
      
      author_celltype_clean %in%
        c(
          "UL",
          "DL"
        ) ~ "EN",
      
      author_celltype_clean == "IN" ~ "IN",
      
      author_celltype_clean == "AST" ~ "AST",
      
      TRUE ~ NA_character_
    )
  )


############################################################
# 7D. Create the annotation object used in the loop
#
# This is the object missing from the previous Step 7.
############################################################

included_annotation <- annotation %>%
  dplyr::filter(
    sample_key %in% sample_information$sample_key,
    !is.na(barcode_core),
    barcode_core != "",
    broad_celltype %in% desired_celltypes
  ) %>%
  dplyr::distinct(
    sample_key,
    barcode_core,
    .keep_all = TRUE
  )


############################################################
# 7E. Confirm every included sample has annotation
############################################################

annotation_count_check <- included_annotation %>%
  dplyr::count(
    sample_key,
    broad_celltype,
    name = "n_annotated_nuclei"
  ) %>%
  tidyr::complete(
    sample_key = sample_information$sample_key,
    broad_celltype = desired_celltypes,
    fill = list(
      n_annotated_nuclei = 0
    )
  ) %>%
  dplyr::arrange(
    factor(
      sample_key,
      levels = sample_information$sample_key
    ),
    factor(
      broad_celltype,
      levels = desired_celltypes
    )
  )

cat(
  "\n========================================\n",
  "Annotated EN, IN, and AST nuclei\n",
  "========================================\n"
)

print(
  annotation_count_check
)

samples_without_annotation <- annotation_count_check %>%
  dplyr::group_by(
    sample_key
  ) %>%
  dplyr::summarise(
    total_annotated = sum(
      n_annotated_nuclei
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    total_annotated == 0
  )

if (nrow(samples_without_annotation) > 0) {
  
  stop(
    "No EN, IN, or AST annotation was found for: ",
    paste(
      samples_without_annotation$sample_key,
      collapse = ", "
    )
  )
}


############################################################
# 7F. Initialize output objects
############################################################

sample_average_expression <- list()

sample_nuclei_numbers <- numeric(0)

sample_matching_summary <- vector(
  mode = "list",
  length = nrow(sample_information)
)


############################################################
# 7G. Load and process every sample
############################################################

for (sample_index in seq_len(nrow(sample_information))) {
  
  current_sample <- sample_information[
    sample_index,
    ,
    drop = FALSE
  ]
  
  current_key <- as.character(
    current_sample$sample_key
  )
  
  current_GSM <- as.character(
    current_sample$GSM
  )
  
  current_cell_line <- as.character(
    current_sample$cell_line
  )
  
  # Correct column name from Step 6
  current_condition <- as.character(
    current_sample$condition_detected
  )
  
  cat(
    "\n========================================\n",
    "Processing ",
    current_key,
    " (",
    current_GSM,
    ")\n",
    "========================================\n",
    sep = ""
  )
  
  
  ##########################################################
  # Select EN, IN, and AST annotations for this sample
  ##########################################################
  
  current_annotation <- included_annotation %>%
    dplyr::filter(
      sample_key == current_key
    )
  
  number_annotated <- nrow(
    current_annotation
  )
  
  if (number_annotated == 0) {
    
    stop(
      "No included EN, IN, or AST annotation rows were found for ",
      current_key,
      "."
    )
  }
  
  
  ##########################################################
  # Locate and load the raw count matrix
  ##########################################################
  
  sample_files <- find_matrix_files(
    GSM = current_GSM,
    raw_directory = raw_directory
  )
  
  counts <- load_human_count_matrix(
    matrix_file = sample_files$matrix,
    barcode_file = sample_files$barcodes,
    feature_file = sample_files$features
  )
  
  if (is.null(colnames(counts))) {
    
    stop(
      "The raw matrix for ",
      current_key,
      " does not have cell barcodes."
    )
  }
  
  
  ##########################################################
  # Match author-annotated barcodes to the raw matrix
  ##########################################################
  
  raw_barcode_core <- extract_barcode_core(
    colnames(counts)
  )
  
  matched_column_indices <- match(
    current_annotation$barcode_core,
    raw_barcode_core
  )
  
  successfully_matched <- !is.na(
    matched_column_indices
  )
  
  matched_annotation <- current_annotation[
    successfully_matched,
    ,
    drop = FALSE
  ]
  
  matched_column_indices <- matched_column_indices[
    successfully_matched
  ]
  
  number_matched <- length(
    matched_column_indices
  )
  
  matching_fraction <- number_matched /
    number_annotated
  
  message(
    "Annotated EN/IN/AST nuclei: ",
    number_annotated
  )
  
  message(
    "Matched to raw matrix: ",
    number_matched,
    " (",
    round(
      100 * matching_fraction,
      digits = 1
    ),
    "%)"
  )
  
  if (
    number_matched < 30 ||
    matching_fraction < 0.5
  ) {
    
    stop(
      "Too few annotated nuclei matched the raw matrix for ",
      current_key,
      "."
    )
  }
  
  
  ##########################################################
  # Subset raw counts to matched annotated nuclei
  ##########################################################
  
  counts <- counts[
    ,
    matched_column_indices,
    drop = FALSE
  ]
  
  unique_cell_names <- paste(
    current_key,
    matched_annotation$barcode_core,
    sep = "__"
  )
  
  colnames(counts) <- unique_cell_names
  
  
  ##########################################################
  # Create Seurat metadata
  ##########################################################
  
  current_metadata <- data.frame(
    
    row.names = unique_cell_names,
    
    sample_key = rep(
      current_key,
      number_matched
    ),
    
    GSM = rep(
      current_GSM,
      number_matched
    ),
    
    cell_line = rep(
      current_cell_line,
      number_matched
    ),
    
    condition_detected = rep(
      current_condition,
      number_matched
    ),
    
    author_celltype = as.character(
      matched_annotation$author_celltype
    ),
    
    broad_celltype = as.character(
      matched_annotation$broad_celltype
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  ##########################################################
  # Create and normalize the sample-level Seurat object
  ##########################################################
  
  current_object <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = current_metadata,
    project = current_key,
    min.cells = 0,
    min.features = 0
  )
  
  Seurat::DefaultAssay(
    current_object
  ) <- "RNA"
  
  current_object <- Seurat::NormalizeData(
    object = current_object,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  
  ##########################################################
  # Count matched nuclei by broad cell type
  ##########################################################
  
  current_cell_counts <- table(
    factor(
      current_object$broad_celltype,
      levels = desired_celltypes
    )
  )
  
  cat(
    "\nMatched nuclei by broad cell type:\n"
  )
  
  print(
    current_cell_counts
  )
  
  missing_current_celltypes <- names(
    current_cell_counts
  )[
    current_cell_counts == 0
  ]
  
  if (length(missing_current_celltypes) > 0) {
    
    stop(
      current_key,
      " is missing the following cell types after barcode matching: ",
      paste(
        missing_current_celltypes,
        collapse = ", "
      )
    )
  }
  
  for (current_celltype in desired_celltypes) {
    
    count_key <- paste(
      current_key,
      current_celltype,
      sep = "__"
    )
    
    sample_nuclei_numbers[
      count_key
    ] <- as.numeric(
      current_cell_counts[
        current_celltype
      ]
    )
  }
  
  
  ##########################################################
  # Calculate average normalized expression
  #
  # AverageExpression converts log-normalized values back
  # into linear normalized expression before averaging.
  ##########################################################
  
  average_result <- tryCatch(
    
    Seurat::AverageExpression(
      object = current_object,
      assays = "RNA",
      layer = "data",
      group.by = "broad_celltype",
      return.seurat = FALSE,
      verbose = FALSE
    ),
    
    error = function(e) {
      
      message(
        "Retrying AverageExpression using slot = 'data'."
      )
      
      Seurat::AverageExpression(
        object = current_object,
        assays = "RNA",
        slot = "data",
        group.by = "broad_celltype",
        return.seurat = FALSE,
        verbose = FALSE
      )
    }
  )
  
  average_matrix <- average_result$RNA
  
  average_matrix <- as.matrix(
    average_matrix
  )
  
  missing_average_celltypes <- setdiff(
    desired_celltypes,
    colnames(average_matrix)
  )
  
  if (length(missing_average_celltypes) > 0) {
    
    stop(
      "AverageExpression is missing ",
      paste(
        missing_average_celltypes,
        collapse = ", "
      ),
      " for ",
      current_key,
      ". Detected columns: ",
      paste(
        colnames(average_matrix),
        collapse = ", "
      )
    )
  }
  
  
  ##########################################################
  # Store expression vectors for this sample
  ##########################################################
  
  sample_average_expression[[current_key]] <- setNames(
    
    lapply(
      desired_celltypes,
      function(current_celltype) {
        
        setNames(
          as.numeric(
            average_matrix[
              ,
              current_celltype
            ]
          ),
          rownames(
            average_matrix
          )
        )
      }
    ),
    
    desired_celltypes
  )
  
  
  ##########################################################
  # Store sample matching summary
  ##########################################################
  
  sample_matching_summary[[sample_index]] <- data.frame(
    
    sample_key = current_key,
    
    GSM = current_GSM,
    
    cell_line = current_cell_line,
    
    condition_detected = current_condition,
    
    n_annotated_EN_IN_AST = number_annotated,
    
    n_matched_EN_IN_AST = number_matched,
    
    matching_fraction = matching_fraction,
    
    n_EN = as.numeric(
      current_cell_counts["EN"]
    ),
    
    n_IN = as.numeric(
      current_cell_counts["IN"]
    ),
    
    n_AST = as.numeric(
      current_cell_counts["AST"]
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  ##########################################################
  # Remove temporary objects before the next sample
  ##########################################################
  
  rm(
    counts,
    current_object,
    average_result,
    average_matrix
  )
  
  invisible(
    gc()
  )
}


############################################################
# 7H. Combine and print sample matching summaries
############################################################

sample_matching_summary <- dplyr::bind_rows(
  sample_matching_summary
)

cat(
  "\n========================================\n",
  "Final barcode-matching summary\n",
  "========================================\n"
)

print(
  sample_matching_summary
)


############################################################
# 7I. Print matched nuclei available for pooling
############################################################

sample_nuclei_table <- data.frame(
  
  sample_celltype = names(
    sample_nuclei_numbers
  ),
  
  n_nuclei = as.numeric(
    sample_nuclei_numbers
  ),
  
  stringsAsFactors = FALSE
)

cat(
  "\n========================================\n",
  "Matched nuclei used for average expression\n",
  "========================================\n"
)

print(
  sample_nuclei_table
)

############################################################
# 8. Validate the sample-level expression results
############################################################

required_step7_objects <- c(
  "sample_average_expression",
  "sample_nuclei_numbers",
  "sample_matching_summary",
  "sample_information"
)

missing_step7_objects <- required_step7_objects[
  !vapply(
    required_step7_objects,
    exists,
    logical(1)
  )
]

if (length(missing_step7_objects) > 0) {
  stop(
    "Missing Step 7 objects: ",
    paste(missing_step7_objects, collapse = ", ")
  )
}

required_sample_information_columns <- c(
  "sample_key",
  "GSM",
  "cell_line",
  "condition_detected"
)

missing_sample_information_columns <- setdiff(
  required_sample_information_columns,
  colnames(sample_information)
)

if (length(missing_sample_information_columns) > 0) {
  stop(
    "sample_information is missing: ",
    paste(missing_sample_information_columns, collapse = ", ")
  )
}

desired_celltypes <- c(
  "EN",
  "IN",
  "AST"
)

condition_order <- c(
  "5m_in_vitro",
  "5m_transplant",
  "6m_transplant",
  "8m_transplant"
)

missing_conditions <- setdiff(
  condition_order,
  unique(sample_information$condition_detected)
)

if (length(missing_conditions) > 0) {
  stop(
    "sample_information is missing these conditions: ",
    paste(missing_conditions, collapse = ", ")
  )
}

missing_sample_expression <- setdiff(
  sample_information$sample_key,
  names(sample_average_expression)
)

if (length(missing_sample_expression) > 0) {
  stop(
    "sample_average_expression is missing: ",
    paste(missing_sample_expression, collapse = ", ")
  )
}






############################################################
# 9. Standardize human gene names
#
# Important:
# Gene names are captured before numeric conversion so that
# as.numeric() cannot remove them.
############################################################

standardize_named_expression_vector <- function(
    expression_vector,
    vector_name
) {
  
  ##########################################################
  # Preserve gene names before numeric conversion
  ##########################################################
  
  original_gene_names <- names(expression_vector)
  
  if (is.null(original_gene_names)) {
    stop(
      vector_name,
      " does not have gene names."
    )
  }
  
  if (length(original_gene_names) != length(expression_vector)) {
    stop(
      vector_name,
      " has inconsistent gene names and expression values."
    )
  }
  
  expression_values <- suppressWarnings(
    as.numeric(expression_vector)
  )
  
  gene_names <- trimws(
    as.character(original_gene_names)
  )
  
  ##########################################################
  # Identify human and mouse prefixes
  ##########################################################
  
  human_prefixed <- grepl(
    "^(hg19|hg38|grch37|grch38|human)[-_:]",
    gene_names,
    ignore.case = TRUE
  )
  
  mouse_prefixed <- grepl(
    "^(mm10|mm39|mouse)[-_:]",
    gene_names,
    ignore.case = TRUE
  )
  
  ##########################################################
  # Retain human genes
  #
  # If both human and mouse prefixes exist, keep only human.
  # If mouse prefixes exist without human prefixes, remove
  # mouse-prefixed genes and retain unprefixed genes.
  # If no mouse prefixes exist, retain all genes.
  ##########################################################
  
  if (any(mouse_prefixed) && any(human_prefixed)) {
    
    genes_to_keep <- human_prefixed
    
  } else if (any(mouse_prefixed)) {
    
    genes_to_keep <- !mouse_prefixed
    
  } else {
    
    genes_to_keep <- rep(
      TRUE,
      length(gene_names)
    )
  }
  
  gene_names <- gene_names[genes_to_keep]
  expression_values <- expression_values[genes_to_keep]
  
  ##########################################################
  # Remove genome-reference prefixes
  ##########################################################
  
  gene_names <- gsub(
    "^(hg19|hg38|grch37|grch38|human)[-_:]",
    "",
    gene_names,
    ignore.case = TRUE
  )
  
  gene_names <- toupper(
    trimws(gene_names)
  )
  
  ##########################################################
  # Remove unusable entries
  ##########################################################
  
  valid_entries <- (
    !is.na(gene_names) &
      gene_names != "" &
      gene_names != "NA" &
      is.finite(expression_values)
  )
  
  gene_names <- gene_names[valid_entries]
  expression_values <- expression_values[valid_entries]
  
  if (length(gene_names) == 0) {
    stop(
      vector_name,
      " contains no usable human genes after cleaning."
    )
  }
  
  ##########################################################
  # Collapse duplicated gene symbols by mean expression
  ##########################################################
  
  cleaned_table <- data.frame(
    gene = gene_names,
    expression = expression_values,
    stringsAsFactors = FALSE
  )
  
  cleaned_table <- cleaned_table %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      expression = mean(
        expression,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(gene)
  
  cleaned_vector <- cleaned_table$expression
  names(cleaned_vector) <- cleaned_table$gene
  
  return(cleaned_vector)
}


############################################################
# 10. Clean every sample-by-cell-type expression vector
############################################################

expression_vector_list <- list()
expression_vector_qc <- list()

qc_index <- 1L

for (current_key in sample_information$sample_key) {
  
  if (!current_key %in% names(sample_average_expression)) {
    stop(
      "No average-expression result was found for ",
      current_key,
      "."
    )
  }
  
  current_sample_result <- sample_average_expression[[current_key]]
  
  missing_current_celltypes <- setdiff(
    desired_celltypes,
    names(current_sample_result)
  )
  
  if (length(missing_current_celltypes) > 0) {
    stop(
      current_key,
      " is missing average expression for: ",
      paste(
        missing_current_celltypes,
        collapse = ", "
      )
    )
  }
  
  for (current_celltype in desired_celltypes) {
    
    vector_key <- paste(
      current_key,
      current_celltype,
      sep = "__"
    )
    
    current_vector <- current_sample_result[[current_celltype]]
    
    cleaned_vector <- standardize_named_expression_vector(
      expression_vector = current_vector,
      vector_name = vector_key
    )
    
    expression_vector_list[[vector_key]] <- cleaned_vector
    
    expression_vector_qc[[qc_index]] <- data.frame(
      vector_key = vector_key,
      n_clean_genes = length(cleaned_vector),
      first_gene = names(cleaned_vector)[1],
      last_gene = names(cleaned_vector)[length(cleaned_vector)],
      stringsAsFactors = FALSE
    )
    
    qc_index <- qc_index + 1L
  }
}

expression_vector_qc <- dplyr::bind_rows(
  expression_vector_qc
)

cat(
  "\n========================================\n",
  "Cleaned expression-vector summary\n",
  "========================================\n"
)

print(
  expression_vector_qc,
  row.names = FALSE
)


############################################################
# 11. Define and validate the common gene universe
############################################################

############################################################
# 11A. Confirm that no vector is empty
############################################################

empty_vectors <- names(expression_vector_list)[
  vapply(
    expression_vector_list,
    length,
    integer(1)
  ) == 0
]

if (length(empty_vectors) > 0) {
  stop(
    "The following expression vectors contain no genes: ",
    paste(
      empty_vectors,
      collapse = ", "
    )
  )
}


############################################################
# 11B. Calculate the gene intersection
############################################################

gene_lists <- lapply(
  expression_vector_list,
  names
)

common_genes <- Reduce(
  intersect,
  gene_lists
)


############################################################
# 11C. Print overlap diagnostics if the intersection fails
############################################################

if (length(common_genes) == 0) {
  
  first_vector_name <- names(expression_vector_list)[1]
  first_gene_list <- names(expression_vector_list[[1]])
  
  overlap_diagnostic <- data.frame(
    vector_key = names(expression_vector_list),
    
    n_genes = vapply(
      expression_vector_list,
      length,
      integer(1)
    ),
    
    overlap_with_first_vector = vapply(
      expression_vector_list,
      function(current_vector) {
        length(
          intersect(
            first_gene_list,
            names(current_vector)
          )
        )
      },
      integer(1)
    ),
    
    stringsAsFactors = FALSE
  )
  
  cat(
    "\nGene overlap relative to ",
    first_vector_name,
    ":\n",
    sep = ""
  )
  
  print(
    overlap_diagnostic,
    row.names = FALSE
  )
  
  stop(
    "No common genes remained after standardization. ",
    "Review the overlap table printed above."
  )
}


############################################################
# 11D. Preserve a stable gene order
############################################################

first_vector_gene_order <- names(
  expression_vector_list[[1]]
)

common_genes <- first_vector_gene_order[
  first_vector_gene_order %in% common_genes
]


############################################################
# 11E. Require a reasonable number of shared genes
############################################################

if (length(common_genes) < 5000) {
  warning(
    "Only ",
    length(common_genes),
    " genes were shared across every sample and cell type."
  )
}

if (length(common_genes) < 500) {
  stop(
    "Only ",
    length(common_genes),
    " common genes remained, which is too few."
  )
}


############################################################
# 11F. Put every vector into the identical gene order
############################################################

expression_vector_list <- lapply(
  expression_vector_list,
  function(current_vector) {
    current_vector[common_genes]
  }
)


############################################################
# 11G. Confirm that every vector now has identical genes
############################################################

identical_gene_order <- vapply(
  expression_vector_list,
  function(current_vector) {
    identical(
      names(current_vector),
      common_genes
    )
  },
  logical(1)
)

if (!all(identical_gene_order)) {
  stop(
    "The expression vectors do not have identical gene order ",
    "after common-gene subsetting."
  )
}

cat(
  "\n========================================\n",
  "Common gene universe\n",
  "========================================\n",
  "Expression vectors: ",
  length(expression_vector_list),
  "\n",
  "Genes shared by every vector: ",
  length(common_genes),
  "\n",
  "Identical gene order in every vector: TRUE\n",
  sep = ""
)







############################################################
# 12. Validate nuclei-number weights
############################################################

required_weight_names <- unlist(
  lapply(
    sample_information$sample_key,
    function(current_key) {
      paste(
        current_key,
        desired_celltypes,
        sep = "__"
      )
    }
  ),
  use.names = FALSE
)

missing_weight_names <- setdiff(
  required_weight_names,
  names(sample_nuclei_numbers)
)

if (length(missing_weight_names) > 0) {
  stop(
    "sample_nuclei_numbers is missing: ",
    paste(missing_weight_names, collapse = ", ")
  )
}

invalid_weight_names <- required_weight_names[
  !is.finite(sample_nuclei_numbers[required_weight_names]) |
    sample_nuclei_numbers[required_weight_names] <= 0
]

if (length(invalid_weight_names) > 0) {
  stop(
    "Invalid nuclei-number weights for: ",
    paste(invalid_weight_names, collapse = ", ")
  )
}


############################################################
# 13. Pool HUES6 and iP822 within each condition
#
# Weighted pooling is used:
#
#   pooled expression =
#     sum(sample mean × sample nuclei) /
#     sum(sample nuclei)
#
# This is equivalent to averaging all nuclei after using the
# same normalization scale for every sample.
############################################################

condition_name_prefix <- c(
  "5m_in_vitro" = "GSE185472_5m_in_vitro",
  "5m_transplant" = "GSE185472_5m_transplant",
  "6m_transplant" = "GSE185472_6m_transplant",
  "8m_transplant" = "GSE185472_8m_transplant"
)

pooled_expression_vectors <- list()

pooling_summary <- list()

pooling_summary_index <- 1L

for (current_celltype in desired_celltypes) {
  
  for (current_condition in condition_order) {
    
    current_sample_keys <- sample_information$sample_key[
      sample_information$condition_detected == current_condition
    ]
    
    if (length(current_sample_keys) == 0) {
      stop(
        "No samples were found for ",
        current_condition,
        "."
      )
    }
    
    current_vector_keys <- paste(
      current_sample_keys,
      current_celltype,
      sep = "__"
    )
    
    missing_vector_keys <- setdiff(
      current_vector_keys,
      names(expression_vector_list)
    )
    
    if (length(missing_vector_keys) > 0) {
      stop(
        "Missing expression vectors: ",
        paste(missing_vector_keys, collapse = ", ")
      )
    }
    
    current_weights <- as.numeric(
      sample_nuclei_numbers[current_vector_keys]
    )
    
    names(current_weights) <- current_sample_keys
    
    current_expression_matrix <- do.call(
      cbind,
      expression_vector_list[current_vector_keys]
    )
    
    current_expression_matrix <- as.matrix(
      current_expression_matrix
    )
    
    storage.mode(current_expression_matrix) <- "double"
    
    rownames(current_expression_matrix) <- common_genes
    
    colnames(current_expression_matrix) <- current_sample_keys
    
    finite_expression <- apply(
      current_expression_matrix,
      1,
      function(x) {
        all(is.finite(x))
      }
    )
    
    if (!all(finite_expression)) {
      stop(
        sum(!finite_expression),
        " genes contained non-finite expression values for ",
        current_condition,
        " ",
        current_celltype,
        "."
      )
    }
    
    pooled_vector <- as.numeric(
      current_expression_matrix %*% current_weights /
        sum(current_weights)
    )
    
    names(pooled_vector) <- common_genes
    
    output_column_name <- paste0(
      unname(condition_name_prefix[current_condition]),
      "_",
      current_celltype
    )
    
    pooled_expression_vectors[[output_column_name]] <- pooled_vector
    
    for (current_sample_key in current_sample_keys) {
      
      current_weight_key <- paste(
        current_sample_key,
        current_celltype,
        sep = "__"
      )
      
      pooling_summary[[pooling_summary_index]] <- data.frame(
        condition = current_condition,
        broad_celltype = current_celltype,
        sample_key = current_sample_key,
        cell_line = sample_information$cell_line[
          sample_information$sample_key == current_sample_key
        ],
        n_nuclei = as.numeric(
          sample_nuclei_numbers[current_weight_key]
        ),
        stringsAsFactors = FALSE
      )
      
      pooling_summary_index <- pooling_summary_index + 1L
    }
  }
}


############################################################
# 14. Set the final expression-column order
############################################################

final_expression_column_order <- c(
  "GSE185472_5m_in_vitro_EN",
  "GSE185472_5m_transplant_EN",
  "GSE185472_6m_transplant_EN",
  "GSE185472_8m_transplant_EN",
  
  "GSE185472_5m_in_vitro_IN",
  "GSE185472_5m_transplant_IN",
  "GSE185472_6m_transplant_IN",
  "GSE185472_8m_transplant_IN",
  
  "GSE185472_5m_in_vitro_AST",
  "GSE185472_5m_transplant_AST",
  "GSE185472_6m_transplant_AST",
  "GSE185472_8m_transplant_AST"
)

missing_final_columns <- setdiff(
  final_expression_column_order,
  names(pooled_expression_vectors)
)

if (length(missing_final_columns) > 0) {
  stop(
    "The following pooled-expression columns are missing: ",
    paste(missing_final_columns, collapse = ", ")
  )
}


############################################################
# 15. Create the final expression table
############################################################

GSE185472_invitro_transplant_EN_IN_AST_expression <- data.frame(
  gene = common_genes,
  stringsAsFactors = FALSE
)

for (current_column in final_expression_column_order) {
  GSE185472_invitro_transplant_EN_IN_AST_expression[[current_column]] <- pooled_expression_vectors[[current_column]]
}


############################################################
# 16. Final quality-control checks
############################################################

if (anyDuplicated(GSE185472_invitro_transplant_EN_IN_AST_expression$gene) > 0) {
  stop(
    "The final expression table contains duplicated genes."
  )
}

remaining_human_prefixes <- sum(
  grepl(
    "^(HG19|HG38|GRCH37|GRCH38)[-_:]",
    GSE185472_invitro_transplant_EN_IN_AST_expression$gene,
    ignore.case = TRUE
  )
)

remaining_mouse_prefixes <- sum(
  grepl(
    "^(MM10|MM39|MOUSE)[-_:]",
    GSE185472_invitro_transplant_EN_IN_AST_expression$gene,
    ignore.case = TRUE
  )
)

if (remaining_human_prefixes > 0) {
  stop(
    remaining_human_prefixes,
    " genome-reference prefixes remain in the final gene column."
  )
}

if (remaining_mouse_prefixes > 0) {
  stop(
    remaining_mouse_prefixes,
    " mouse-prefixed genes remain in the final gene column."
  )
}

expression_only_matrix <- as.matrix(
  GSE185472_invitro_transplant_EN_IN_AST_expression[
    ,
    final_expression_column_order,
    drop = FALSE
  ]
)

if (any(!is.finite(expression_only_matrix))) {
  stop(
    "The final expression table contains non-finite values."
  )
}


############################################################
# 17. Print nuclei used for weighted pooling
############################################################

pooling_summary <- dplyr::bind_rows(
  pooling_summary
)

pooling_summary <- pooling_summary %>%
  dplyr::arrange(
    factor(
      broad_celltype,
      levels = desired_celltypes
    ),
    factor(
      condition,
      levels = condition_order
    ),
    cell_line
  )

cat(
  "\n========================================\n",
  "Nuclei used for weighted pooling\n",
  "========================================\n"
)

print(
  pooling_summary
)

pooled_nuclei_summary <- pooling_summary %>%
  dplyr::group_by(
    condition,
    broad_celltype
  ) %>%
  dplyr::summarise(
    n_cell_lines = dplyr::n_distinct(cell_line),
    total_nuclei = sum(n_nuclei),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    factor(
      broad_celltype,
      levels = desired_celltypes
    ),
    factor(
      condition,
      levels = condition_order
    )
  )

cat(
  "\n========================================\n",
  "Pooled nuclei totals\n",
  "========================================\n"
)

print(
  pooled_nuclei_summary
)


############################################################
# 18. Save the final CSV
############################################################

final_csv_filename <- paste0(
  "GSE185472_5m_invitro_",
  "5m_6m_8m_transplant_",
  "EN_IN_AST_expression.csv"
)

write.csv(
  GSE185472_invitro_transplant_EN_IN_AST_expression,
  file = final_csv_filename,
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 19. Print the final result directly in R
############################################################

cat(
  "\n========================================\n",
  "Final expression table\n",
  "========================================\n"
)

print(
  head(
    GSE185472_invitro_transplant_EN_IN_AST_expression
  )
)

cat(
  "\nDimensions: ",
  nrow(GSE185472_invitro_transplant_EN_IN_AST_expression),
  " genes × ",
  ncol(GSE185472_invitro_transplant_EN_IN_AST_expression) - 1,
  " expression profiles\n",
  sep = ""
)

cat(
  "\nFinal columns:\n"
)

print(
  colnames(
    GSE185472_invitro_transplant_EN_IN_AST_expression
  )
)

cat(
  "\nCSV saved to:\n",
  normalizePath(final_csv_filename),
  "\n",
  sep = ""
)







#compare using person correlation

xeno_expression <-as.data.frame(xeno_expression)

GSE190815_expression <-as.data.frame(GSE190815_expression)

Park_iMicro_expression <-as.data.frame(Park_iMicro_expression)

Gage_MG_expression <-as.data.frame(Gage_MG_expression)

CHOOSE_CTRL_expression <-as.data.frame(CHOOSE_CTRL_expression)

GSE168408_broad_age_celltype_average_expression <-
  as.data.frame(GSE168408_broad_age_celltype_average_expression)

GSE185472_5m_invitro_5m_6m_8m_transplant_EN_IN_AST_expression <- as.data.frame(GSE185472_5m_invitro_5m_6m_8m_transplant_EN_IN_AST_expression)

head(xeno_expression)
head(GSE168408_broad_age_celltype_average_expression)#reference
head(GSE190815_expression)
head(Gage_MG_expression)
head(Park_iMicro_expression)
head(CHOOSE_CTRL_expression)
head(GSE185472_5m_invitro_5m_6m_8m_transplant_EN_IN_AST_expression)



############################################################
# Publication workflow
# Model versus human developmental-age Spearman correlation
#
# - Human reference: GSE168408
# - Gene symbols used exactly as supplied
# - No marker/HVG/technical-gene filtering
# - Zero values retained
# - One identical shared-gene set per cell type
# - Nothing written to disk
############################################################

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("Install the pheatmap package before running this workflow.")
}

############################################################
# 1. Register input tables
############################################################

raw_expression_tables <- list(
  
  reference =
    GSE168408_broad_age_celltype_average_expression,
  
  HuMiNAX =
    xeno_expression,
  
  GSE190815 =
    GSE190815_expression,
  
  GSE185472 =
    GSE185472_5m_invitro_5m_6m_8m_transplant_EN_IN_AST_expression,
  
  Gage =
    Gage_MG_expression,
  
  CHOOSE =
    CHOOSE_CTRL_expression,
  
  Park =
    Park_iMicro_expression
)


############################################################
# 2. Clean each input table once
#
# This step:
# - keeps gene symbols exactly as provided
# - removes blank gene names
# - converts expression columns to numeric
# - collapses duplicate gene symbols by their mean
#
# It does not select or remove genes based on expression.
############################################################

clean_expression_table <- function(
    expression_table,
    object_name
) {
  
  expression_table <- data.frame(
    expression_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!"gene" %in% colnames(expression_table)) {
    stop(
      object_name,
      " does not contain a column named 'gene'."
    )
  }
  
  expression_columns <- setdiff(
    colnames(expression_table),
    "gene"
  )
  
  if (length(expression_columns) == 0) {
    stop(
      object_name,
      " contains no expression columns."
    )
  }
  
  expression_table$gene <- trimws(
    as.character(expression_table$gene)
  )
  
  expression_table <- expression_table[
    !is.na(expression_table$gene) &
      expression_table$gene != "",
    c("gene", expression_columns),
    drop = FALSE
  ]
  
  expression_table[, expression_columns] <- lapply(
    expression_table[
      ,
      expression_columns,
      drop = FALSE
    ],
    function(values) {
      suppressWarnings(
        as.numeric(values)
      )
    }
  )
  
  cleaned_table <- stats::aggregate(
    expression_table[
      ,
      expression_columns,
      drop = FALSE
    ],
    by = list(
      gene = expression_table$gene
    ),
    FUN = function(values) {
      
      values <- values[
        is.finite(values)
      ]
      
      if (length(values) == 0) {
        NA_real_
      } else {
        mean(values)
      }
    }
  )
  
  colnames(cleaned_table) <- c(
    "gene",
    expression_columns
  )
  
  cleaned_table
}

cleaned_expression_tables <- setNames(
  lapply(
    names(raw_expression_tables),
    function(object_name) {
      
      clean_expression_table(
        expression_table =
          getElement(
            raw_expression_tables,
            object_name
          ),
        object_name = object_name
      )
    }
  ),
  names(raw_expression_tables)
)

get_clean_table <- function(dataset_name) {
  
  getElement(
    cleaned_expression_tables,
    dataset_name
  )
}

############################################################
# 3. Define every model column explicitly
############################################################

model_specifications <- data.frame(
  
  cell_type = c(
    rep("EN", 9),
    rep("IN", 8),
    rep("AST", 9),
    rep("MG", 6)
  ),
  
  dataset = c(
    
    # EN
    "HuMiNAX",
    "HuMiNAX",
    "GSE190815",
    "GSE190815",
    "CHOOSE",
    "GSE185472",
    "GSE185472",
    "GSE185472",
    "GSE185472",
    
    # IN
    "HuMiNAX",
    "HuMiNAX",
    "GSE190815",
    "CHOOSE",
    "GSE185472",
    "GSE185472",
    "GSE185472",
    "GSE185472",
    
    # AST
    "HuMiNAX",
    "HuMiNAX",
    "GSE190815",
    "GSE190815",
    "CHOOSE",
    "GSE185472",
    "GSE185472",
    "GSE185472",
    "GSE185472",
    
    # MG
    "HuMiNAX",
    "HuMiNAX",
    "Gage",
    "Gage",
    "Gage",
    "Park"
  ),
  
  expression_column = c(
    
    # EN
    "3R_EN",
    "4R_EN",
    "hCS_EN",
    "t-hCS_EN",
    "WT_organoid_EN",
    "GSE185472_5m_in_vitro_EN",
    "GSE185472_5m_transplant_EN",
    "GSE185472_6m_transplant_EN",
    "GSE185472_8m_transplant_EN",
    
    # IN
    "3R_IN",
    "4R_IN",
    "hCS_IN",
    "WT_organoid_IN",
    "GSE185472_5m_in_vitro_IN",
    "GSE185472_5m_transplant_IN",
    "GSE185472_6m_transplant_IN",
    "GSE185472_8m_transplant_IN",
    
    # AST
    "3R_AST",
    "4R_AST",
    "hCS_AST",
    "t-hCS_AST",
    "WT_organoid_AST",
    "GSE185472_5m_in_vitro_AST",
    "GSE185472_5m_transplant_AST",
    "GSE185472_6m_transplant_AST",
    "GSE185472_8m_transplant_AST",
    
    # MG
    "3R_MG",
    "4R_MG",
    "Gage_6w_MG",
    "Gage_12w_MG",
    "Gage_24w_MG",
    "Park_iMicro"
  ),
  
  heatmap_label = c(
    
    # EN: 9 rows
    "3RWT 4mo",
    "4RWT 4mo",
    "hCS in vitro 8mo",
    "hCS transplant 8mo",
    "in vitro 4mo",
    "in vitro 5mo",
    "transplant 5mo",
    "transplant 6mo",
    "transplant 8mo",
    
    # IN: 8 rows
    "3RWT 4mo",
    "4RWT 4mo",
    "hCS in vitro 8mo",
    "in vitro 4mo",
    "in vitro 5mo",
    "transplant 5mo",
    "transplant 6mo",
    "transplant 8mo",
    
    # AST: 9 rows
    "3RWT 4mo",
    "4RWT 4mo",
    "hCS in vitro 8mo",
    "hCS transplant 8mo",
    "in vitro 4mo",
    "in vitro 5mo",
    "transplant 5mo",
    "transplant 6mo",
    "transplant 8mo",
    
    # MG: 6 rows
    "3RWT 4mo",
    "4RWT 4mo",
    "transplant 6w",
    "transplant 12w",
    "transplant 24w",
    "MC-HBOs"
  ),

  stringsAsFactors = FALSE,
  check.names = FALSE
)



############################################################
# 4. Define human developmental-age columns
############################################################

cell_types <- c(
  "EN",
  "IN",
  "AST",
  "MG"
)

age_keys <- c(
  "Prenatal",
  "Neonatal",
  "Early_infancy",
  "Late_infancy",
  "Early_childhood",
  "Childhood",
  "Adolescence",
  "Adult"
)

age_labels <- c(
  "Prenatal",
  "Neonatal",
  "Early infancy",
  "Late infancy",
  "Early childhood",
  "Childhood",
  "Adolescence",
  "Adult"
)

reference_clean <- get_clean_table(
  "reference"
)

############################################################
# 5. Validate every requested column before analysis
############################################################

expected_reference_columns <- unlist(
  lapply(
    cell_types,
    function(cell_type) {
      
      paste0(
        age_keys,
        "_",
        cell_type
      )
    }
  ),
  use.names = FALSE
)

missing_reference_columns <- setdiff(
  expected_reference_columns,
  colnames(reference_clean)
)

if (length(missing_reference_columns) > 0) {
  
  stop(
    "The reference table is missing these columns: ",
    paste(
      missing_reference_columns,
      collapse = ", "
    )
  )
}

for (
  row_number in seq_len(
    nrow(model_specifications)
  )
) {
  
  dataset_name <-
    model_specifications$dataset[
      row_number
    ]
  
  expression_column <-
    model_specifications$expression_column[
      row_number
    ]
  
  if (
    !dataset_name %in%
    names(cleaned_expression_tables)
  ) {
    
    stop(
      "Dataset '",
      dataset_name,
      "' was not registered."
    )
  }
  
  current_table <- get_clean_table(
    dataset_name
  )
  
  if (
    !expression_column %in%
    colnames(current_table)
  ) {
    
    stop(
      "Column '",
      expression_column,
      "' is missing from dataset '",
      dataset_name,
      "'."
    )
  }
}

############################################################
# 6. Calculate one Spearman correlation matrix per cell type
#
# For each cell type, the common gene set is the intersection
# across:
#
# - every included model for that cell type
# - all eight human developmental-age reference columns
#
# A gene with NA, NaN, or Inf in any one input is removed from
# every comparison for that cell type.
#
# Zero-expression values are retained.
############################################################

calculate_celltype_correlations <- function(
    cell_type
) {
  
  specifications <- model_specifications[
    model_specifications$cell_type ==
      cell_type,
    ,
    drop = FALSE
  ]
  
  if (nrow(specifications) == 0) {
    
    stop(
      "No model columns were specified for ",
      cell_type,
      "."
    )
  }
  
  reference_columns <- paste0(
    age_keys,
    "_",
    cell_type
  )
  
  gene_lists <- c(
    
    list(
      reference_clean$gene
    ),
    
    lapply(
      seq_len(
        nrow(specifications)
      ),
      function(row_number) {
        
        current_table <- get_clean_table(
          specifications$dataset[
            row_number
          ]
        )
        
        current_table$gene
      }
    )
  )
  
  shared_gene_set <- Reduce(
    intersect,
    gene_lists
  )
  
  # Preserve the ordering of the reference gene table.
  common_genes <- reference_clean$gene[
    reference_clean$gene %in%
      shared_gene_set
  ]
  
  if (length(common_genes) < 2) {
    
    stop(
      cell_type,
      " has fewer than two shared genes."
    )
  }
  
  reference_rows <- match(
    common_genes,
    reference_clean$gene
  )
  
  reference_matrix <- as.matrix(
    reference_clean[
      reference_rows,
      reference_columns,
      drop = FALSE
    ]
  )
  
  storage.mode(
    reference_matrix
  ) <- "double"
  
  rownames(
    reference_matrix
  ) <- common_genes
  
  colnames(
    reference_matrix
  ) <- age_labels
  
  model_matrix <- do.call(
    cbind,
    lapply(
      seq_len(
        nrow(specifications)
      ),
      function(row_number) {
        
        dataset_name <-
          specifications$dataset[
            row_number
          ]
        
        expression_column <-
          specifications$expression_column[
            row_number
          ]
        
        current_table <- get_clean_table(
          dataset_name
        )
        
        model_rows <- match(
          common_genes,
          current_table$gene
        )
        
        as.numeric(
          current_table[
            model_rows,
            expression_column,
            drop = TRUE
          ]
        )
      }
    )
  )
  
  storage.mode(
    model_matrix
  ) <- "double"
  
  rownames(
    model_matrix
  ) <- common_genes
  
  colnames(
    model_matrix
  ) <- specifications$heatmap_label
  
  ##########################################################
  # Enforce complete values across every model and age
  ##########################################################
  
  combined_input_matrix <- cbind(
    reference_matrix,
    model_matrix
  )
  
  complete_gene_rows <- apply(
    combined_input_matrix,
    1,
    function(values) {
      all(
        is.finite(values)
      )
    }
  )
  
  reference_matrix <- reference_matrix[
    complete_gene_rows,
    ,
    drop = FALSE
  ]
  
  model_matrix <- model_matrix[
    complete_gene_rows,
    ,
    drop = FALSE
  ]
  
  common_genes <- common_genes[
    complete_gene_rows
  ]
  
  if (length(common_genes) < 2) {
    
    stop(
      cell_type,
      " has fewer than two complete shared genes."
    )
  }
  
  ##########################################################
  # Make sure no input column is completely constant
  ##########################################################
  
  constant_reference_columns <- colnames(
    reference_matrix
  )[
    apply(
      reference_matrix,
      2,
      function(values) {
        length(
          unique(values)
        ) < 2
      }
    )
  ]
  
  constant_model_columns <- colnames(
    model_matrix
  )[
    apply(
      model_matrix,
      2,
      function(values) {
        length(
          unique(values)
        ) < 2
      }
    )
  ]
  
  constant_columns <- c(
    constant_reference_columns,
    constant_model_columns
  )
  
  if (length(constant_columns) > 0) {
    
    stop(
      "Spearman correlation cannot be calculated because ",
      "these columns are constant across the common genes: ",
      paste(
        unique(constant_columns),
        collapse = ", "
      )
    )
  }
  
  ##########################################################
  # Models are rows; human ages are columns
  ##########################################################
  
  correlation_matrix <- stats::cor(
    model_matrix,
    reference_matrix,
    method = "spearman",
    use = "everything"
  )
  
  list(
    cell_type = cell_type,
    correlation = correlation_matrix,
    common_genes = common_genes,
    n_genes = length(common_genes)
  )
}

correlation_results <- setNames(
  lapply(
    cell_types,
    calculate_celltype_correlations
  ),
  cell_types
)

############################################################
# 7. Print common-gene counts
############################################################

common_gene_summary <- data.frame(
  
  cell_type = cell_types,
  
  common_genes_used = vapply(
    correlation_results,
    function(result) {
      
      getElement(
        result,
        "n_genes"
      )
    },
    numeric(1)
  ),
  
  same_genes_in_every_available_comparison =
    TRUE,
  
  stringsAsFactors = FALSE
)

cat(
  "\n========================================\n"
)

cat(
  "COMMON GENE INPUT SUMMARY\n"
)

cat(
  "========================================\n"
)

print(
  common_gene_summary,
  row.names = FALSE
)

############################################################
# 8. Print every correlation matrix
############################################################

for (cell_type in cell_types) {
  
  result <- getElement(
    correlation_results,
    cell_type
  )
  
  correlation_matrix <- getElement(
    result,
    "correlation"
  )
  
  n_genes <- getElement(
    result,
    "n_genes"
  )
  
  cat(
    "\n========================================\n"
  )
  
  cat(
    cell_type,
    ": ",
    format(
      n_genes,
      big.mark = ","
    ),
    " identical shared genes\n",
    sep = ""
  )
  
  cat(
    "========================================\n"
  )
  
  print(
    round(
      correlation_matrix,
      3
    )
  )
}









############################################################
# 9. Draw all four cell-type heatmaps
#
# - Every heatmap contains the same rows in the same order.
# - Missing cell types are displayed in gray.
# - White spaces separate different studies.
# - EN, IN, AST, and MG each have their own color scale.
# - The numerical color ranges are never shared.
############################################################

if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop(
    "Install gridExtra before running Step 9:\n",
    "install.packages('gridExtra')"
  )
}

############################################################
# 9A. Define common row order
#
# This controls only the row order.
# It does not control the color scale.
############################################################

common_row_information <- unique(
  model_specifications[
    ,
    c(
      "heatmap_label",
      "dataset"
    ),
    drop = FALSE
  ]
)

common_heatmap_rows <-
  common_row_information$heatmap_label

common_heatmap_studies <-
  common_row_information$dataset

############################################################
# 9B. Add white gaps between studies
############################################################

study_gap_positions <- which(
  common_heatmap_studies[
    seq_len(
      length(common_heatmap_studies) - 1
    )
  ] !=
    common_heatmap_studies[
      seq.int(
        2,
        length(common_heatmap_studies)
      )
    ]
)

############################################################
# 9C. Original blue-to-red palette
#
# The palette colors are the same, but the numerical mapping
# is calculated independently inside each cell-type plot.
############################################################




heatmap_colors <- grDevices::colorRampPalette(
  c(
    "#67a9cf", # Original Blue
    "#FFFFF0", # White
    "#a6d96a", # Light Green
    "#1a9641", # Medium Green
    "#00441b"  # Dark Green
    
  )
)(100)





############################################################
# 9D. Build one independently scaled cell-type heatmap
############################################################

if (!requireNamespace("gtable", quietly = TRUE)) {
  stop(
    "Install gtable before running Step 9:\n",
    "install.packages('gtable')"
  )
}

make_celltype_heatmap <- function(
    result
) {
  
  ##########################################################
  # Retrieve results for the current cell type
  ##########################################################
  
  cell_type <- getElement(
    result,
    "cell_type"
  )
  
  correlation_matrix <- getElement(
    result,
    "correlation"
  )
  
  ##########################################################
  # Validate row and column names
  ##########################################################
  
  if (is.null(rownames(correlation_matrix))) {
    stop(
      cell_type,
      " correlation matrix does not have row names."
    )
  }
  
  if (is.null(colnames(correlation_matrix))) {
    stop(
      cell_type,
      " correlation matrix does not have column names."
    )
  }
  
  if (anyDuplicated(rownames(correlation_matrix)) > 0) {
    stop(
      cell_type,
      " correlation matrix contains duplicated row names."
    )
  }
  
  unexpected_rows <- setdiff(
    rownames(correlation_matrix),
    common_heatmap_rows
  )
  
  if (length(unexpected_rows) > 0) {
    stop(
      cell_type,
      " contains unexpected model labels: ",
      paste(
        unexpected_rows,
        collapse = ", "
      )
    )
  }
  
  ##########################################################
  # Create identical row structure for all cell types
  #
  # Missing cell types remain NA and are displayed in gray.
  ##########################################################
  
  expanded_correlation_matrix <- matrix(
    NA_real_,
    nrow = length(common_heatmap_rows),
    ncol = ncol(correlation_matrix),
    dimnames = list(
      common_heatmap_rows,
      colnames(correlation_matrix)
    )
  )
  
  available_rows <- intersect(
    common_heatmap_rows,
    rownames(correlation_matrix)
  )
  
  expanded_correlation_matrix[
    available_rows,
    
  ] <- correlation_matrix[
    available_rows,
    ,
    drop = FALSE
  ]
  
  ##########################################################
  # Calculate the exact independent scale for this cell type
  #
  # No rounding to 0.05 is performed.
  # Each panel therefore uses its own true minimum and maximum.
  ##########################################################
  
  celltype_minimum <- min(
    correlation_matrix,
    na.rm = TRUE
  )
  
  celltype_maximum <- max(
    correlation_matrix,
    na.rm = TRUE
  )
  
  if (!is.finite(celltype_minimum)) {
    stop(
      cell_type,
      " does not contain a finite minimum correlation."
    )
  }
  
  if (!is.finite(celltype_maximum)) {
    stop(
      cell_type,
      " does not contain a finite maximum correlation."
    )
  }
  
  if (celltype_minimum == celltype_maximum) {
    celltype_minimum <-
      celltype_minimum - 0.01
    
    celltype_maximum <-
      celltype_maximum + 0.01
  }
  
  ##########################################################
  # Independent color breaks and legend for this cell type
  ##########################################################
  
  celltype_breaks <- seq(
    celltype_minimum,
    celltype_maximum,
    length.out =
      length(heatmap_colors) + 1
  )
  
  celltype_legend_breaks <- seq(
    celltype_minimum,
    celltype_maximum,
    length.out = 5
  )
  
  
  
  
  ##########################################################
  # Build the heatmap
  #
  # White cell borders create:
  # - thin vertical lines between developmental stages
  # - thin horizontal lines between datasets
  #
  # gaps_row still creates larger spaces between studies.
  ##########################################################
  
  heatmap_object <- pheatmap::pheatmap(
    expanded_correlation_matrix,
    
    color = heatmap_colors,
    breaks = celltype_breaks,
    
    na_col = "#BDBDBD",
    
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    
    gaps_row = rep(
      study_gap_positions,
      each = 3
    ),
    
    border_color = "white",
    
    display_numbers = FALSE,
    
    fontsize = 9,
    fontsize_row = 9,
    fontsize_col = 9,
    
    angle_col = 45,
    
    legend = TRUE,
    legend_breaks =
      celltype_legend_breaks,
    legend_labels = sprintf(
      "%.2f",
      celltype_legend_breaks
    ),
    
    main = cell_type,
    
    silent = TRUE
  )
  
  ##########################################################
  # Return the completed heatmap
  ##########################################################
  
  heatmap_gtable <- getElement(
    heatmap_object,
    "gtable"
  )
  
  heatmap_gtable
  
  
}




############################################################
# 9E. Build the four heatmaps independently
############################################################

EN_heatmap <- make_celltype_heatmap(
  getElement(
    correlation_results,
    "EN"
  )
)

IN_heatmap <- make_celltype_heatmap(
  getElement(
    correlation_results,
    "IN"
  )
)

AST_heatmap <- make_celltype_heatmap(
  getElement(
    correlation_results,
    "AST"
  )
)

MG_heatmap <- make_celltype_heatmap(
  getElement(
    correlation_results,
    "MG"
  )
)

############################################################
# 9F. Put the four independently scaled heatmaps together
#
# Each panel keeps its own separate color bar.
############################################################

grid::grid.newpage()

gridExtra::grid.arrange(
  EN_heatmap,
  IN_heatmap,
  AST_heatmap,
  MG_heatmap,
  ncol = 2,
  top = grid::textGrob(
    "Model versus human developmental-age Spearman correlation",
    gp = grid::gpar(
      fontsize = 16,
      fontface = "bold"
    )
  )
)








