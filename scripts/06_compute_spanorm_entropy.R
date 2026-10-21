#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

SAVE_OBJECTS <- FALSE

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(here)
})

source(here::here("R", "shannon_entropy.R"))

dir_input_objects <- here::here("results", "05_spanorm_normalization", "objects")

dir_results <- here::here("results", "06_compute_spanorm_entropy")
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")
dir_objects <- file.path(dir_results, "objects")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)
if (isTRUE(SAVE_OBJECTS)) {
  dir.create(dir_objects, recursive = TRUE, showWarnings = FALSE)
}


all_rds_files <- list.files(dir_input_objects, pattern = "_seurat\\.rds$", full.names = FALSE)
all_samples   <- gsub("_seurat\\.rds$", "", all_rds_files)

sample_numbers <- suppressWarnings(as.numeric(gsub("[^0-9]", "", all_samples)))
if (!any(is.na(sample_numbers))) {
  all_samples <- all_samples[order(sample_numbers)]
}

if (identical(SAMPLES_TO_RUN, "all")) {
  samples_to_process <- all_samples
} else {
  missing_samples <- setdiff(SAMPLES_TO_RUN, all_samples)
  if (length(missing_samples) > 0) {
    stop("The following requested samples were not found in ", dir_input_objects, ": ",
         paste(missing_samples, collapse = ", "))
  }
  samples_to_process <- SAMPLES_TO_RUN
}

message("=== Starting Analysis: 06_compute_spanorm_entropy ===")
message("Samples selected to process (", length(samples_to_process), "): ",
        paste(samples_to_process, collapse = ", "))


for (sample_name in samples_to_process) {
  message("\n------------------------------------------------------------")
  message("Processing: ", sample_name)
  message("------------------------------------------------------------")
  
  rds_file <- file.path(dir_input_objects, paste0(sample_name, "_seurat.rds"))
  if (!file.exists(rds_file)) {
    warning("SpaNorm Seurat object not found for ", sample_name, ": ", rds_file, ". Skipping...")
    next
  }
  
  message("-> Loading SpaNorm Seurat object from: ", basename(rds_file))
  seurat_obj <- readRDS(rds_file)
  
  message("-> Extracting SpaNorm normalized values from 'data' layer...")
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    log_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "data")
  } else {
    log_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", slot = "data")
  }
  
  message("-> Reversing log2 transformation (2^x - 1)...")
  unlogged_mat <- 2^log_mat - 1
  unlogged_mat[unlogged_mat < 0] <- 0
  
  unlogged_sparse <- Matrix::drop0(Matrix::Matrix(unlogged_mat, sparse = TRUE))
  rm(log_mat, unlogged_mat)
  gc(verbose = FALSE)
  
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    seurat_obj <- Seurat::SetAssayData(
      object   = seurat_obj,
      assay    = "Spatial",
      layer    = "unlogged_spanorm",
      new.data = unlogged_sparse
    )
  } else {
    seurat_obj <- Seurat::SetAssayData(
      object   = seurat_obj,
      assay    = "Spatial",
      slot     = "unlogged_spanorm",
      new.data = unlogged_sparse
    )
  }
  
  message("-> Computing Shannon entropy on unlogged SpaNorm values...")
  seurat_obj <- calculate_shannon_entropy(
    seurat_obj = seurat_obj,
    assay      = "Spatial",
    layer      = "unlogged_spanorm",
    col_name   = "spanorm_entropy"
  )
  
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    tryCatch({
      seurat_obj[["Spatial"]]$unlogged_spanorm <- NULL
    }, error = function(e) NULL)
  }
  rm(unlogged_sparse)
  gc(verbose = FALSE)
  
  csv_file <- file.path(dir_tables, paste0(sample_name, "_spanorm_entropy.csv"))
  message("-> Saving summary table to: ", csv_file)
  
  entropy_df <- data.frame(
    spot_id         = colnames(seurat_obj),
    spanorm_entropy = seurat_obj$spanorm_entropy,
    entropy         = seurat_obj$spanorm_entropy,  
    nCount_RNA      = seurat_obj$nCount_Spatial,
    nFeature_RNA    = seurat_obj$nFeature_Spatial,
    sample          = sample_name,
    stringsAsFactors = FALSE
  )
  write.csv(entropy_df, file = csv_file, row.names = FALSE)
  
  fig_file <- file.path(dir_figures, paste0(sample_name, "_spatial_spanorm_entropy.png"))
  message("-> Generating dual spatial plot: ", fig_file)
  
  p1 <- Seurat::SpatialPlot(seurat_obj, pt.size.factor = 0) +
    ggplot2::ggtitle(paste0(sample_name, " - H&E Histology")) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "none"
    )
  
  p2 <- Seurat::SpatialFeaturePlot(seurat_obj, features = "spanorm_entropy") +
    ggplot2::ggtitle(paste0(sample_name, " - SpaNorm Shannon Entropy")) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  dual_plot <- p1 | p2
  
  ggplot2::ggsave(
    filename = fig_file,
    plot     = dual_plot,
    width    = 12,
    height   = 6,
    dpi      = 300
  )
  
  if (isTRUE(SAVE_OBJECTS)) {
    out_rds <- file.path(dir_objects, paste0(sample_name, "_seurat.rds"))
    message("-> Saving updated Seurat object to: ", out_rds)
    saveRDS(seurat_obj, file = out_rds, compress = FALSE)
  }
  
  rm(seurat_obj, entropy_df, dual_plot, p1, p2)
  gc(verbose = FALSE)
  
  message("-> Done with ", sample_name, "!")
}

