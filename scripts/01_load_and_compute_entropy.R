# Purpose: Load the 10%-filtered Seurat objects from Script 00, calculate
#          Plug-in Shannon entropy per spot, and save summary tables and plots.
# Input:   results/00_load_and_filter_seurat/objects/filtered/<sample>_seurat.rds
# Output:  results/01_load_and_compute_entropy/tables/<sample>_entropy.csv
#          results/01_load_and_compute_entropy/figures/<sample>_spatial_entropy.png

# 1. Dependencies and configuration

#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

library(Seurat)
library(ggplot2)
library(patchwork)
library(here)

source(here::here("R", "shannon_entropy.R"))

dir_input_objects <- here::here("results", "00_load_and_filter_seurat", "objects", "filtered")

dir_results <- here::here("results", "01_load_and_compute_entropy")
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

all_rds_files <- list.files(dir_input_objects, pattern = "_seurat\\.rds$", full.names = FALSE)
all_samples   <- gsub("_seurat\\.rds$", "", all_rds_files)

sample_numbers <- as.numeric(gsub("sample", "", all_samples))
all_samples    <- all_samples[order(sample_numbers)]

if (identical(SAMPLES_TO_RUN, "all")) {
  samples_to_process <- all_samples
} else {
  samples_to_process <- SAMPLES_TO_RUN
}

# 2. Processing loop

for (sample_name in samples_to_process) {
  
  # A. Load filtered Seurat object
  message("Processing: ", sample_name)
  
  rds_file <- file.path(dir_input_objects, paste0(sample_name, "_seurat.rds"))
  if (!file.exists(rds_file)) {
    warning("Filtered Seurat .rds file not found for ", sample_name, ": ", rds_file, ". Skipping...")
    next
  }
  
  message("Loading filtered Seurat object from: ", basename(rds_file))
  seurat_obj <- readRDS(rds_file)
  
  # B. Calculate Shannon entropy
  
  message("Calculating Shannon entropy")
  seurat_obj <- calculate_shannon_entropy(
    seurat_obj = seurat_obj,
    assay      = "Spatial",
    layer      = "counts",
    col_name   = "shannon_entropy"
  )

  # C. Save summary table 
  
  csv_file <- file.path(dir_tables, paste0(sample_name, "_entropy.csv"))
  message("Saving summary table to: ", csv_file)
  
  entropy_df <- data.frame(
    spot_id      = colnames(seurat_obj),
    entropy      = seurat_obj$shannon_entropy,
    nCount_RNA   = seurat_obj$nCount_Spatial,
    nFeature_RNA = seurat_obj$nFeature_Spatial,
    sample       = sample_name,
    stringsAsFactors = FALSE
  )
  write.csv(entropy_df, file = csv_file, row.names = FALSE)
  
  # D. Save spatial plots
  
  fig_file <- file.path(dir_figures, paste0(sample_name, "_spatial_entropy.png"))
  message("Generating dual spatial plot: ", fig_file)
  
  p1 <- Seurat::SpatialPlot(seurat_obj, pt.size.factor = 0) +
    ggplot2::ggtitle(paste0(sample_name, " - H&E Histology")) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "none"
    )
  
  p2 <- Seurat::SpatialFeaturePlot(seurat_obj, features = "shannon_entropy") +
    ggplot2::ggtitle(paste0(sample_name, " - Shannon Entropy")) +
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
  
  # E. Clean up memory
  
  rm(seurat_obj, entropy_df, dual_plot, p1, p2)
  gc(verbose = FALSE)
  
  message("Done with ", sample_name)
}
