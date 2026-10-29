# Purpose: Calculate the 2024 paper meningioma stemness module score using Seurat::AddModuleScore,
#          and generate a spatial feature plot and summary table for each sample.
# Input:   Seurat objects from results/00_load_and_filter_seurat/objects/filtered/<sample>_seurat.rds
# Output:  results/10_compute_module_score/tables/<sample>_module_score.csv
#          results/10_compute_module_score/figures/<sample>_spatial_module_score.png

# 1. Configuration and dependencies

#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

# Stemness signature from 2024 paper
MODULE_NAME  <- "module_score"
MODULE_GENES <- c("NOTCH3", "THY1", "CXCR4", "PROM1", "CDH2", "NES")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(here)
})

dir_input_objects <- here::here("results", "00_load_and_filter_seurat", "objects", "filtered")

dir_results <- here::here("results", "10_compute_module_score")
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

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

message("=== Starting Analysis: 10_compute_module_score ===")
message("Samples selected to process (", length(samples_to_process), "): ",
        paste(samples_to_process, collapse = ", "))

# 2. Processing loop

for (sample_name in samples_to_process) {
  message("Processing: ", sample_name)
  
  # A. Load Seurat object
  rds_file <- file.path(dir_input_objects, paste0(sample_name, "_seurat.rds"))
  if (!file.exists(rds_file)) {
    warning("Seurat .rds file not found for ", sample_name, ": ", rds_file, ". Skipping...")
    next
  }
  
  message("-> Loading Seurat object from: ", basename(rds_file))
  seurat_obj <- readRDS(rds_file)
  
  # B. Normalize data
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  
  # C. Filter module features present in the dataset
  genes_present <- intersect(MODULE_GENES, rownames(seurat_obj))
  
  if (length(genes_present) == 0) {
    warning("None of the module genes are present in ", sample_name, ". Skipping...")
    next
  }
  
  message("-> Module genes present (", length(genes_present), "/", length(MODULE_GENES), "): ",
          paste(genes_present, collapse = ", "))
  
  # D. Compute module score
  # AddModuleScore splits all genes into 24 expression bins and samples 'ctrl'
  # control genes per bin (default 100). If a filtered sample has fewer genes,
  # a bin may have fewer than 100 genes. We cap 'ctrl' based on total gene count.
  ctrl_genes <- min(100, floor(nrow(seurat_obj) / 30))
  message("-> Using ctrl = ", ctrl_genes, " control genes per feature (", nrow(seurat_obj), " total genes in object)")
  
  seurat_obj <- Seurat::AddModuleScore(
    object   = seurat_obj,
    features = list(genes_present),
    name     = MODULE_NAME,
    ctrl     = ctrl_genes,
    seed     = 1
  )
  seurat_obj[[MODULE_NAME]] <- seurat_obj[[paste0(MODULE_NAME, "1"), drop = TRUE]]
  
  # E. Save summary table
  csv_file <- file.path(dir_tables, paste0(sample_name, "_module_score.csv"))
  message("-> Saving summary table to: ", basename(csv_file))
  
  score_df <- data.frame(
    spot_id      = colnames(seurat_obj),
    module_score = seurat_obj[[MODULE_NAME, drop = TRUE]],
    nCount_RNA   = seurat_obj$nCount_Spatial,
    nFeature_RNA = seurat_obj$nFeature_Spatial,
    sample       = sample_name,
    stringsAsFactors = FALSE
  )
  write.csv(score_df, file = csv_file, row.names = FALSE)
  
  # F. Save spatial plot
  fig_file <- file.path(dir_figures, paste0(sample_name, "_spatial_module_score.png"))
  message("-> Generating spatial plot: ", basename(fig_file))
  
  p1 <- Seurat::SpatialPlot(seurat_obj, pt.size.factor = 0) +
    ggplot2::ggtitle(paste0(sample_name, " - H&E")) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "none"
    )
  
  p2 <- Seurat::SpatialFeaturePlot(seurat_obj, features = MODULE_NAME) +
    ggplot2::ggtitle(paste0(sample_name, " - 2024 Paper Signature")) +
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
  
  # G. Clean up memory
  rm(seurat_obj, score_df, dual_plot, p1, p2)
  gc(verbose = FALSE)
  
  message("-> Done with ", sample_name, "!")
}

message("\n=== Analysis complete for script 10 ===")
