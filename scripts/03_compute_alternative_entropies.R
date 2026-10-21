#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(here)
})

source(here::here("R", "chao_shen_entropy.R"))
source(here::here("R", "chao_wang_jost_entropy.R"))

dir_input_objects <- here::here("results", "01_load_and_compute_entropy", "objects")

dir_results <- here::here("results", "03_compute_alternative_entropies")
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

message("=== Starting Analysis: 03_compute_alternative_entropies ===")
message("Samples selected: ", paste(samples_to_process, collapse = ", "))


for (sample_name in samples_to_process) {
  message("\n------------------------------------------------------------")
  message("Processing: ", sample_name)
  message("------------------------------------------------------------")
  

  rds_file <- file.path(dir_input_objects, paste0(sample_name, "_seurat.rds"))
  if (!file.exists(rds_file)) {
    warning("Seurat .rds file not found for ", sample_name, ": ", rds_file, ". Skipping...")
    next
  }
  
  message("-> Loading Seurat object from: ", basename(rds_file))
  seurat_obj <- readRDS(rds_file)
  
  message("-> Calculating Chao-Shen entropy per spot...")
  seurat_obj <- calculate_chao_shen_entropy(
    seurat_obj = seurat_obj,
    assay      = "Spatial",
    layer      = "counts",
    col_name   = "chao_shen_entropy"
  )
  
  message("-> Calculating Chao-Wang-Jost entropy per spot...")
  seurat_obj <- calculate_cwj_entropy(
    seurat_obj = seurat_obj,
    assay      = "Spatial",
    layer      = "counts",
    col_name   = "cwj_entropy"
  )
  
  csv_file <- file.path(dir_tables, paste0(sample_name, "_alternative_entropies.csv"))
  message("-> Saving summary table to: ", csv_file)
  
  entropy_df <- data.frame(
    spot_id           = colnames(seurat_obj),
    chao_shen_entropy = seurat_obj$chao_shen_entropy,
    cwj_entropy       = seurat_obj$cwj_entropy,
    nCount_RNA        = seurat_obj$nCount_Spatial,
    nFeature_RNA      = seurat_obj$nFeature_Spatial,
    sample            = sample_name,
    stringsAsFactors  = FALSE
  )
  write.csv(entropy_df, file = csv_file, row.names = FALSE)
  
  fig_file <- file.path(dir_figures, paste0(sample_name, "_spatial_alternative_entropies.png"))
  message("-> Generating spatial plot: ", fig_file)
  
  p1 <- Seurat::SpatialPlot(seurat_obj, pt.size.factor = 0) +
    ggplot2::ggtitle(paste0(sample_name, " - H&E Histology")) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "none"
    )
  
  p2 <- Seurat::SpatialFeaturePlot(seurat_obj, features = "chao_shen_entropy") +
    ggplot2::ggtitle(paste0(sample_name, " - Chao-Shen Entropy")) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  p3 <- Seurat::SpatialFeaturePlot(seurat_obj, features = "cwj_entropy") +
    ggplot2::ggtitle(paste0(sample_name, " - Chao-Wang-Jost Entropy")) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  combined_plot <- p1 | p2 | p3
  
  ggplot2::ggsave(
    filename = fig_file,
    plot     = combined_plot,
    width    = 18,
    height   = 6,
    dpi      = 300
  )
  
  rm(seurat_obj, entropy_df, combined_plot, p1, p2, p3)
  gc(verbose = FALSE)
  
  message("-> Done with ", sample_name, "!")
}
