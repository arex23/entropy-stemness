#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(here)
})

source(here::here("R", "shannon_entropy.R"))

dir_raw <- here::here("data", "raw")

dir_results <- here::here("results", "01_load_and_compute_entropy")
dir_objects <- file.path(dir_results, "objects")
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")

dir.create(dir_objects, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

all_sample_dirs <- list.dirs(dir_raw, full.names = FALSE, recursive = FALSE)
all_sample_dirs <- grep("^sample", all_sample_dirs, value = TRUE)

sample_numbers  <- as.numeric(gsub("sample", "", all_sample_dirs))
all_sample_dirs <- all_sample_dirs[order(sample_numbers)]

if (identical(SAMPLES_TO_RUN, "all")) {
  samples_to_process <- all_sample_dirs
} else {
  samples_to_process <- SAMPLES_TO_RUN
}

message("=== Starting Analysis: 01_load_and_compute_entropy ===")
message("Samples selected to process: ", paste(samples_to_process, collapse = ", "))

for (sample_name in samples_to_process) {
  message("\n------------------------------------------------------------")
  message("Processing: ", sample_name)
  message("------------------------------------------------------------")
  
  sample_dir <- file.path(dir_raw, sample_name)
  
  if (!dir.exists(sample_dir)) {
    warning("Directory not found for ", sample_name, ". Skipping...")
    next
  }
  
  h5_files <- list.files(sample_dir, pattern = "\\.h5$", full.names = FALSE)
  if (length(h5_files) == 0) {
    warning("No .h5 matrix file found in ", sample_dir, ". Skipping...")
    next
  }
  h5_filename <- h5_files[1]
  
  message("-> Loading 10x Visium data from: ", h5_filename)
  seurat_obj <- Seurat::Load10X_Spatial(
    data.dir      = sample_dir,
    filename      = h5_filename,
    filter.matrix = TRUE  
  )

  positions <- read.csv(file.path(sample_dir, "spatial", "tissue_positions.csv"))
  in_tissue <- positions$in_tissue[match(colnames(seurat_obj), positions$barcode)]
  if (anyNA(in_tissue)) {
    stop("Missing tissue flags for spots in ", sample_name)
  }
  tissue_spots <- colnames(seurat_obj)[in_tissue == 1]
  message("-> Keeping ", length(tissue_spots), " of ", ncol(seurat_obj), " spots marked in tissue")
  seurat_obj <- subset(seurat_obj, cells = tissue_spots)
  
  seurat_obj$sample <- sample_name
  
  rds_file <- file.path(dir_objects, paste0(sample_name, "_seurat.rds"))
  message("-> Saving plain Seurat object to: ", rds_file)
  saveRDS(seurat_obj, file = rds_file)
  
  message("-> Calculating Shannon entropy per spot...")
  seurat_obj <- calculate_shannon_entropy(
    seurat_obj = seurat_obj,
    assay      = "Spatial",
    layer      = "counts",
    col_name   = "shannon_entropy"
  )
  
  csv_file <- file.path(dir_tables, paste0(sample_name, "_entropy.csv"))
  message("-> Saving summary table to: ", csv_file)
  
  entropy_df <- data.frame(
    spot_id      = colnames(seurat_obj),
    entropy      = seurat_obj$shannon_entropy,
    nCount_RNA   = seurat_obj$nCount_Spatial,
    nFeature_RNA = seurat_obj$nFeature_Spatial,
    sample       = sample_name,
    stringsAsFactors = FALSE
  )
  write.csv(entropy_df, file = csv_file, row.names = FALSE)
  
  fig_file <- file.path(dir_figures, paste0(sample_name, "_spatial_entropy.png"))
  message("-> Generating dual spatial plot: ", fig_file)
  
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
  
  rm(seurat_obj, entropy_df, dual_plot, p1, p2)
  gc(verbose = FALSE)
  
  message("-> Done with ", sample_name, "!")
}


