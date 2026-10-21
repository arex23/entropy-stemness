SAMPLES_TO_RUN <- "sample1"
#SAMPLES_TO_RUN <- "all"

MIN_SPOT_PROP <- 0.10

SAMPLE_P      <- 0.25
DF_TPS        <- 6
ADJ_METHOD    <- "logpac"
RANDOM_SEED   <- 36
SAVE_COMPRESS <- FALSE 

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(SpaNorm)
  library(here)
})

dir_input_objects <- here::here("results", "01_load_and_compute_entropy", "objects")

dir_results <- here::here("results", "05_spanorm_normalization")
dir_objects <- file.path(dir_results, "objects")

dir.create(dir_objects, recursive = TRUE, showWarnings = FALSE)


all_rds_files <- list.files(dir_input_objects, pattern = "_seurat\\.rds$", full.names = FALSE)
all_samples   <- gsub("_seurat\\.rds$", "", all_rds_files)

sample_nums <- suppressWarnings(as.numeric(gsub("[^0-9]", "", all_samples)))
if (!any(is.na(sample_nums))) {
  all_samples <- all_samples[order(sample_nums)]
}

if (length(SAMPLES_TO_RUN) == 1 && SAMPLES_TO_RUN == "all") {
  samples_to_process <- all_samples
} else {
  missing_samples <- setdiff(SAMPLES_TO_RUN, all_samples)
  if (length(missing_samples) > 0) {
    stop("The following requested samples were not found in ", dir_input_objects, ": ",
         paste(missing_samples, collapse = ", "))
  }
  samples_to_process <- SAMPLES_TO_RUN
}

message("Samples queued for SpaNorm normalization (", length(samples_to_process), "): ",
        paste(samples_to_process, collapse = ", "))


for (sample_name in samples_to_process) {
  message("\n==================================================================")
  message("Processing ", sample_name, "...")
  message("==================================================================")
  
  rds_in_file <- file.path(dir_input_objects, paste0(sample_name, "_seurat.rds"))
  
  message("-> Loading Seurat object from: ", rds_in_file)
  seurat_obj <- readRDS(rds_in_file)
  
  counts_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  n_genes_total <- nrow(counts_mat)
  n_spots_total <- ncol(counts_mat)
  
  spot_prop  <- Matrix::rowSums(counts_mat > 0) / n_spots_total
  keep_genes <- names(spot_prop)[spot_prop >= MIN_SPOT_PROP]
  if (length(keep_genes) == 0) {
    stop("No genes passed the filter for ", sample_name)
  }
  
  seurat_obj <- subset(seurat_obj, features = keep_genes)
  counts_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  rm(spot_prop, keep_genes)
  
  spot_totals <- Matrix::colSums(counts_mat)
  if (any(spot_totals == 0)) {
    zero_spots <- colnames(counts_mat)[spot_totals == 0]
    message(
      "-> Removing ", length(zero_spots), " zero-count spot(s) in ", sample_name, ": ",
      paste(zero_spots, collapse = ", ")
    )
    keep_spots <- colnames(counts_mat)[spot_totals > 0]
    if (length(keep_spots) == 0) {
      stop("No spots with positive counts left in ", sample_name)
    }
    rm(counts_mat)
    seurat_obj <- subset(seurat_obj, cells = keep_spots)
    counts_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  }
  rm(spot_totals)
  
  n_genes_kept <- nrow(counts_mat)
  n_spots_kept <- ncol(counts_mat)
  message(sprintf(
    "-> Keeping %d / %d genes across %d / %d spots (min spot proportion >= %.2f)",
    n_genes_kept, n_genes_total, n_spots_kept, n_spots_total, MIN_SPOT_PROP
  ))
  
  coords_df <- Seurat::GetTissueCoordinates(seurat_obj)
  if ("cell" %in% colnames(coords_df)) {
    rownames(coords_df) <- coords_df$cell
  }
  coords_mat <- as.matrix(coords_df[colnames(counts_mat), c("x", "y"), drop = FALSE])
  
  spe <- SpatialExperiment::SpatialExperiment(
    assays        = list(counts = counts_mat),
    spatialCoords = coords_mat
  )
  rm(counts_mat, coords_df, coords_mat)
  
  message("-> Running SpaNorm (sample.p = ", SAMPLE_P,
          ", df.tps = ", DF_TPS,
          ", adj.method = '", ADJ_METHOD, "')...")
  set.seed(RANDOM_SEED)
  spe <- SpaNorm::SpaNorm(
    spe        = spe,
    sample.p   = SAMPLE_P,
    df.tps     = DF_TPS,
    adj.method = ADJ_METHOD,
    verbose    = TRUE
  )
  
  seurat_obj <- Seurat::SetAssayData(
    object   = seurat_obj,
    assay    = "Spatial",
    layer    = "data",
    new.data = SingleCellExperiment::logcounts(spe)
  )
  seurat_obj@misc$SpaNorm <- S4Vectors::metadata(spe)$SpaNorm
  
  rm(spe)
  gc(verbose = FALSE)
  
  rds_out_file <- file.path(dir_objects, paste0(sample_name, "_seurat.rds"))
  message("-> Saving SpaNorm-normalized Seurat object to: ", rds_out_file)
  saveRDS(seurat_obj, file = rds_out_file, compress = SAVE_COMPRESS)
  
  rm(seurat_obj)
  gc(verbose = FALSE)
  
  message("-> Done with ", sample_name, "!")
}

