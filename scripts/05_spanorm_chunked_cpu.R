# Purpose: Normalize each sample with SpaNorm on the CPU, using gene chunks
#          for count adjustment, and record elapsed normalization time.
# Input:   results/01_load_and_compute_entropy/objects/sampleX_seurat.rds
# Output:  results/05_spanorm_chunked_cpu/objects/sampleX_seurat.rds
#          results/05_spanorm_chunked_cpu/tables/timings_<run>.csv
#          results/05_spanorm_chunked_cpu/tables/session_<run>.txt
# Run from the project root in the gpu-meningioma-entropy Conda environment.
# Run the CPU and GPU scripts separately, in fresh R sessions.

# Keep these settings identical in the CPU and GPU scripts.
#SAMPLES_TO_RUN <- "sample1"  # Or "all", or c("sample1", "sample2").
#SAMPLES_TO_RUN <- "all"
SAMPLES_TO_RUN <- paste0("sample", 8:22)
MIN_SPOT_PROP <- 0.10
SAMPLE_P <- 0.25
DF_TPS <- 6
ADJ_METHOD <- "logpac"
RANDOM_SEED <- 36
ADJUSTMENT_BLOCK_MIB <- 128
ADJUSTMENT_WORKERS <- 2L  # CPU workers for step (2/2) logpac adjustment; use 1L for serial.
SAVE_COMPRESS <- FALSE    # FALSE avoids ~4-5 min of gzip compression on the 1.1 GiB dense matrix.
BACKEND <- "cpu"

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(DelayedArray)
  library(BiocParallel)
  library(SpaNorm)
  library(here)
})

if (packageVersion("SpaNorm") < "1.7.15") {
  stop("Use the GitHub SpaNorm installation in gpu-meningioma-entropy.")
}

# A DelayedArray wrapper tells SpaNorm to adjust counts in gene blocks. The
# underlying counts remain sparse and in RAM; no HDF5 files are needed.
# This budget describes ONE dense matrix per adjustment block, not total RAM.
# SpaNorm still holds the final matrix and the sampled fitting data in RAM;
# its input checks can also temporarily realize the count matrix.
DelayedArray::setAutoBlockSize(ADJUSTMENT_BLOCK_MIB * 1024^2)
bpparam <- if (ADJUSTMENT_WORKERS > 1L) {
  BiocParallel::MulticoreParam(workers = ADJUSTMENT_WORKERS)
} else {
  BiocParallel::SerialParam()
}

dir_input <- here::here("results", "01_load_and_compute_entropy", "objects")
dir_results <- here::here("results", "05_spanorm_chunked_cpu")
dir_objects <- file.path(dir_results, "objects")
dir_tables <- file.path(dir_results, "tables")

all_files <- list.files(dir_input, pattern = "^sample[0-9]+_seurat\\.rds$")
all_samples <- sub("_seurat\\.rds$", "", all_files)
all_samples <- all_samples[order(as.integer(sub("sample", "", all_samples)))]
if (length(all_samples) == 0) stop("No input objects found in ", dir_input)

samples_to_process <- if (identical(SAMPLES_TO_RUN, "all")) all_samples else SAMPLES_TO_RUN
missing_samples <- setdiff(samples_to_process, all_samples)
if (length(missing_samples) > 0) {
  stop("Input objects missing for: ", paste(missing_samples, collapse = ", "))
}

dir.create(dir_objects, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
timing_file <- file.path(dir_tables, paste0("timings_", run_id, ".csv"))
capture.output(sessionInfo(), file = file.path(dir_tables, paste0("session_", run_id, ".txt")))
timings <- data.frame()

for (sample_name in samples_to_process) {
  message("\nProcessing ", sample_name, " with backend = ", BACKEND)
  input_file <- file.path(dir_input, paste0(sample_name, "_seurat.rds"))
  seurat_obj <- readRDS(input_file)
  counts_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  n_genes_before <- nrow(counts_mat)
  n_spots_before <- ncol(counts_mat)

  # Retain the same filter as script 05: genes detected in >=10% of spots.
  # This changes the gene set, so use the same filter in both comparisons.
  keep_genes <- rownames(counts_mat)[
    Matrix::rowSums(counts_mat > 0) / ncol(counts_mat) >= MIN_SPOT_PROP
  ]
  if (length(keep_genes) == 0) stop("No genes passed the filter for ", sample_name)
  rm(counts_mat)
  seurat_obj <- subset(seurat_obj, features = keep_genes)
  counts_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")

  # Remove spots with zero total counts after gene filtering. SpaNorm requires
  # positive library sizes, and zero-count spots have undefined Shannon entropy.
  spot_totals <- Matrix::colSums(counts_mat)
  if (any(spot_totals == 0)) {
    zero_spots <- colnames(counts_mat)[spot_totals == 0]
    message(
      "Removing ", length(zero_spots), " zero-count spot(s) in ", sample_name, ": ",
      paste(zero_spots, collapse = ", ")
    )
    keep_spots <- colnames(counts_mat)[spot_totals > 0]
    if (length(keep_spots) == 0) stop("No spots with positive counts left in ", sample_name)
    rm(counts_mat)
    seurat_obj <- subset(seurat_obj, cells = keep_spots)
    counts_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  }
  rm(spot_totals)
  n_genes <- nrow(counts_mat)
  n_spots <- ncol(counts_mat)

  # Align image coordinates to the count matrix by spot barcode.
  coords_df <- Seurat::GetTissueCoordinates(seurat_obj)
  if ("cell" %in% colnames(coords_df)) rownames(coords_df) <- coords_df$cell
  if (!all(c("x", "y") %in% colnames(coords_df)) ||
      !all(colnames(counts_mat) %in% rownames(coords_df))) {
    stop("Cannot align x/y coordinates to spots in ", sample_name)
  }
  coords_mat <- as.matrix(coords_df[colnames(counts_mat), c("x", "y"), drop = FALSE])
  if (any(!is.finite(coords_mat)) || any(apply(coords_mat, 2, function(x) diff(range(x))) == 0)) {
    stop("Coordinates are missing or have no spatial range in ", sample_name)
  }

  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = DelayedArray::DelayedArray(counts_mat)),
    spatialCoords = coords_mat
  )
  rm(counts_mat, coords_df, coords_mat, keep_genes)
  gc(verbose = FALSE)
  message(
    "Keeping ", n_genes, "/", n_genes_before, " genes across ",
    n_spots, "/", n_spots_before, " spots"
  )

  # Fit ONE spatial model for the complete sample. SpaNorm then applies that
  # model in gene chunks, preserving the shared library-size and dispersion
  # calculations. Splitting samples into independent SpaNorm calls would
  # change those calculations and would not be the same analysis.
  set.seed(RANDOM_SEED)
  elapsed_seconds <- system.time({
    spe <- SpaNorm::SpaNorm(
      spe = spe,
      sample.p = SAMPLE_P,
      df.tps = DF_TPS,
      adj.method = ADJ_METHOD,
      backend = BACKEND,
      BPPARAM = bpparam,
      overwrite = TRUE,
      verbose = TRUE
    )
  })[["elapsed"]]

  # logpac is log2(percentile-adjusted counts + 1), not raw counts or entropy.
  seurat_obj <- Seurat::SetAssayData(
    seurat_obj, assay = "Spatial", layer = "data",
    new.data = SingleCellExperiment::logcounts(spe)
  )
  seurat_obj@misc$SpaNorm <- S4Vectors::metadata(spe)$SpaNorm
  rm(spe)
  gc(verbose = FALSE)

  # Timing includes fitting and count adjustment; loading, filtering, and
  # saving are excluded. A new CSV per run preserves repeated measurements.
  timings <- rbind(timings, data.frame(
    sample = sample_name, backend = BACKEND, elapsed_seconds = elapsed_seconds,
    genes_before = n_genes_before, genes_used = n_genes,
    spots_before = n_spots_before, spots = n_spots,
    min_spot_prop = MIN_SPOT_PROP, sample_p = SAMPLE_P, df_tps = DF_TPS,
    adj_method = ADJ_METHOD, seed = RANDOM_SEED,
    adjustment_block_mib = ADJUSTMENT_BLOCK_MIB,
    adjustment_workers = ADJUSTMENT_WORKERS,
    spanorm_version = as.character(packageVersion("SpaNorm")),
    openblas_threads = Sys.getenv("OPENBLAS_NUM_THREADS", unset = "default")
  ))
  write.csv(timings, timing_file, row.names = FALSE)
  # Repeated runs replace this backend's saved object, but keep their timings.
  saveRDS(
    seurat_obj,
    file.path(dir_objects, paste0(sample_name, "_seurat.rds")),
    compress = SAVE_COMPRESS
  )
  rm(seurat_obj)
  gc(verbose = FALSE)
  message(sample_name, ": ", round(elapsed_seconds, 1), " seconds for SpaNorm")
}

message("Timing results: ", timing_file)
