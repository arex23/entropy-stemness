cwj_entropy_from_counts <- function(counts) {
  if (!is.numeric(counts) || any(!is.finite(counts)) ||
      any(counts < 0) || any(counts != floor(counts))) {
    stop("CWJ entropy requires finite, nonnegative integer counts.")
  }

  counts <- counts[counts > 0]
  if (length(counts) == 0) return(NA_real_)

  n <- sum(counts)
  if (!is.finite(n)) stop("The total count must be finite.")
  if (n == 1) return(0)

  observed_nats <- sum((counts / n) * (digamma(n) - digamma(counts)))

  f1 <- sum(counts == 1)
  f2 <- sum(counts == 2)
  if (f1 == 0) return(observed_nats / log(2))

  if (f2 > 0) {
    A <- (2 * f2) / ((n - 1) * f1 + 2 * f2)
  } else {
    A <- 2 / ((n - 1) * (f1 - 1) + 2)
  }
  if (A == 1) return(observed_nats / log(2))
  if (!is.finite(A) || A <= 0 || A > 1) {
    stop("The CWJ discovery parameter could not be calculated reliably.")
  }

  b <- 1 - A
  tail_integral <- stats::integrate(
    function(u) exp(-u) / (n * A + b * n * (-expm1(-u / n))),
    lower = 0,
    upper = Inf,
    rel.tol = 1e-8,
    abs.tol = 0,
    subdivisions = 200L,
    stop.on.error = TRUE
  )
  unseen_nats <- (f1 / n) * b * tail_integral$value

  (observed_nats + unseen_nats) / log(2)
}

calculate_cwj_entropy <- function(seurat_obj, assay = NULL, layer = "counts", col_name = "cwj_entropy") {
  # Ensure necessary packages are available
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The Seurat package is required but not installed.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required but not installed.")
  }
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }
  
  expr_mat <- tryCatch({
    if (utils::packageVersion("SeuratObject") >= "5.0.0") {
      Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
    } else {
      Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer)
    }
  }, error = function(e) {
    Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer)
  })
  
  if (is.null(expr_mat) || nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    stop(sprintf("No data found in assay '%s' and layer/slot '%s'.", assay, layer))
  }
  
  full_entropy <- rep(NA_real_, ncol(expr_mat))
  names(full_entropy) <- colnames(expr_mat)

  for (j in seq_len(ncol(expr_mat))) {
    if (inherits(expr_mat, "dgCMatrix")) {
      first <- expr_mat@p[j] + 1L
      last <- expr_mat@p[j + 1L]
      counts <- if (first <= last) expr_mat@x[seq.int(first, last)] else numeric(0)
    } else {
      counts <- as.numeric(expr_mat[, j])
    }
    full_entropy[j] <- cwj_entropy_from_counts(counts)
  }
  
  seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = full_entropy, col.name = col_name)
  
  return(seurat_obj)
}
