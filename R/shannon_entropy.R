calculate_shannon_entropy <- function(seurat_obj, assay = NULL, layer = "counts", col_name = "shannon_entropy") {
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
  col_sums <- Matrix::colSums(expr_mat)
  
  if (inherits(expr_mat, "dgCMatrix")) {
    
    counts_per_cell <- diff(expr_mat@p)
    
    rep_col_sums <- rep(col_sums, counts_per_cell)
    
    p <- expr_mat@x / rep_col_sums
    
    entropy_components <- - (p * log2(p))
    
    col_indices <- rep(seq_len(ncol(expr_mat)), counts_per_cell)
    cell_entropies <- tapply(entropy_components, col_indices, sum)
    
    full_entropy <- rep(NA_real_, ncol(expr_mat))
    names(full_entropy) <- colnames(expr_mat)
    
    valid_indices <- as.numeric(names(cell_entropies))
    full_entropy[valid_indices] <- cell_entropies
    
  } else {
    full_entropy <- apply(expr_mat, 2, function(x) {
      x <- x[x > 0]
      if (length(x) == 0) return(NA_real_)
      p <- x / sum(x)
      -sum(p * log2(p))
    })
  }
  
  seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = full_entropy, col.name = col_name)
  
  return(seurat_obj)
}
