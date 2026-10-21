calculate_chao_shen_entropy <- function(seurat_obj, assay = NULL, layer = "counts", col_name = "chao_shen_entropy") {
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
    col_indices     <- rep(seq_len(ncol(expr_mat)), counts_per_cell)
    
    is_singleton <- as.numeric(expr_mat@x == 1)
    f1_per_col   <- numeric(ncol(expr_mat))
    f1_tab       <- tapply(is_singleton, col_indices, sum)
    f1_per_col[as.numeric(names(f1_tab))] <- f1_tab
    
    f1_adj   <- ifelse(f1_per_col == col_sums & col_sums > 0, f1_per_col - 1, f1_per_col)
    coverage <- ifelse(col_sums > 0, 1 - (f1_adj / col_sums), NA_real_)
    
    rep_col_sums <- rep(col_sums, counts_per_cell)
    rep_coverage <- rep(coverage, counts_per_cell)
    p_adj        <- (expr_mat@x / rep_col_sums) * rep_coverage
    
    denom <- 1 - (1 - p_adj)^rep_col_sums
    entropy_components <- ifelse(
      p_adj > 0 & denom > 0,
      -(p_adj * log2(p_adj)) / denom,
      0
    )
    
    cell_entropies <- tapply(entropy_components, col_indices, sum)
    
    full_entropy <- rep(NA_real_, ncol(expr_mat))
    names(full_entropy) <- colnames(expr_mat)
    
    valid_indices <- as.numeric(names(cell_entropies))
    full_entropy[valid_indices] <- cell_entropies
    
  } else {
    full_entropy <- apply(expr_mat, 2, function(x) {
      x <- x[x > 0]
      if (length(x) == 0) return(NA_real_)
      n  <- sum(x)
      f1 <- sum(x == 1)
      if (f1 == n) f1 <- n - 1
      C_hat <- 1 - (f1 / n)
      p_adj <- C_hat * (x / n)
      denom <- 1 - (1 - p_adj)^n
      comp  <- ifelse(p_adj > 0 & denom > 0, -(p_adj * log2(p_adj)) / denom, 0)
      sum(comp)
    })
  }
  
  seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = full_entropy, col.name = col_name)
  
  return(seurat_obj)
}
