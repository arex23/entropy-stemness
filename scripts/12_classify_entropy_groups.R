#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

GROUP_PERCENTAGE <- 15

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(here)
})

if (!is.numeric(GROUP_PERCENTAGE) || length(GROUP_PERCENTAGE) != 1 ||
    GROUP_PERCENTAGE <= 0 || GROUP_PERCENTAGE >= 50) {
  stop("GROUP_PERCENTAGE must be a single number strictly between 0 and 50 (e.g., 20).")
}

pct_prop   <- GROUP_PERCENTAGE / 100
pct_folder <- paste0("top_bottom_", GROUP_PERCENTAGE, "pct")

dir_input_objects <- here::here("results", "05_spanorm_normalization", "objects")
dir_input_entropy <- here::here("results", "06_compute_spanorm_entropy", "tables")

dir_results <- here::here("results", "12_classify_entropy_groups", pct_folder)
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

all_csv_files <- list.files(dir_input_entropy, pattern = "_spanorm_entropy\\.csv$", full.names = FALSE)
all_samples   <- gsub("_spanorm_entropy\\.csv$", "", all_csv_files)

sample_numbers <- suppressWarnings(as.numeric(gsub("[^0-9]", "", all_samples)))
if (!any(is.na(sample_numbers))) {
  all_samples <- all_samples[order(sample_numbers)]
}

if (identical(SAMPLES_TO_RUN, "all")) {
  samples_to_process <- all_samples
} else {
  missing_samples <- setdiff(SAMPLES_TO_RUN, all_samples)
  if (length(missing_samples) > 0) {
    stop("The following requested samples were not found in ", dir_input_entropy, ": ",
         paste(missing_samples, collapse = ", "))
  }
  samples_to_process <- SAMPLES_TO_RUN
}

message("=== Starting Analysis: 12_classify_entropy_groups ===")
message("Quantile cutoff: Top/Bottom ", GROUP_PERCENTAGE, "% (output folder: ", pct_folder, ")")
message("Samples selected (", length(samples_to_process), "): ",
        paste(samples_to_process, collapse = ", "))

group_colors <- c(
  "High"   = "#d7191c",
  "Medium" = "#d9d9d9",
  "Low"    = "#2c7bb6"
)

summary_list <- list()

for (sample_name in samples_to_process) {
  message("\nProcessing: ", sample_name)
  
  rds_file <- file.path(dir_input_objects, paste0(sample_name, "_seurat.rds"))
  csv_file <- file.path(dir_input_entropy, paste0(sample_name, "_spanorm_entropy.csv"))
  
  if (!file.exists(rds_file) || !file.exists(csv_file)) {
    warning("Missing input Seurat object or SpaNorm entropy table for ", sample_name, ". Skipping...")
    next
  }
  
  seurat_obj <- readRDS(rds_file)
  entropy_df <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  ent_col <- if ("spanorm_entropy" %in% colnames(entropy_df)) "spanorm_entropy" else "entropy"
  
  spot_ids <- colnames(seurat_obj)
  matched_entropy <- entropy_df[[ent_col]][match(spot_ids, entropy_df$spot_id)]
  
  if (all(is.na(matched_entropy))) {
    warning("Could not match spot IDs for ", sample_name, ". Skipping...")
    rm(seurat_obj, entropy_df)
    next
  }
  
  seurat_obj$spanorm_entropy <- matched_entropy
  
  valid_ent <- matched_entropy[is.finite(matched_entropy)]
  low_cutoff  <- unname(quantile(valid_ent, probs = pct_prop,     na.rm = TRUE))
  high_cutoff <- unname(quantile(valid_ent, probs = 1 - pct_prop, na.rm = TRUE))
  
  entropy_group <- rep(NA_character_, length(matched_entropy))
  entropy_group[is.finite(matched_entropy)] <- "Medium"
  entropy_group[!is.na(matched_entropy) & matched_entropy <= low_cutoff]  <- "Low"
  entropy_group[!is.na(matched_entropy) & matched_entropy >= high_cutoff] <- "High"
  
  seurat_obj$entropy_group <- factor(entropy_group, levels = c("High", "Medium", "Low"))
  
  n_high   <- sum(entropy_group == "High",   na.rm = TRUE)
  n_medium <- sum(entropy_group == "Medium", na.rm = TRUE)
  n_low    <- sum(entropy_group == "Low",    na.rm = TRUE)
  
  message(sprintf(
    "-> Cutoffs: Low <= %.3f (%d spots) | High >= %.3f (%d spots) | Medium: %d spots",
    low_cutoff, n_low, high_cutoff, n_high, n_medium
  ))
  
  out_csv <- file.path(dir_tables, paste0(sample_name, "_entropy_groups.csv"))
  groups_df <- data.frame(
    spot_id         = spot_ids,
    spanorm_entropy = matched_entropy,
    entropy_group   = entropy_group,
    nCount_RNA      = seurat_obj$nCount_Spatial,
    nFeature_RNA    = seurat_obj$nFeature_Spatial,
    sample          = sample_name,
    stringsAsFactors = FALSE
  )
  write.csv(groups_df, file = out_csv, row.names = FALSE)
  message("-> Saved group assignments to: ", basename(out_csv))
  
  summary_list[[sample_name]] <- data.frame(
    sample           = sample_name,
    group_percentage = GROUP_PERCENTAGE,
    n_spots_total    = length(valid_ent),
    n_high           = n_high,
    n_medium         = n_medium,
    n_low            = n_low,
    low_cutoff       = low_cutoff,
    high_cutoff      = high_cutoff,
    mean_high_ent    = mean(matched_entropy[entropy_group == "High"],   na.rm = TRUE),
    mean_low_ent     = mean(matched_entropy[entropy_group == "Low"],    na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  fig_file <- file.path(dir_figures, paste0(sample_name, "_spatial_entropy_groups.png"))
  message("-> Generating spatial plot: ", basename(fig_file))
  
  p1 <- Seurat::SpatialPlot(seurat_obj, pt.size.factor = 0) +
    ggplot2::ggtitle(paste0(sample_name, " - H&E Histology")) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
      legend.position = "none"
    )
  
  p2 <- Seurat::SpatialFeaturePlot(seurat_obj, features = "spanorm_entropy") +
    ggplot2::ggtitle(paste0(sample_name, "SpaNorm Entropy")) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13)
    )
  
  p3 <- Seurat::SpatialDimPlot(
    seurat_obj,
    group.by = "entropy_group",
    cols     = group_colors
  ) +
    ggplot2::ggtitle(paste0(sample_name, "Top/Bottom ", GROUP_PERCENTAGE, "% Groups")) +
    ggplot2::labs(fill = "Entropy Group") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13)
    )
  
  combined_plot <- p1 | p2 | p3
  
  ggplot2::ggsave(
    filename = fig_file,
    plot     = combined_plot,
    width    = 18,
    height   = 6,
    dpi      = 300
  )
  
  rm(seurat_obj, entropy_df, groups_df, combined_plot, p1, p2, p3)
  gc(verbose = FALSE)
  
  message("-> Done with ", sample_name, "!")
}

if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  
  summary_csv <- file.path(dir_tables, "entropy_groups_summary.csv")
  write.csv(summary_df, file = summary_csv, row.names = FALSE)
  message("\nSaved cohort group summary table to: ", summary_csv)
}

