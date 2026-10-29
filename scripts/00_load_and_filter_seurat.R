# Purpose: Load 10x Visium spatial transcriptomics data into Seurat objects,
#          save the unfiltered objects, apply a 10% spot-detection gene filter
#          (matching SpaNorm's default cutoff), analyze which genes and
#          meningioma stemness markers are lost, and save the filtered objects.
# Input:   Raw 10x Visium directories in data/raw/<sample>/
# Output:  results/00_load_and_filter_seurat/objects/unfiltered/<sample>_seurat.rds
#          results/00_load_and_filter_seurat/objects/filtered/<sample>_seurat.rds
#          results/00_load_and_filter_seurat/tables/filtering_summary_by_sample.csv
#          results/00_load_and_filter_seurat/tables/stemness_markers_filtering_detail.csv
#          results/00_load_and_filter_seurat/tables/stemness_markers_cohort_summary.csv
#          results/00_load_and_filter_seurat/figures/cohort_genes_kept_vs_lost.png
#          results/00_load_and_filter_seurat/figures/stemness_markers_filtering_heatmap.png

# 1. Configuration and dependencies

SAMPLES_TO_RUN <- "sample1"
#SAMPLES_TO_RUN <- "all"

# Minimum proportion of spots in which a gene must be detected (count > 0)
MIN_SPOT_PROP <- 0.10

# Literature-curated meningioma stemness-associated markers to track during filtering
meningioma_stemness_all <- c(
  "NOTCH3",
  "THY1",
  "HEY1",
  "CXCR4",
  "NANOG",
  "POU5F1",
  "SPP1",
  "PROM1",
  "CDH2",
  "NES",
  "SOX2",
  "ACKR3",
  "CD44",
  "JAG1",
  "JAG2",
  "DLL1",
  "DLL4",
  "FN1",
  "PDGFRB",
  "RGS5",
  "CSPG4",
  "CD248",
  "ABCC9",
  "GJA4",
  "SDC2",
  "ACTA2",
  "SERPING1",
  "CLU",
  "CXCL12",
  "CXCL11",
  "VIM"
)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(here)
})

dir_raw <- here::here("data", "raw")

dir_results            <- here::here("results", "00_load_and_filter_seurat")
dir_objects_unfiltered <- file.path(dir_results, "objects", "unfiltered")
dir_objects_filtered   <- file.path(dir_results, "objects", "filtered")
dir_tables             <- file.path(dir_results, "tables")
dir_figures            <- file.path(dir_results, "figures")

dir.create(dir_objects_unfiltered, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_objects_filtered,   recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables,             recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures,            recursive = TRUE, showWarnings = FALSE)

all_sample_dirs <- list.dirs(dir_raw, full.names = FALSE, recursive = FALSE)
all_sample_dirs <- grep("^sample", all_sample_dirs, value = TRUE)

sample_numbers  <- as.numeric(gsub("sample", "", all_sample_dirs))
all_sample_dirs <- all_sample_dirs[order(sample_numbers)]

if (identical(SAMPLES_TO_RUN, "all")) {
  samples_to_process <- all_sample_dirs
} else {
  samples_to_process <- SAMPLES_TO_RUN
}

message("=== Starting Analysis: 00_load_and_filter_seurat ===")
message("Minimum spot detection proportion: ", MIN_SPOT_PROP)
message("Samples selected: ", paste(samples_to_process, collapse = ", "))

sample_summary_list  <- list()
stemness_detail_list <- list()

# 2. Processing loop

for (sample_name in samples_to_process) {
  message("\n------------------------------------------------------------")
  message("Processing: ", sample_name)
  message("------------------------------------------------------------")

  # A. Load raw 10x Visium data into a Seurat object
  sample_dir <- file.path(dir_raw, sample_name)

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

  # Keep only spots matched to the tissue image coordinates
  seurat_obj <- subset(seurat_obj, cells = Cells(seurat_obj[[Images(seurat_obj)[1]]]))
  seurat_obj$sample <- sample_name

  # B. Save the unfiltered Seurat object (used as input by 05_spanorm_normalization.R)
  rds_unfiltered <- file.path(dir_objects_unfiltered, paste0(sample_name, "_seurat.rds"))
  message("-> Saving UNFILTERED Seurat object to: ", rds_unfiltered)
  saveRDS(seurat_obj, file = rds_unfiltered)

  # C. Evaluate gene detection proportions across spots
  counts_mat    <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  n_genes_total <- nrow(counts_mat)
  n_spots_total <- ncol(counts_mat)

  gene_spots_detected <- Matrix::rowSums(counts_mat > 0)
  gene_total_counts   <- Matrix::rowSums(counts_mat)
  spot_prop           <- gene_spots_detected / n_spots_total

  keep_genes          <- names(spot_prop)[spot_prop >= MIN_SPOT_PROP]
  lost_genes          <- names(spot_prop)[spot_prop <  MIN_SPOT_PROP]
  zero_count_genes    <- names(spot_prop)[gene_spots_detected == 0]
  filtered_rare_genes <- names(spot_prop)[gene_spots_detected > 0 & spot_prop < MIN_SPOT_PROP]

  n_genes_detected        <- sum(gene_spots_detected > 0)
  n_genes_zero            <- length(zero_count_genes)
  n_genes_filtered_active <- length(filtered_rare_genes)

  if (length(keep_genes) == 0) {
    stop("No genes passed the filter for ", sample_name)
  }

  # D. Check status of the literature-curated stemness markers
  marker_status_df <- data.frame(
    sample         = sample_name,
    gene           = meningioma_stemness_all,
    in_dataset     = meningioma_stemness_all %in% rownames(counts_mat),
    spots_detected = NA_integer_,
    n_spots_total  = n_spots_total,
    spot_prop      = NA_real_,
    total_counts   = NA_real_,
    status         = "Not in dataset",
    stringsAsFactors = FALSE
  )

  for (i in seq_along(meningioma_stemness_all)) {
    g <- meningioma_stemness_all[i]
    if (g %in% rownames(counts_mat)) {
      prop_g   <- unname(spot_prop[g])
      spots_g  <- as.integer(unname(gene_spots_detected[g]))
      counts_g <- unname(gene_total_counts[g])

      marker_status_df$spots_detected[i] <- spots_g
      marker_status_df$spot_prop[i]      <- prop_g
      marker_status_df$total_counts[i]   <- counts_g
      if (spots_g == 0) {
        marker_status_df$status[i] <- "Zero counts (0 spots)"
      } else if (prop_g >= MIN_SPOT_PROP) {
        marker_status_df$status[i] <- "Kept"
      } else {
        marker_status_df$status[i] <- "Filtered out"
      }
    }
  }

  stemness_detail_list[[sample_name]] <- marker_status_df

  stemness_kept    <- marker_status_df$gene[marker_status_df$status == "Kept"]
  stemness_lost    <- marker_status_df$gene[marker_status_df$status == "Filtered out"]
  stemness_zero    <- marker_status_df$gene[marker_status_df$status == "Zero counts (0 spots)"]
  stemness_missing <- marker_status_df$gene[marker_status_df$status == "Not in dataset"]

  total_counts_before <- sum(gene_total_counts)
  counts_lost_active  <- sum(gene_total_counts[filtered_rare_genes])

  message(sprintf(
    "-> Total panel genes: %d | Detected in >=1 spot: %d | Zero counts everywhere: %d (%.1f%% of panel)",
    n_genes_total,
    n_genes_detected,
    n_genes_zero,
    100 * n_genes_zero / n_genes_total
  ))
  message(sprintf(
    "-> Breakdown of %d removed genes: %d had 0 counts (no entropy effect) + %d were active rare genes (<10%% spots)",
    length(lost_genes),
    n_genes_zero,
    n_genes_filtered_active
  ))
  message(sprintf(
    "-> Entropy-relevant filter impact: keeping %d / %d active genes (%.1f%% of active genes removed, representing %.2f%% of total UMI counts)",
    length(keep_genes),
    n_genes_detected,
    100 * n_genes_filtered_active / n_genes_detected,
    100 * counts_lost_active / total_counts_before
  ))
  message(sprintf(
    "-> Stemness markers KEPT (>=10%% spots, %d/%d): %s",
    length(stemness_kept),
    length(meningioma_stemness_all),
    ifelse(length(stemness_kept) > 0, paste(stemness_kept, collapse = ", "), "None")
  ))
  message(sprintf(
    "-> Stemness markers FILTERED OUT (active in <10%% spots, %d/%d): %s",
    length(stemness_lost),
    length(meningioma_stemness_all),
    ifelse(length(stemness_lost) > 0, paste(stemness_lost, collapse = ", "), "None")
  ))
  if (length(stemness_zero) > 0) {
    message(sprintf(
      "-> Stemness markers with ZERO COUNTS (0 spots, %d/%d): %s",
      length(stemness_zero),
      length(meningioma_stemness_all),
      paste(stemness_zero, collapse = ", ")
    ))
  }
  if (length(stemness_missing) > 0) {
    message(sprintf(
      "-> Stemness markers NOT IN DATASET (%d/%d): %s",
      length(stemness_missing),
      length(meningioma_stemness_all),
      paste(stemness_missing, collapse = ", ")
    ))
  }

  # E. Subset Seurat object to kept genes and remove any resulting zero-count spots
  rm(counts_mat, gene_spots_detected, gene_total_counts, spot_prop)

  seurat_obj <- subset(seurat_obj, features = keep_genes)
  counts_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")

  spot_totals <- Matrix::colSums(counts_mat)
  if (any(spot_totals == 0)) {
    zero_spots <- colnames(counts_mat)[spot_totals == 0]
    message(
      "-> Removing ", length(zero_spots), " zero-count spot(s) after gene filtering in ",
      sample_name, ": ", paste(zero_spots, collapse = ", ")
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

  n_genes_kept       <- nrow(counts_mat)
  n_spots_kept       <- ncol(counts_mat)
  total_counts_after <- sum(counts_mat)
  rm(counts_mat)

  # F. Save the filtered Seurat object
  rds_filtered <- file.path(dir_objects_filtered, paste0(sample_name, "_seurat.rds"))
  message("-> Saving FILTERED Seurat object to: ", rds_filtered)
  saveRDS(seurat_obj, file = rds_filtered)

  # G. Record per-sample summary metrics
  sample_summary_list[[sample_name]] <- data.frame(
    sample                      = sample_name,
    n_spots_before              = n_spots_total,
    n_spots_after               = n_spots_kept,
    n_genes_total               = n_genes_total,
    n_genes_zero_counts         = n_genes_zero,
    n_genes_detected_nonzero    = n_genes_detected,
    n_genes_kept                = n_genes_kept,
    n_genes_lost_total          = length(lost_genes),
    n_genes_lost_active_nonzero = n_genes_filtered_active,
    pct_total_genes_lost        = 100 * length(lost_genes) / n_genes_total,
    pct_lost_that_were_zero     = 100 * n_genes_zero / length(lost_genes),
    pct_active_genes_lost       = 100 * n_genes_filtered_active / n_genes_detected,
    total_counts_before         = total_counts_before,
    total_counts_after          = total_counts_after,
    pct_counts_lost             = 100 * (total_counts_before - total_counts_after) / total_counts_before,
    n_stemness_kept             = length(stemness_kept),
    n_stemness_lost_active      = length(stemness_lost),
    n_stemness_zero_counts      = length(stemness_zero),
    n_stemness_missing          = length(stemness_missing),
    stemness_genes_kept         = paste(stemness_kept, collapse = ";"),
    stemness_genes_lost_active  = paste(stemness_lost, collapse = ";"),
    stemness_genes_zero_counts  = paste(stemness_zero, collapse = ";"),
    stringsAsFactors            = FALSE
  )

  rm(seurat_obj, keep_genes, lost_genes, zero_count_genes, filtered_rare_genes)
  gc(verbose = FALSE)

  message("-> Done with ", sample_name, "!")
}

# 3. Save summary tables and diagnostic figures

if (length(sample_summary_list) > 0) {
  summary_df <- do.call(rbind, sample_summary_list)
  rownames(summary_df) <- NULL

  summary_csv <- file.path(dir_tables, "filtering_summary_by_sample.csv")
  write.csv(summary_df, file = summary_csv, row.names = FALSE)
  message("\n-> Saved sample filtering summary table to: ", summary_csv)

  stemness_detail_df <- do.call(rbind, stemness_detail_list)
  rownames(stemness_detail_df) <- NULL

  stemness_detail_csv <- file.path(dir_tables, "stemness_markers_filtering_detail.csv")
  write.csv(stemness_detail_df, file = stemness_detail_csv, row.names = FALSE)
  message("-> Saved per-sample stemness marker detail table to: ", stemness_detail_csv)

  # Cohort-level summary for each stemness marker
  marker_cohort_rows <- list()
  for (g in meningioma_stemness_all) {
    sub_g <- stemness_detail_df[stemness_detail_df$gene == g, ]
    marker_cohort_rows[[g]] <- data.frame(
      gene                  = g,
      n_samples_evaluated   = nrow(sub_g),
      n_samples_kept        = sum(sub_g$status == "Kept"),
      n_samples_filtered    = sum(sub_g$status == "Filtered out"),
      n_samples_zero_counts = sum(sub_g$status == "Zero counts (0 spots)"),
      n_samples_not_in_data = sum(sub_g$status == "Not in dataset"),
      mean_spot_prop        = mean(sub_g$spot_prop, na.rm = TRUE),
      min_spot_prop         = ifelse(all(is.na(sub_g$spot_prop)), NA_real_, min(sub_g$spot_prop, na.rm = TRUE)),
      max_spot_prop         = ifelse(all(is.na(sub_g$spot_prop)), NA_real_, max(sub_g$spot_prop, na.rm = TRUE)),
      stringsAsFactors      = FALSE
    )
  }
  marker_cohort_df <- do.call(rbind, marker_cohort_rows)
  rownames(marker_cohort_df) <- NULL

  marker_cohort_csv <- file.path(dir_tables, "stemness_markers_cohort_summary.csv")
  write.csv(marker_cohort_df, file = marker_cohort_csv, row.names = FALSE)
  message("-> Saved cohort stemness marker summary to: ", marker_cohort_csv)

  # Figure 1: Bar plot of kept vs filtered-active vs zero-count genes per sample
  summary_df$sample <- factor(summary_df$sample, levels = summary_df$sample)

  genes_plot_df <- rbind(
    data.frame(
      sample = summary_df$sample,
      group  = "Kept (>= 10% spots)",
      count  = summary_df$n_genes_kept
    ),
    data.frame(
      sample = summary_df$sample,
      group  = "Filtered active (0% < spots < 10%)",
      count  = summary_df$n_genes_lost_active_nonzero
    ),
    data.frame(
      sample = summary_df$sample,
      group  = "Zero counts (0 spots, no entropy effect)",
      count  = summary_df$n_genes_zero_counts
    )
  )
  genes_plot_df$group <- factor(
    genes_plot_df$group,
    levels = c(
      "Zero counts (0 spots, no entropy effect)",
      "Filtered active (0% < spots < 10%)",
      "Kept (>= 10% spots)"
    )
  )

  p_genes <- ggplot(genes_plot_df, aes(x = sample, y = count, fill = group)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = c(
      "Kept (>= 10% spots)"                      = "#2b8f5c",
      "Filtered active (0% < spots < 10%)"       = "#d95f02",
      "Zero counts (0 spots, no entropy effect)" = "#bdbdbd"
    )) +
    labs(
      title    = "Gene Filtering Breakdown Across Samples (10% Spot Detection Cutoff)",
      subtitle = "Distinguishing active rare genes (<10% spots) from unexpressed zero-count genes (0 spots)",
      x        = "Sample",
      y        = "Number of Genes",
      fill     = "Gene Category"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  fig_genes_path <- file.path(dir_figures, "cohort_genes_kept_vs_lost.png")
  ggsave(
    filename = fig_genes_path,
    plot     = p_genes,
    width    = 12,
    height   = 6,
    dpi      = 300
  )
  message("-> Saved gene filtering bar plot to: ", fig_genes_path)

  # Figure 2: Stemness marker status heatmap across processed samples
  stemness_detail_df$sample <- factor(stemness_detail_df$sample, levels = levels(summary_df$sample))
  stemness_detail_df$gene   <- factor(stemness_detail_df$gene, levels = rev(meningioma_stemness_all))
  stemness_detail_df$status <- factor(
    stemness_detail_df$status,
    levels = c("Kept", "Filtered out", "Zero counts (0 spots)", "Not in dataset")
  )
  stemness_detail_df$pct_label <- ifelse(
    is.na(stemness_detail_df$spot_prop),
    "NA",
    sprintf("%.0f%%", 100 * stemness_detail_df$spot_prop)
  )

  p_stemness <- ggplot(stemness_detail_df, aes(x = sample, y = gene, fill = status)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = pct_label), size = 2.7, color = "black") +
    scale_fill_manual(
      values = c(
        "Kept"                  = "#74c476",
        "Filtered out"          = "#fb6a4a",
        "Zero counts (0 spots)" = "#fdd0a2",
        "Not in dataset"        = "#d9d9d9"
      ),
      drop = FALSE
    ) +
    labs(
      title    = "Meningioma Stemness Markers: 10% Spot-Detection Filter Status",
      subtitle = "Cell labels show the percentage of tissue spots where each marker is detected (count > 0)",
      x        = "Sample",
      y        = "Stemness Marker Gene",
      fill     = "Filter Status"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      panel.grid      = element_blank()
    )

  fig_stemness_path <- file.path(dir_figures, "stemness_markers_filtering_heatmap.png")
  ggsave(
    filename = fig_stemness_path,
    plot     = p_stemness,
    width    = max(8, 3 + 0.55 * nrow(summary_df)),
    height   = 9,
    dpi      = 300
  )
  message("-> Saved stemness marker filtering heatmap to: ", fig_stemness_path)
}
