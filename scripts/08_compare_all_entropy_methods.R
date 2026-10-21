#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(here)
})

dir_objects_01 <- here::here("results", "01_load_and_compute_entropy", "objects")
dir_tables_01  <- here::here("results", "01_load_and_compute_entropy", "tables")
dir_tables_03  <- here::here("results", "03_compute_alternative_entropies", "tables")
dir_tables_06  <- here::here("results", "06_compute_spanorm_entropy", "tables")

dir_results <- here::here("results", "08_compare_all_entropy_methods")
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

all_csv_01  <- list.files(dir_tables_01, pattern = "_entropy\\.csv$", full.names = FALSE)
all_samples <- gsub("_entropy\\.csv$", "", all_csv_01)

sample_numbers <- suppressWarnings(as.numeric(gsub("[^0-9]", "", all_samples)))
if (!any(is.na(sample_numbers))) {
  all_samples <- all_samples[order(sample_numbers)]
}

if (identical(SAMPLES_TO_RUN, "all")) {
  samples_to_process <- all_samples
} else {
  samples_to_process <- SAMPLES_TO_RUN
}

message("=== Starting Analysis: 08_compare_all_entropy_methods ===")
message("Samples selected (", length(samples_to_process), "): ",
        paste(samples_to_process, collapse = ", "))

summary_list <- list()

method_levels <- c(
  "Plug-in (Raw)",
  "Chao-Shen",
  "Chao-Wang-Jost (CWJ)",
  "SpaNorm"
)

method_colors <- c(
  "Plug-in (Raw)"        = "#2b5c8f",
  "Chao-Shen"            = "#d95f02",
  "Chao-Wang-Jost (CWJ)" = "#1b9e77",
  "SpaNorm"              = "#7570b3"
)

spatial_palette <- rev(RColorBrewer::brewer.pal(11, "Spectral"))


for (sample_name in samples_to_process) {
  message("\n------------------------------------------------------------")
  message("Processing: ", sample_name)
  message("------------------------------------------------------------")
  
  rds_path    <- file.path(dir_objects_01, paste0(sample_name, "_seurat.rds"))
  csv_01_path <- file.path(dir_tables_01,  paste0(sample_name, "_entropy.csv"))
  csv_03_path <- file.path(dir_tables_03,  paste0(sample_name, "_alternative_entropies.csv"))
  csv_06_path <- file.path(dir_tables_06,  paste0(sample_name, "_spanorm_entropy.csv"))
  
  if (!file.exists(csv_06_path)) {
    csv_06_path <- file.path(dir_tables_06, paste0(sample_name, "_entropy.csv"))
  }
  
  missing_inputs <- c(
    if (!file.exists(rds_path))    "Script 01 Seurat RDS",
    if (!file.exists(csv_01_path)) "Script 01 Plug-in CSV",
    if (!file.exists(csv_03_path)) "Script 03 Chao-Shen/CWJ CSV",
    if (!file.exists(csv_06_path)) "Script 06 SpaNorm CSV"
  )
  
  if (length(missing_inputs) > 0) {
    warning("Skipping ", sample_name, " because required inputs are missing: ",
            paste(missing_inputs, collapse = ", "))
    next
  }
  
  df_01 <- read.csv(csv_01_path, stringsAsFactors = FALSE)
  df_03 <- read.csv(csv_03_path, stringsAsFactors = FALSE)
  df_06 <- read.csv(csv_06_path, stringsAsFactors = FALSE)
  
  spanorm_col <- if ("spanorm_entropy" %in% colnames(df_06)) "spanorm_entropy" else "entropy"
  
  merged_df <- data.frame(
    spot_id           = df_01$spot_id,
    nCount_RNA        = df_01$nCount_RNA,
    nFeature_RNA      = df_01$nFeature_RNA,
    plugin_entropy    = df_01$entropy,
    chao_shen_entropy = df_03$chao_shen_entropy[match(df_01$spot_id, df_03$spot_id)],
    cwj_entropy       = df_03$cwj_entropy[match(df_01$spot_id, df_03$spot_id)],
    spanorm_entropy   = df_06[[spanorm_col]][match(df_01$spot_id, df_06$spot_id)],
    stringsAsFactors  = FALSE
  )
  
  valid_mask <- merged_df$nCount_RNA > 0 &
    !is.na(merged_df$plugin_entropy) &
    !is.na(merged_df$chao_shen_entropy) &
    !is.na(merged_df$cwj_entropy) &
    !is.na(merged_df$spanorm_entropy)
  
  valid_df <- merged_df[valid_mask, ]
  n_spots  <- nrow(valid_df)
  
  message("-> Valid spots across all 4 methods: ", n_spots)
  
  cor_pi_ncount  <- cor.test(valid_df$nCount_RNA, valid_df$plugin_entropy,    method = "spearman", exact = FALSE)
  cor_cs_ncount  <- cor.test(valid_df$nCount_RNA, valid_df$chao_shen_entropy, method = "spearman", exact = FALSE)
  cor_cwj_ncount <- cor.test(valid_df$nCount_RNA, valid_df$cwj_entropy,       method = "spearman", exact = FALSE)
  cor_spa_ncount <- cor.test(valid_df$nCount_RNA, valid_df$spanorm_entropy,   method = "spearman", exact = FALSE)
  
  cor_pi_nfeat   <- cor.test(valid_df$nFeature_RNA, valid_df$plugin_entropy,    method = "spearman", exact = FALSE)
  cor_cs_nfeat   <- cor.test(valid_df$nFeature_RNA, valid_df$chao_shen_entropy, method = "spearman", exact = FALSE)
  cor_cwj_nfeat  <- cor.test(valid_df$nFeature_RNA, valid_df$cwj_entropy,       method = "spearman", exact = FALSE)
  cor_spa_nfeat  <- cor.test(valid_df$nFeature_RNA, valid_df$spanorm_entropy,   method = "spearman", exact = FALSE)
  
  rho_pi_nc  <- unname(cor_pi_ncount$estimate)
  rho_cs_nc  <- unname(cor_cs_ncount$estimate)
  rho_cwj_nc <- unname(cor_cwj_ncount$estimate)
  rho_spa_nc <- unname(cor_spa_ncount$estimate)
  
  message(sprintf(
    "-> Spearman rho (nCount_RNA): Plug-in = %.3f | Chao-Shen = %.3f | CWJ = %.3f | SpaNorm = %.3f",
    rho_pi_nc, rho_cs_nc, rho_cwj_nc, rho_spa_nc
  ))
  
  summary_list[[sample_name]] <- data.frame(
    sample                     = sample_name,
    n_spots                    = n_spots,
    mean_nCount                = mean(valid_df$nCount_RNA),
    mean_nFeature              = mean(valid_df$nFeature_RNA),
    mean_plugin_entropy        = mean(valid_df$plugin_entropy),
    sd_plugin_entropy          = sd(valid_df$plugin_entropy),
    mean_chao_shen_entropy     = mean(valid_df$chao_shen_entropy),
    sd_chao_shen_entropy       = sd(valid_df$chao_shen_entropy),
    mean_cwj_entropy           = mean(valid_df$cwj_entropy),
    sd_cwj_entropy             = sd(valid_df$cwj_entropy),
    mean_spanorm_entropy       = mean(valid_df$spanorm_entropy),
    sd_spanorm_entropy         = sd(valid_df$spanorm_entropy),
    rho_nCount_plugin          = rho_pi_nc,
    rho_nCount_chao_shen       = rho_cs_nc,
    rho_nCount_cwj             = rho_cwj_nc,
    rho_nCount_spanorm         = rho_spa_nc,
    rho_nFeature_plugin        = unname(cor_pi_nfeat$estimate),
    rho_nFeature_chao_shen     = unname(cor_cs_nfeat$estimate),
    rho_nFeature_cwj           = unname(cor_cwj_nfeat$estimate),
    rho_nFeature_spanorm       = unname(cor_spa_nfeat$estimate),
    stringsAsFactors           = FALSE
  )
  
  message("-> Loading Seurat object for 4-panel spatial comparison plot...")
  seurat_obj <- readRDS(rds_path)
  
  spot_match <- match(colnames(seurat_obj), merged_df$spot_id)
  seurat_obj$plugin_entropy    <- merged_df$plugin_entropy[spot_match]
  seurat_obj$chao_shen_entropy <- merged_df$chao_shen_entropy[spot_match]
  seurat_obj$cwj_entropy       <- merged_df$cwj_entropy[spot_match]
  seurat_obj$spanorm_entropy   <- merged_df$spanorm_entropy[spot_match]
  
  all_entropy_values <- c(
    valid_df$plugin_entropy,
    valid_df$chao_shen_entropy,
    valid_df$cwj_entropy,
    valid_df$spanorm_entropy
  )
  common_limits <- range(all_entropy_values, na.rm = TRUE)
  
  make_spatial_panel <- function(feature_col, panel_title, rho_val) {
    Seurat::SpatialFeaturePlot(seurat_obj, features = feature_col) +
      ggplot2::scale_fill_gradientn(
        colours = spatial_palette,
        limits  = common_limits,
        name    = "Entropy\n(bits)"
      ) +
      ggplot2::labs(
        title    = panel_title,
        subtitle = sprintf("Mean = %.2f bits | rho(nCount) = %.3f",
                           mean(valid_df[[feature_col]], na.rm = TRUE), rho_val)
      ) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
        legend.position = "right"
      )
  }
  
  p_plugin  <- make_spatial_panel("plugin_entropy",    "1. Plug-in Entropy (Raw)",        rho_pi_nc)
  p_cs      <- make_spatial_panel("chao_shen_entropy", "2. Chao-Shen Entropy",            rho_cs_nc)
  p_cwj     <- make_spatial_panel("cwj_entropy",       "3. Chao-Wang-Jost (CWJ) Entropy", rho_cwj_nc)
  p_spanorm <- make_spatial_panel("spanorm_entropy",   "4. SpaNorm Entropy",              rho_spa_nc)
  
  spatial_grid <- ((p_plugin | p_cs) / (p_cwj | p_spanorm)) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title    = paste0(sample_name, ": Spatial Comparison of Entropy Methods (Shared Scale)"),
      subtitle = sprintf("Common color scale across all 4 methods: [%.2f, %.2f] bits (n = %d spots)",
                         common_limits[1], common_limits[2], n_spots),
      theme    = ggplot2::theme(
        plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12)
      )
    )
  
  fig_spatial_path <- file.path(dir_figures, paste0(sample_name, "_spatial_4methods_comparison.png"))
  ggplot2::ggsave(
    filename = fig_spatial_path,
    plot     = spatial_grid,
    width    = 13,
    height   = 11,
    dpi      = 300
  )
  message("-> Saved 4-panel spatial comparison plot to: ", fig_spatial_path)
  
  rm(seurat_obj, merged_df, valid_df, spatial_grid, p_plugin, p_cs, p_cwj, p_spanorm)
  gc(verbose = FALSE)
}


if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  
  sample_csv_path <- file.path(dir_tables, "all_methods_sample_summary.csv")
  write.csv(summary_df, file = sample_csv_path, row.names = FALSE)
  message("\n-> Saved per-sample summary table to: ", sample_csv_path)
  
  cohort_overview <- data.frame(
    method                  = method_levels,
    mean_entropy_across_cohort = c(
      mean(summary_df$mean_plugin_entropy),
      mean(summary_df$mean_chao_shen_entropy),
      mean(summary_df$mean_cwj_entropy),
      mean(summary_df$mean_spanorm_entropy)
    ),
    mean_rho_nCount         = c(
      mean(summary_df$rho_nCount_plugin),
      mean(summary_df$rho_nCount_chao_shen),
      mean(summary_df$rho_nCount_cwj),
      mean(summary_df$rho_nCount_spanorm)
    ),
    sd_rho_nCount           = c(
      sd(summary_df$rho_nCount_plugin),
      sd(summary_df$rho_nCount_chao_shen),
      sd(summary_df$rho_nCount_cwj),
      sd(summary_df$rho_nCount_spanorm)
    ),
    median_rho_nCount       = c(
      median(summary_df$rho_nCount_plugin),
      median(summary_df$rho_nCount_chao_shen),
      median(summary_df$rho_nCount_cwj),
      median(summary_df$rho_nCount_spanorm)
    ),
    min_rho_nCount          = c(
      min(summary_df$rho_nCount_plugin),
      min(summary_df$rho_nCount_chao_shen),
      min(summary_df$rho_nCount_cwj),
      min(summary_df$rho_nCount_spanorm)
    ),
    max_rho_nCount          = c(
      max(summary_df$rho_nCount_plugin),
      max(summary_df$rho_nCount_chao_shen),
      max(summary_df$rho_nCount_cwj),
      max(summary_df$rho_nCount_spanorm)
    ),
    mean_abs_rho_nCount     = c(
      mean(abs(summary_df$rho_nCount_plugin)),
      mean(abs(summary_df$rho_nCount_chao_shen)),
      mean(abs(summary_df$rho_nCount_cwj)),
      mean(abs(summary_df$rho_nCount_spanorm))
    ),
    mean_rho_nFeature       = c(
      mean(summary_df$rho_nFeature_plugin),
      mean(summary_df$rho_nFeature_chao_shen),
      mean(summary_df$rho_nFeature_cwj),
      mean(summary_df$rho_nFeature_spanorm)
    ),
    stringsAsFactors = FALSE
  )
  
  cohort_csv_path <- file.path(dir_tables, "all_methods_cohort_overview.csv")
  write.csv(cohort_overview, file = cohort_csv_path, row.names = FALSE)
  message("-> Saved cohort-wide overview table to: ", cohort_csv_path)
  
  summary_df$sample <- factor(summary_df$sample, levels = summary_df$sample)
  
  plot_df <- rbind(
    data.frame(
      sample = summary_df$sample,
      method = "Plug-in (Raw)",
      rho    = summary_df$rho_nCount_plugin
    ),
    data.frame(
      sample = summary_df$sample,
      method = "Chao-Shen",
      rho    = summary_df$rho_nCount_chao_shen
    ),
    data.frame(
      sample = summary_df$sample,
      method = "Chao-Wang-Jost (CWJ)",
      rho    = summary_df$rho_nCount_cwj
    ),
    data.frame(
      sample = summary_df$sample,
      method = "SpaNorm",
      rho    = summary_df$rho_nCount_spanorm
    )
  )
  
  plot_df$method <- factor(plot_df$method, levels = method_levels)
  
  p_barplot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample, y = rho, fill = method)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.85), width = 0.78) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
    ggplot2::scale_fill_manual(values = method_colors) +
    ggplot2::labs(
      title    = "Comparison of Sequencing Depth Correlation (nCount_RNA) Across All Entropy Methods",
      subtitle = "Spearman rank correlation (rho) between spot entropy and total UMI counts per sample",
      x        = "Sample",
      y        = "Spearman Rank Correlation with nCount_RNA (rho)",
      fill     = "Entropy Method"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
  
  barplot_path <- file.path(dir_figures, "cohort_ncount_correlation_barplot.png")
  ggplot2::ggsave(
    filename = barplot_path,
    plot     = p_barplot,
    width    = 15,
    height   = 6.5,
    dpi      = 300
  )
  message("-> Saved comprehensive cohort correlation bar plot to: ", barplot_path)
  
  if (nrow(summary_df) > 1) {
    cohort_overview$method <- factor(cohort_overview$method, levels = method_levels)
    
    p_mean_bar <- ggplot2::ggplot(cohort_overview, ggplot2::aes(x = method, y = mean_rho_nCount, fill = method)) +
      ggplot2::geom_col(width = 0.6, alpha = 0.9) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = mean_rho_nCount - sd_rho_nCount,
          ymax = mean_rho_nCount + sd_rho_nCount
        ),
        width = 0.2,
        linewidth = 0.7
      ) +
      ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
      ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
      ggplot2::scale_fill_manual(values = method_colors) +
      ggplot2::labs(
        title    = "Mean Correlation Across Cohort (+/- SD)",
        subtitle = "Average Spearman rho with nCount_RNA",
        x        = NULL,
        y        = "Mean Spearman rho"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        axis.text.x     = ggplot2::element_text(angle = 20, hjust = 1),
        legend.position = "none"
      )
    
    p_box <- ggplot2::ggplot(plot_df, ggplot2::aes(x = method, y = rho, fill = method)) +
      ggplot2::geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
      ggplot2::geom_jitter(width = 0.12, size = 2, alpha = 0.85, color = "black") +
      ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.6, linetype = "dashed") +
      ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
      ggplot2::scale_fill_manual(values = method_colors) +
      ggplot2::labs(
        title    = "Distribution of Sample Correlations",
        subtitle = "Each point represents one sample (n = 22)",
        x        = NULL,
        y        = "Spearman rho with nCount_RNA"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        axis.text.x     = ggplot2::element_text(angle = 20, hjust = 1),
        legend.position = "none"
      )
    
    p_summary_combined <- p_mean_bar | p_box
    
    summary_fig_path <- file.path(dir_figures, "cohort_ncount_correlation_summary_boxplot.png")
    ggplot2::ggsave(
      filename = summary_fig_path,
      plot     = p_summary_combined,
      width    = 13,
      height   = 6,
      dpi      = 300
    )
    message("-> Saved cohort correlation summary figure to: ", summary_fig_path)
  }
}

