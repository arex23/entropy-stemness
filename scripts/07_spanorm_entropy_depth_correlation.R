#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(here)
})

dir_input_spanorm <- here::here("results", "06_compute_spanorm_entropy", "tables")
dir_input_raw     <- here::here("results", "01_load_and_compute_entropy", "tables")

dir_results <- here::here("results", "07_spanorm_entropy_depth_correlation")
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)


all_csv_files <- list.files(dir_input_spanorm, pattern = "_spanorm_entropy\\.csv$", full.names = FALSE)
if (length(all_csv_files) == 0) {
  all_csv_files <- list.files(dir_input_spanorm, pattern = "_entropy\\.csv$", full.names = FALSE)
}

all_samples <- gsub("(_spanorm)?_entropy\\.csv$", "", all_csv_files)

sample_numbers <- suppressWarnings(as.numeric(gsub("[^0-9]", "", all_samples)))
if (!any(is.na(sample_numbers))) {
  all_samples <- all_samples[order(sample_numbers)]
}

if (identical(SAMPLES_TO_RUN, "all")) {
  samples_to_process <- all_samples
} else {
  missing_samples <- setdiff(SAMPLES_TO_RUN, all_samples)
  if (length(missing_samples) > 0) {
    stop("The following requested samples were not found in ", dir_input_spanorm, ": ",
         paste(missing_samples, collapse = ", "))
  }
  samples_to_process <- SAMPLES_TO_RUN
}

message("=== Starting Analysis: 07_spanorm_entropy_depth_correlation ===")
message("Samples selected (", length(samples_to_process), "): ",
        paste(samples_to_process, collapse = ", "))

summary_list <- list()


for (sample_name in samples_to_process) {
  message("\n------------------------------------------------------------")
  message("Analyzing: ", sample_name)
  message("------------------------------------------------------------")
  
  spanorm_csv_path <- file.path(dir_input_spanorm, paste0(sample_name, "_spanorm_entropy.csv"))
  if (!file.exists(spanorm_csv_path)) {
    spanorm_csv_path <- file.path(dir_input_spanorm, paste0(sample_name, "_entropy.csv"))
  }
  if (!file.exists(spanorm_csv_path)) {
    warning("SpaNorm entropy CSV not found for ", sample_name, ": ", spanorm_csv_path, ". Skipping...")
    next
  }
  
  df_spanorm <- read.csv(spanorm_csv_path, stringsAsFactors = FALSE)
  spanorm_col <- if ("spanorm_entropy" %in% colnames(df_spanorm)) "spanorm_entropy" else "entropy"
  df <- df_spanorm
  
  raw_csv_path <- file.path(dir_input_raw, paste0(sample_name, "_entropy.csv"))
  if (file.exists(raw_csv_path)) {
    df_raw <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
    df <- merge(
      df,
      df_raw[, c("spot_id", "entropy")],
      by = "spot_id",
      all.x = TRUE
    )
    names(df)[names(df) == "entropy.y"] <- "raw_entropy"
    names(df)[names(df) == "entropy.x"] <- "spanorm_entropy_alias"
  } else {
    df$raw_entropy <- NA_real_
  }
  
  valid_df <- df[!is.na(df[[spanorm_col]]) & df$nCount_RNA > 0, ]
  n_spots  <- nrow(valid_df)
  
  message("-> Spots with valid SpaNorm entropy: ", n_spots)
  
  cor_spanorm_ncount <- cor.test(valid_df$nCount_RNA, valid_df[[spanorm_col]], method = "spearman", exact = FALSE)
  cor_spanorm_nfeat  <- cor.test(valid_df$nFeature_RNA, valid_df[[spanorm_col]], method = "spearman", exact = FALSE)
  
  rho_spanorm_ncount <- unname(cor_spanorm_ncount$estimate)
  rho_spanorm_nfeat  <- unname(cor_spanorm_nfeat$estimate)
  
  if (!all(is.na(valid_df$raw_entropy))) {
    cor_raw_ncount <- cor.test(valid_df$nCount_RNA, valid_df$raw_entropy, method = "spearman", exact = FALSE)
    cor_raw_nfeat  <- cor.test(valid_df$nFeature_RNA, valid_df$raw_entropy, method = "spearman", exact = FALSE)
    rho_raw_ncount <- unname(cor_raw_ncount$estimate)
    p_raw_ncount   <- cor_raw_ncount$p.value
    rho_raw_nfeat  <- unname(cor_raw_nfeat$estimate)
    p_raw_nfeat    <- cor_raw_nfeat$p.value
  } else {
    rho_raw_ncount <- NA_real_
    p_raw_ncount   <- NA_real_
    rho_raw_nfeat  <- NA_real_
    p_raw_nfeat    <- NA_real_
  }
  
  message(sprintf("-> Spearman rho with nCount_RNA:   SpaNorm = %.3f | Raw Baseline = %.3f",
                  rho_spanorm_ncount, rho_raw_ncount))
  message(sprintf("-> Spearman rho with nFeature_RNA: SpaNorm = %.3f | Raw Baseline = %.3f",
                  rho_spanorm_nfeat, rho_raw_nfeat))
  
  summary_list[[sample_name]] <- data.frame(
    sample                        = sample_name,
    n_spots                       = n_spots,
    mean_spanorm_entropy          = mean(valid_df[[spanorm_col]]),
    sd_spanorm_entropy            = sd(valid_df[[spanorm_col]]),
    mean_nCount                   = mean(valid_df$nCount_RNA),
    mean_nFeature                 = mean(valid_df$nFeature_RNA),
    spearman_rho_spanorm_nCount   = rho_spanorm_ncount,
    spearman_p_spanorm_nCount     = cor_spanorm_ncount$p.value,
    spearman_rho_raw_nCount       = rho_raw_ncount,
    spearman_p_raw_nCount         = p_raw_ncount,
    delta_rho_nCount              = rho_spanorm_ncount - rho_raw_ncount,
    spearman_rho_spanorm_nFeature = rho_spanorm_nfeat,
    spearman_p_spanorm_nFeature   = cor_spanorm_nfeat$p.value,
    spearman_rho_raw_nFeature     = rho_raw_nfeat,
    spearman_p_raw_nFeature       = p_raw_nfeat,
    stringsAsFactors              = FALSE
  )
  
  sub_ncount <- if (!is.na(rho_raw_ncount)) {
    sprintf("Spearman rho = %.3f (vs Raw Baseline rho = %.3f)", rho_spanorm_ncount, rho_raw_ncount)
  } else {
    sprintf("Spearman rho = %.3f (n = %d spots)", rho_spanorm_ncount, n_spots)
  }
  
  sub_nfeat <- if (!is.na(rho_raw_nfeat)) {
    sprintf("Spearman rho = %.3f (vs Raw Baseline rho = %.3f)", rho_spanorm_nfeat, rho_raw_nfeat)
  } else {
    sprintf("Spearman rho = %.3f", rho_spanorm_nfeat)
  }
  
  p_ncount <- ggplot(valid_df, aes(x = nCount_RNA, y = .data[[spanorm_col]])) +
    geom_point(alpha = 0.25, size = 0.8, color = "#7570b3") +
    geom_smooth(method = "lm", formula = y ~ log10(x), color = "#d95f02", se = FALSE, linewidth = 1) +
    labs(
      title    = paste0(sample_name, ": SpaNorm Entropy vs Total Counts"),
      subtitle = sub_ncount,
      x        = "Total UMI Counts per Spot (nCount_RNA)",
      y        = "SpaNorm Shannon Entropy (bits)"
    ) +
    theme_bw(base_size = 12)
  
  p_nfeat <- ggplot(valid_df, aes(x = nFeature_RNA, y = .data[[spanorm_col]])) +
    geom_point(alpha = 0.25, size = 0.8, color = "#1b9e77") +
    stat_function(
      fun      = log2,
      aes(color = "Theoretical Max: log2(nFeature)"),
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_color_manual(name = NULL, values = c("Theoretical Max: log2(nFeature)" = "#e7298a")) +
    labs(
      title    = paste0(sample_name, ": SpaNorm Entropy vs Detected Genes"),
      subtitle = sub_nfeat,
      x        = "Detected Genes per Spot (nFeature_RNA)",
      y        = "SpaNorm Shannon Entropy (bits)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  
  combined_fig <- p_ncount | p_nfeat
  
  fig_path <- file.path(dir_figures, paste0(sample_name, "_spanorm_entropy_vs_depth.png"))
  ggsave(
    filename = fig_path,
    plot     = combined_fig,
    width    = 13,
    height   = 6,
    dpi      = 300
  )
  message("-> Saved diagnostic plot to: ", fig_path)
}


if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  
  summary_csv_path <- file.path(dir_tables, "correlation_summary.csv")
  write.csv(summary_df, file = summary_csv_path, row.names = FALSE)
  message("\n-> Saved correlation summary table to: ", summary_csv_path)
  
  if (nrow(summary_df) > 1) {
    summary_df$sample <- factor(summary_df$sample, levels = summary_df$sample)
    
    plot_df_comparison <- rbind(
      data.frame(
        sample    = summary_df$sample,
        estimator = "Raw Data Entropy",
        rho       = summary_df$spearman_rho_raw_nCount
      ),
      data.frame(
        sample    = summary_df$sample,
        estimator = "SpaNorm Entropy",
        rho       = summary_df$spearman_rho_spanorm_nCount
      )
    )
    plot_df_comparison <- plot_df_comparison[!is.na(plot_df_comparison$rho), ]
    plot_df_comparison$estimator <- factor(
      plot_df_comparison$estimator,
      levels = c("Raw Data Entropy", "SpaNorm Entropy")
    )
    
    p_cohort_comparison <- ggplot(plot_df_comparison, aes(x = sample, y = rho, fill = estimator)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
      scale_fill_manual(values = c(
        "Raw Data Entropy" = "#2b5c8f",
        "SpaNorm Entropy"  = "#7570b3"
      )) +
      labs(
        title    = "Cohort Overview: Spearman Correlation with Total Counts (nCount_RNA)",
        subtitle = "Comparing Raw Data Entropy vs. SpaNorm Normalized Entropy across samples",
        x        = "Sample",
        y        = "Spearman Rank Correlation (rho)",
        fill     = "Entropy Calculation"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "top"
      )
    
    cohort_comp_fig_path <- file.path(dir_figures, "cohort_spanorm_correlations_overview.png")
    ggsave(
      filename = cohort_comp_fig_path,
      plot     = p_cohort_comparison,
      width    = 13,
      height   = 6,
      dpi      = 300
    )
    message("-> Saved cohort comparison plot to: ", cohort_comp_fig_path)
    
    plot_df_metrics <- rbind(
      data.frame(
        sample = summary_df$sample,
        metric = "nCount_RNA (Spearman rho)",
        value  = summary_df$spearman_rho_spanorm_nCount
      ),
      data.frame(
        sample = summary_df$sample,
        metric = "nFeature_RNA (Spearman rho)",
        value  = summary_df$spearman_rho_spanorm_nFeature
      )
    )
    
    p_cohort_metrics <- ggplot(plot_df_metrics, aes(x = sample, y = value, fill = metric)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
      scale_fill_manual(values = c(
        "nCount_RNA (Spearman rho)"   = "#7570b3",
        "nFeature_RNA (Spearman rho)" = "#1b9e77"
      )) +
      labs(
        title    = "Cohort Overview: SpaNorm Entropy vs Technical Variables",
        subtitle = "Comparing residual correlation with library size (nCount_RNA) vs detected genes (nFeature_RNA)",
        x        = "Sample",
        y        = "Spearman Rank Correlation (rho)",
        fill     = "Technical Variable"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "top"
      )
    
    cohort_metrics_fig_path <- file.path(dir_figures, "cohort_spanorm_depth_metrics_overview.png")
    ggsave(
      filename = cohort_metrics_fig_path,
      plot     = p_cohort_metrics,
      width    = 12,
      height   = 6,
      dpi      = 300
    )
    message("-> Saved cohort SpaNorm metrics plot to: ", cohort_metrics_fig_path)
  }
}
