#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(here)
})

dir_input_alt <- here::here("results", "03_compute_alternative_entropies", "tables")
dir_input_raw <- here::here("results", "01_load_and_compute_entropy", "tables")

dir_results <- here::here("results", "04_alternative_entropy_depth_correlation")
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)


all_csv_files <- list.files(dir_input_alt, pattern = "_alternative_entropies\\.csv$", full.names = FALSE)
all_samples   <- gsub("_alternative_entropies\\.csv$", "", all_csv_files)

sample_numbers <- as.numeric(gsub("sample", "", all_samples))
all_samples    <- all_samples[order(sample_numbers)]

if (identical(SAMPLES_TO_RUN, "all")) {
  samples_to_process <- all_samples
} else {
  samples_to_process <- SAMPLES_TO_RUN
}

message("=== Starting Analysis: 04_alternative_entropy_depth_correlation ===")
message("Samples selected: ", paste(samples_to_process, collapse = ", "))

summary_list <- list()


for (sample_name in samples_to_process) {
  message("\n------------------------------------------------------------")
  message("Analyzing: ", sample_name)
  message("------------------------------------------------------------")
  
  alt_csv_path <- file.path(dir_input_alt, paste0(sample_name, "_alternative_entropies.csv"))
  if (!file.exists(alt_csv_path)) {
    warning("Alternative entropy CSV not found for ", sample_name, ". Skipping...")
    next
  }
  
  df_alt <- read.csv(alt_csv_path, stringsAsFactors = FALSE)
  
  raw_csv_path <- file.path(dir_input_raw, paste0(sample_name, "_entropy.csv"))
  if (file.exists(raw_csv_path)) {
    df_raw <- read.csv(raw_csv_path, stringsAsFactors = FALSE)
    df <- merge(
      df_alt,
      df_raw[, c("spot_id", "entropy")],
      by = "spot_id",
      all.x = TRUE
    )
    names(df)[names(df) == "entropy"] <- "plugin_entropy"
  } else {
    df <- df_alt
    df$plugin_entropy <- NA_real_
  }
  
  valid_df <- df[!is.na(df$chao_shen_entropy) & !is.na(df$cwj_entropy) & df$nCount_RNA > 0, ]
  n_spots  <- nrow(valid_df)
  
  message("-> Spots with valid alternative entropies: ", n_spots)
  
  cor_cs_ncount <- cor.test(valid_df$nCount_RNA, valid_df$chao_shen_entropy, method = "spearman", exact = FALSE)
  cor_cs_nfeat  <- cor.test(valid_df$nFeature_RNA, valid_df$chao_shen_entropy, method = "spearman", exact = FALSE)
  
  cor_cwj_ncount <- cor.test(valid_df$nCount_RNA, valid_df$cwj_entropy, method = "spearman", exact = FALSE)
  cor_cwj_nfeat  <- cor.test(valid_df$nFeature_RNA, valid_df$cwj_entropy, method = "spearman", exact = FALSE)
  
  if (!all(is.na(valid_df$plugin_entropy))) {
    cor_pi_ncount <- cor.test(valid_df$nCount_RNA, valid_df$plugin_entropy, method = "spearman", exact = FALSE)
    cor_pi_nfeat  <- cor.test(valid_df$nFeature_RNA, valid_df$plugin_entropy, method = "spearman", exact = FALSE)
    rho_pi_ncount <- unname(cor_pi_ncount$estimate)
    rho_pi_nfeat  <- unname(cor_pi_nfeat$estimate)
  } else {
    rho_pi_ncount <- NA_real_
    rho_pi_nfeat  <- NA_real_
  }
  
  message(sprintf("-> Spearman rho with nCount_RNA:   Plug-in = %.3f | Chao-Shen = %.3f | CWJ = %.3f",
                  rho_pi_ncount, cor_cs_ncount$estimate, cor_cwj_ncount$estimate))
  message(sprintf("-> Spearman rho with nFeature_RNA: Plug-in = %.3f | Chao-Shen = %.3f | CWJ = %.3f",
                  rho_pi_nfeat, cor_cs_nfeat$estimate, cor_cwj_nfeat$estimate))
  
  summary_list[[sample_name]] <- data.frame(
    sample                  = sample_name,
    n_spots                 = n_spots,
    mean_chao_shen          = mean(valid_df$chao_shen_entropy),
    mean_cwj                = mean(valid_df$cwj_entropy),
    mean_nCount             = mean(valid_df$nCount_RNA),
    mean_nFeature           = mean(valid_df$nFeature_RNA),
    spearman_rho_plugin_nCount  = rho_pi_ncount,
    spearman_rho_cs_nCount      = unname(cor_cs_ncount$estimate),
    spearman_p_cs_nCount        = cor_cs_ncount$p.value,
    spearman_rho_cwj_nCount     = unname(cor_cwj_ncount$estimate),
    spearman_p_cwj_nCount       = cor_cwj_ncount$p.value,
    spearman_rho_plugin_nFeature = rho_pi_nfeat,
    spearman_rho_cs_nFeature    = unname(cor_cs_nfeat$estimate),
    spearman_p_cs_nFeature      = cor_cs_nfeat$p.value,
    spearman_rho_cwj_nFeature   = unname(cor_cwj_nfeat$estimate),
    spearman_p_cwj_nFeature     = cor_cwj_nfeat$p.value,
    stringsAsFactors        = FALSE
  )
  
  p_cs_ncount <- ggplot(valid_df, aes(x = nCount_RNA, y = chao_shen_entropy)) +
    geom_point(alpha = 0.25, size = 0.8, color = "#d95f02") +
    geom_smooth(method = "lm", formula = y ~ log10(x), color = "#7570b3", se = FALSE, linewidth = 1) +
    labs(
      title    = paste0(sample_name, ": Chao-Shen vs Total Counts"),
      subtitle = sprintf("Spearman rho = %.3f (vs Plug-in rho = %.3f)",
                         cor_cs_ncount$estimate, rho_pi_ncount),
      x        = "Total UMI Counts per Spot (nCount_RNA)",
      y        = "Chao-Shen Entropy (bits)"
    ) +
    theme_bw(base_size = 11)
  
  p_cs_nfeat <- ggplot(valid_df, aes(x = nFeature_RNA, y = chao_shen_entropy)) +
    geom_point(alpha = 0.25, size = 0.8, color = "#d95f02") +
    stat_function(
      fun      = log2,
      aes(color = "Observed-gene bound: log2(nFeature)"),
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_color_manual(name = NULL, values = c("Observed-gene bound: log2(nFeature)" = "#e7298a")) +
    labs(
      title    = paste0(sample_name, ": Chao-Shen vs Detected Genes"),
      subtitle = sprintf("Spearman rho = %.3f (vs Plug-in rho = %.3f)",
                         cor_cs_nfeat$estimate, rho_pi_nfeat),
      x        = "Detected Genes per Spot (nFeature_RNA)",
      y        = "Chao-Shen Entropy (bits)"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
  
  p_cwj_ncount <- ggplot(valid_df, aes(x = nCount_RNA, y = cwj_entropy)) +
    geom_point(alpha = 0.25, size = 0.8, color = "#1b9e77") +
    geom_smooth(method = "lm", formula = y ~ log10(x), color = "#7570b3", se = FALSE, linewidth = 1) +
    labs(
      title    = paste0(sample_name, ": Chao-Wang-Jost vs Total Counts"),
      subtitle = sprintf("Spearman rho = %.3f (vs Plug-in rho = %.3f)",
                         cor_cwj_ncount$estimate, rho_pi_ncount),
      x        = "Total UMI Counts per Spot (nCount_RNA)",
      y        = "CWJ Entropy (bits)"
    ) +
    theme_bw(base_size = 11)
  
  p_cwj_nfeat <- ggplot(valid_df, aes(x = nFeature_RNA, y = cwj_entropy)) +
    geom_point(alpha = 0.25, size = 0.8, color = "#1b9e77") +
    stat_function(
      fun      = log2,
      aes(color = "Observed-gene bound: log2(nFeature)"),
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_color_manual(name = NULL, values = c("Observed-gene bound: log2(nFeature)" = "#e7298a")) +
    labs(
      title    = paste0(sample_name, ": Chao-Wang-Jost vs Detected Genes"),
      subtitle = sprintf("Spearman rho = %.3f (vs Plug-in rho = %.3f)",
                         cor_cwj_nfeat$estimate, rho_pi_nfeat),
      x        = "Detected Genes per Spot (nFeature_RNA)",
      y        = "CWJ Entropy (bits)"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  
  grid_plot <- (p_cs_ncount | p_cs_nfeat) / (p_cwj_ncount | p_cwj_nfeat)
  
  fig_path <- file.path(dir_figures, paste0(sample_name, "_alternative_entropy_vs_depth.png"))
  ggsave(
    filename = fig_path,
    plot     = grid_plot,
    width    = 13,
    height   = 10,
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
    
    plot_df <- rbind(
      data.frame(
        sample    = summary_df$sample,
        estimator = "Plug-in (Baseline)",
        rho       = summary_df$spearman_rho_plugin_nCount
      ),
      data.frame(
        sample    = summary_df$sample,
        estimator = "Chao-Shen (2003)",
        rho       = summary_df$spearman_rho_cs_nCount
      ),
      data.frame(
        sample    = summary_df$sample,
        estimator = "Chao-Wang-Jost (2013)",
        rho       = summary_df$spearman_rho_cwj_nCount
      )
    )
    
    plot_df <- plot_df[!is.na(plot_df$rho), ]
    plot_df$estimator <- factor(
      plot_df$estimator,
      levels = c("Plug-in (Baseline)", "Chao-Shen (2003)", "Chao-Wang-Jost (2013)")
    )
    
    p_cohort <- ggplot(plot_df, aes(x = sample, y = rho, fill = estimator)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
      scale_fill_manual(values = c(
        "Plug-in (Baseline)"     = "#2b5c8f",
        "Chao-Shen (2003)"       = "#d95f02",
        "Chao-Wang-Jost (2013)"  = "#1b9e77"
      )) +
      labs(
        title    = "Cohort Overview: Spearman Correlation with Total Counts (nCount_RNA)",
        subtitle = "Comparing Plug-in, Chao-Shen, and Chao-Wang-Jost entropy estimators across all samples",
        x        = "Sample",
        y        = "Spearman Rank Correlation (rho)",
        fill     = "Entropy Estimator"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "top"
      )
    
    cohort_fig_path <- file.path(dir_figures, "cohort_alternative_correlations_overview.png")
    ggsave(
      filename = cohort_fig_path,
      plot     = p_cohort,
      width    = 13,
      height   = 6,
      dpi      = 300
    )
    message("-> Saved cohort overview plot to: ", cohort_fig_path)
  }
}

