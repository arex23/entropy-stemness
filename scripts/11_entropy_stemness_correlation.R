# 1. Configuration and dependencies

#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(here)
})

dir_tables_01 <- here::here("results", "01_load_and_compute_entropy", "tables")
dir_tables_03 <- here::here("results", "03_compute_alternative_entropies", "tables")
dir_tables_06 <- here::here("results", "06_compute_spanorm_entropy", "tables")
dir_tables_10 <- here::here("results", "10_compute_module_score", "tables")

dir_results <- here::here("results", "11_entropy_stemness_correlation")
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

all_csv_10  <- list.files(dir_tables_10, pattern = "_module_score\\.csv$", full.names = FALSE)
all_samples <- gsub("_module_score\\.csv$", "", all_csv_10)

sample_numbers <- suppressWarnings(as.numeric(gsub("[^0-9]", "", all_samples)))
if (!any(is.na(sample_numbers))) {
  all_samples <- all_samples[order(sample_numbers)]
}

if (identical(SAMPLES_TO_RUN, "all")) {
  samples_to_process <- all_samples
} else {
  missing_samples <- setdiff(SAMPLES_TO_RUN, all_samples)
  if (length(missing_samples) > 0) {
    stop("The following requested samples were not found in ", dir_tables_10, ": ",
         paste(missing_samples, collapse = ", "))
  }
  samples_to_process <- SAMPLES_TO_RUN
}

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

# 2. Per-sample processing loop

for (sample_name in samples_to_process) {
  message("Analyzing: ", sample_name)
  
  csv_01_path <- file.path(dir_tables_01, paste0(sample_name, "_entropy.csv"))
  csv_03_path <- file.path(dir_tables_03, paste0(sample_name, "_alternative_entropies.csv"))
  csv_06_path <- file.path(dir_tables_06, paste0(sample_name, "_spanorm_entropy.csv"))
  csv_10_path <- file.path(dir_tables_10, paste0(sample_name, "_module_score.csv"))
  
  if (!file.exists(csv_06_path)) {
    csv_06_path <- file.path(dir_tables_06, paste0(sample_name, "_entropy.csv"))
  }
  
  missing_inputs <- c(
    if (!file.exists(csv_01_path)) "Script 01 Plug-in CSV",
    if (!file.exists(csv_03_path)) "Script 03 Chao-Shen/CWJ CSV",
    if (!file.exists(csv_06_path)) "Script 06 SpaNorm CSV",
    if (!file.exists(csv_10_path)) "Script 10 Module Score CSV"
  )
  
  if (length(missing_inputs) > 0) {
    warning("Skipping ", sample_name, " because required inputs are missing: ",
            paste(missing_inputs, collapse = ", "))
    next
  }
  
  df_01 <- read.csv(csv_01_path, stringsAsFactors = FALSE)
  df_03 <- read.csv(csv_03_path, stringsAsFactors = FALSE)
  df_06 <- read.csv(csv_06_path, stringsAsFactors = FALSE)
  df_10 <- read.csv(csv_10_path, stringsAsFactors = FALSE)
  
  spanorm_col <- if ("spanorm_entropy" %in% colnames(df_06)) "spanorm_entropy" else "entropy"
  score_col   <- if ("module_score" %in% colnames(df_10)) "module_score" else "score_2024_paper"
  
  merged_df <- data.frame(
    spot_id           = df_10$spot_id,
    nCount_RNA        = df_10$nCount_RNA,
    nFeature_RNA      = df_10$nFeature_RNA,
    module_score      = df_10[[score_col]],
    plugin_entropy    = df_01$entropy[match(df_10$spot_id, df_01$spot_id)],
    chao_shen_entropy = df_03$chao_shen_entropy[match(df_10$spot_id, df_03$spot_id)],
    cwj_entropy       = df_03$cwj_entropy[match(df_10$spot_id, df_03$spot_id)],
    spanorm_entropy   = df_06[[spanorm_col]][match(df_10$spot_id, df_06$spot_id)],
    stringsAsFactors  = FALSE
  )
  
  valid_mask <- merged_df$nCount_RNA > 0 &
    !is.na(merged_df$module_score) &
    !is.na(merged_df$plugin_entropy) &
    !is.na(merged_df$chao_shen_entropy) &
    !is.na(merged_df$cwj_entropy) &
    !is.na(merged_df$spanorm_entropy)
  
  valid_df <- merged_df[valid_mask, ]
  n_spots  <- nrow(valid_df)
  
  message("-> Spots with valid entropy and module score: ", n_spots)
  
  cor_pi  <- cor.test(valid_df$module_score, valid_df$plugin_entropy,    method = "spearman", exact = FALSE)
  cor_cs  <- cor.test(valid_df$module_score, valid_df$chao_shen_entropy, method = "spearman", exact = FALSE)
  cor_cwj <- cor.test(valid_df$module_score, valid_df$cwj_entropy,       method = "spearman", exact = FALSE)
  cor_spa <- cor.test(valid_df$module_score, valid_df$spanorm_entropy,   method = "spearman", exact = FALSE)
  
  rho_pi  <- unname(cor_pi$estimate)
  rho_cs  <- unname(cor_cs$estimate)
  rho_cwj <- unname(cor_cwj$estimate)
  rho_spa <- unname(cor_spa$estimate)
  
  message(sprintf(
    "-> Spearman rho (Module Score): Plug-in = %.3f | Chao-Shen = %.3f | CWJ = %.3f | SpaNorm = %.3f",
    rho_pi, rho_cs, rho_cwj, rho_spa
  ))
  
  summary_list[[sample_name]] <- data.frame(
    sample                 = sample_name,
    n_spots                = n_spots,
    mean_module_score      = mean(valid_df$module_score),
    sd_module_score        = sd(valid_df$module_score),
    mean_plugin_entropy    = mean(valid_df$plugin_entropy),
    mean_chao_shen_entropy = mean(valid_df$chao_shen_entropy),
    mean_cwj_entropy       = mean(valid_df$cwj_entropy),
    mean_spanorm_entropy   = mean(valid_df$spanorm_entropy),
    spearman_rho_plugin    = rho_pi,
    spearman_p_plugin      = cor_pi$p.value,
    spearman_rho_chao_shen = rho_cs,
    spearman_p_chao_shen   = cor_cs$p.value,
    spearman_rho_cwj       = rho_cwj,
    spearman_p_cwj         = cor_cwj$p.value,
    spearman_rho_spanorm   = rho_spa,
    spearman_p_spanorm     = cor_spa$p.value,
    stringsAsFactors       = FALSE
  )
  
  p_pi <- ggplot(valid_df, aes(x = module_score, y = plugin_entropy)) +
    geom_point(alpha = 0.25, size = 0.8, color = "#2b5c8f") +
    geom_smooth(method = "lm", formula = y ~ x, color = "#e7298a", se = FALSE, linewidth = 1) +
    labs(
      title    = paste0(sample_name, ": Plug-in vs Stemness"),
      subtitle = sprintf("Spearman rho = %.3f (p = %.2e)", rho_pi, cor_pi$p.value),
      x        = "Stemness Module Score",
      y        = "Plug-in Shannon Entropy (bits)"
    ) +
    theme_bw(base_size = 11)
  
  p_cs <- ggplot(valid_df, aes(x = module_score, y = chao_shen_entropy)) +
    geom_point(alpha = 0.25, size = 0.8, color = "#d95f02") +
    geom_smooth(method = "lm", formula = y ~ x, color = "#e7298a", se = FALSE, linewidth = 1) +
    labs(
      title    = paste0(sample_name, ": Chao-Shen vs Stemness"),
      subtitle = sprintf("Spearman rho = %.3f (p = %.2e)", rho_cs, cor_cs$p.value),
      x        = "Stemness Module Score",
      y        = "Chao-Shen Entropy (bits)"
    ) +
    theme_bw(base_size = 11)
  
  p_cwj <- ggplot(valid_df, aes(x = module_score, y = cwj_entropy)) +
    geom_point(alpha = 0.25, size = 0.8, color = "#1b9e77") +
    geom_smooth(method = "lm", formula = y ~ x, color = "#e7298a", se = FALSE, linewidth = 1) +
    labs(
      title    = paste0(sample_name, ": Chao-Wang-Jost vs Stemness"),
      subtitle = sprintf("Spearman rho = %.3f (p = %.2e)", rho_cwj, cor_cwj$p.value),
      x        = "Stemness Module Score",
      y        = "CWJ Entropy (bits)"
    ) +
    theme_bw(base_size = 11)
  
  p_spa <- ggplot(valid_df, aes(x = module_score, y = spanorm_entropy)) +
    geom_point(alpha = 0.25, size = 0.8, color = "#7570b3") +
    geom_smooth(method = "lm", formula = y ~ x, color = "#e7298a", se = FALSE, linewidth = 1) +
    labs(
      title    = paste0(sample_name, ": SpaNorm vs Stemness"),
      subtitle = sprintf("Spearman rho = %.3f (p = %.2e)", rho_spa, cor_spa$p.value),
      x        = "Stemness Module Score",
      y        = "SpaNorm Entropy (bits)"
    ) +
    theme_bw(base_size = 11)
  
  grid_plot <- ((p_pi | p_cs) / (p_cwj | p_spa)) +
    patchwork::plot_annotation(
      title    = paste0(sample_name, ": Shannon Entropy vs. Stemness Module Score (n = ", n_spots, " spots)"),
      theme    = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14)
      )
    )
  
  fig_path <- file.path(dir_figures, paste0(sample_name, "_entropy_vs_stemness.png"))
  ggsave(
    filename = fig_path,
    plot     = grid_plot,
    width    = 13,
    height   = 10,
    dpi      = 300
  )
  message("-> Saved diagnostic plot to: ", fig_path)
}

# 3. Cohort-wide summary tables and overview plots

if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  
  summary_csv_path <- file.path(dir_tables, "correlation_summary.csv")
  write.csv(summary_df, file = summary_csv_path, row.names = FALSE)
  message("\n-> Saved per-sample correlation summary table to: ", summary_csv_path)
  
  cohort_overview <- data.frame(
    method        = method_levels,
    mean_rho      = c(
      mean(summary_df$spearman_rho_plugin),
      mean(summary_df$spearman_rho_chao_shen),
      mean(summary_df$spearman_rho_cwj),
      mean(summary_df$spearman_rho_spanorm)
    ),
    sd_rho        = c(
      sd(summary_df$spearman_rho_plugin),
      sd(summary_df$spearman_rho_chao_shen),
      sd(summary_df$spearman_rho_cwj),
      sd(summary_df$spearman_rho_spanorm)
    ),
    median_rho    = c(
      median(summary_df$spearman_rho_plugin),
      median(summary_df$spearman_rho_chao_shen),
      median(summary_df$spearman_rho_cwj),
      median(summary_df$spearman_rho_spanorm)
    ),
    min_rho       = c(
      min(summary_df$spearman_rho_plugin),
      min(summary_df$spearman_rho_chao_shen),
      min(summary_df$spearman_rho_cwj),
      min(summary_df$spearman_rho_spanorm)
    ),
    max_rho       = c(
      max(summary_df$spearman_rho_plugin),
      max(summary_df$spearman_rho_chao_shen),
      max(summary_df$spearman_rho_cwj),
      max(summary_df$spearman_rho_spanorm)
    ),
    mean_abs_rho  = c(
      mean(abs(summary_df$spearman_rho_plugin)),
      mean(abs(summary_df$spearman_rho_chao_shen)),
      mean(abs(summary_df$spearman_rho_cwj)),
      mean(abs(summary_df$spearman_rho_spanorm))
    ),
    stringsAsFactors = FALSE
  )
  
  cohort_csv_path <- file.path(dir_tables, "cohort_stemness_correlation_overview.csv")
  write.csv(cohort_overview, file = cohort_csv_path, row.names = FALSE)
  message("-> Saved cohort-wide overview table to: ", cohort_csv_path)
  
  if (nrow(summary_df) > 1) {
    summary_df$sample <- factor(summary_df$sample, levels = summary_df$sample)
    
    plot_df <- rbind(
      data.frame(
        sample = summary_df$sample,
        method = "Plug-in",
        rho    = summary_df$spearman_rho_plugin
      ),
      data.frame(
        sample = summary_df$sample,
        method = "Chao-Shen",
        rho    = summary_df$spearman_rho_chao_shen
      ),
      data.frame(
        sample = summary_df$sample,
        method = "Chao-Wang-Jost",
        rho    = summary_df$spearman_rho_cwj
      ),
      data.frame(
        sample = summary_df$sample,
        method = "SpaNorm",
        rho    = summary_df$spearman_rho_spanorm
      )
    )
    
    plot_df$method <- factor(plot_df$method, levels = method_levels)
    
    p_cohort <- ggplot(plot_df, aes(x = sample, y = rho, fill = method)) +
      geom_col(position = position_dodge(width = 0.85), width = 0.78) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
      scale_fill_manual(values = method_colors) +
      labs(
        title    = "Spearman Correlation Between Entropy and Stemness Module Score",
        subtitle = "Evaluated across all processed samples for all four entropy estimators",
        x        = "Sample",
        y        = "Spearman Rank Correlation with Stemness Score (rho)",
        fill     = "Entropy Method"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "top"
      )
    
    cohort_fig_path <- file.path(dir_figures, "cohort_stemness_correlations_overview.png")
    ggsave(
      filename = cohort_fig_path,
      plot     = p_cohort,
      width    = 15,
      height   = 6.5,
      dpi      = 300
    )
    message("-> Saved cohort correlation overview bar plot to: ", cohort_fig_path)
    
    cohort_overview$method <- factor(cohort_overview$method, levels = method_levels)
    
    p_mean_bar <- ggplot(cohort_overview, aes(x = method, y = mean_rho, fill = method)) +
      geom_col(width = 0.6, alpha = 0.9) +
      geom_errorbar(
        aes(
          ymin = mean_rho - sd_rho,
          ymax = mean_rho + sd_rho
        ),
        width = 0.2,
        linewidth = 0.7
      ) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
      scale_fill_manual(values = method_colors) +
      labs(
        title    = "Mean Stemness Correlation Across Cohort",
        subtitle = "Average Spearman rho with Stemness Module Score",
        x        = NULL,
        y        = "Mean Spearman rho"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x     = element_text(angle = 20, hjust = 1),
        legend.position = "none"
      )
    
    p_box <- ggplot(plot_df, aes(x = method, y = rho, fill = method)) +
      geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.12, size = 2, alpha = 0.85, color = "black") +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.6, linetype = "dashed") +
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
      scale_fill_manual(values = method_colors) +
      labs(
        title    = "Distribution of Sample Stemness Correlations",
        subtitle = paste0("Each point represents one sample (n = ", nrow(summary_df), ")"),
        x        = NULL,
        y        = "Spearman rho with Stemness Module Score"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x     = element_text(angle = 20, hjust = 1),
        legend.position = "none"
      )
    
    p_summary_combined <- p_mean_bar | p_box
    
    summary_fig_path <- file.path(dir_figures, "cohort_stemness_correlation_summary_boxplot.png")
    ggsave(
      filename = summary_fig_path,
      plot     = p_summary_combined,
      width    = 13,
      height   = 6,
      dpi      = 300
    )
    message("-> Saved cohort correlation summary boxplot to: ", summary_fig_path)
  }
}

message("\n=== Analysis complete for script 11 ===")
