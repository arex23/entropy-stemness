#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(here)
})

dir_input_tables <- here::here("results", "01_load_and_compute_entropy", "tables")

dir_results <- here::here("results", "02_entropy_depth_correlation")
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

all_csv_files <- list.files(dir_input_tables, pattern = "_entropy\\.csv$", full.names = FALSE)
all_samples   <- gsub("_entropy\\.csv$", "", all_csv_files)

sample_numbers <- as.numeric(gsub("sample", "", all_samples))
all_samples    <- all_samples[order(sample_numbers)]

if (identical(SAMPLES_TO_RUN, "all")) {
  samples_to_process <- all_samples
} else {
  samples_to_process <- SAMPLES_TO_RUN
}

message("=== Starting Analysis: 02_entropy_depth_correlation ===")
message("Samples selected: ", paste(samples_to_process, collapse = ", "))

summary_list <- list()

for (sample_name in samples_to_process) {
  message("\n------------------------------------------------------------")
  message("Analyzing: ", sample_name)
  message("------------------------------------------------------------")
  
  csv_path <- file.path(dir_input_tables, paste0(sample_name, "_entropy.csv"))
  
  if (!file.exists(csv_path)) {
    warning("CSV file not found for ", sample_name, ": ", csv_path, ". Skipping...")
    next
  }
  
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  
  positions <- read.csv(here::here("data", "raw", sample_name, "spatial", "tissue_positions.csv"))
  in_tissue <- positions$in_tissue[match(df$spot_id, positions$barcode)]
  if (anyNA(in_tissue)) {
    stop("Missing tissue flags for spots in ", sample_name)
  }

  valid_df <- df[in_tissue == 1 & !is.na(df$entropy) & df$nCount_RNA > 0, ]
  n_spots  <- nrow(valid_df)
  
  message("-> Tissue spots with valid entropy: ", n_spots, " (out of ", nrow(df), " total spots)")
  
  cor_spearman_ncount <- cor.test(valid_df$nCount_RNA, valid_df$entropy, method = "spearman", exact = FALSE)
  
  cor_spearman_nfeat  <- cor.test(valid_df$nFeature_RNA, valid_df$entropy, method = "spearman", exact = FALSE)
  
  message(sprintf("-> Entropy vs nCount_RNA:   Spearman rho = %.3f (p = %.2e)",
                  cor_spearman_ncount$estimate, cor_spearman_ncount$p.value))
  message(sprintf("-> Entropy vs nFeature_RNA: Spearman rho = %.3f (p = %.2e)",
                  cor_spearman_nfeat$estimate, cor_spearman_nfeat$p.value))
  
  summary_list[[sample_name]] <- data.frame(
    sample                = sample_name,
    n_spots               = n_spots,
    mean_entropy          = mean(valid_df$entropy),
    sd_entropy            = sd(valid_df$entropy),
    mean_nCount           = mean(valid_df$nCount_RNA),
    mean_nFeature         = mean(valid_df$nFeature_RNA),
    spearman_rho_nCount   = unname(cor_spearman_ncount$estimate),
    spearman_p_nCount     = cor_spearman_ncount$p.value,
    spearman_rho_nFeature = unname(cor_spearman_nfeat$estimate),
    spearman_p_nFeature   = cor_spearman_nfeat$p.value,
    stringsAsFactors      = FALSE
  )
  
  p_ncount <- ggplot(valid_df, aes(x = nCount_RNA, y = entropy)) +
    geom_point(alpha = 0.25, size = 0.8, color = "#2b5c8f") +
    geom_smooth(method = "lm", formula = y ~ log10(x), color = "#d95f02", se = FALSE, linewidth = 1) +
    labs(
      title    = paste0(sample_name, ": Entropy vs Total Counts (nCount_RNA)"),
      subtitle = sprintf("Spearman rho = %.3f (n = %d spots)",
                         cor_spearman_ncount$estimate, n_spots),
      x        = "Total UMI Counts per Spot (nCount_RNA)",
      y        = "Shannon Entropy (bits)"
    ) +
    theme_bw(base_size = 12)
  
  p_nfeat <- ggplot(valid_df, aes(x = nFeature_RNA, y = entropy)) +
    geom_point(alpha = 0.25, size = 0.8, color = "#2b8f5c") +
    stat_function(
      fun      = log2,
      aes(color = "Theoretical Max: log2(nFeature)"),
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_color_manual(name = NULL, values = c("Theoretical Max: log2(nFeature)" = "#e7298a")) +
    labs(
      title    = paste0(sample_name, ": Entropy vs Detected Genes (nFeature_RNA)"),
      subtitle = sprintf("Spearman rho = %.3f",
                         cor_spearman_nfeat$estimate),
      x        = "Detected Genes per Spot (nFeature_RNA)",
      y        = "Shannon Entropy (bits)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  
  combined_fig <- p_ncount | p_nfeat
  
  fig_path <- file.path(dir_figures, paste0(sample_name, "_entropy_vs_depth.png"))
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
    
    plot_df <- rbind(
      data.frame(
        sample = summary_df$sample,
        metric = "nCount_RNA (Spearman rho)",
        value  = summary_df$spearman_rho_nCount
      ),
      data.frame(
        sample = summary_df$sample,
        metric = "nFeature_RNA (Spearman rho)",
        value  = summary_df$spearman_rho_nFeature
      )
    )
    
    p_cohort <- ggplot(plot_df, aes(x = sample, y = value, fill = metric)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_manual(values = c(
        "nCount_RNA (Spearman rho)"   = "#2b5c8f",
        "nFeature_RNA (Spearman rho)" = "#2b8f5c"
      )) +
      labs(
        title    = "Cohort Overview: Spearman Correlation with Plug-in Shannon Entropy",
        subtitle = "Evaluated across all processed samples",
        x        = "Sample",
        y        = "Spearman Rank Correlation (rho)",
        fill     = "Technical / Coverage Variable"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "top"
      )
    
    cohort_fig_path <- file.path(dir_figures, "cohort_correlations_overview.png")
    ggsave(
      filename = cohort_fig_path,
      plot     = p_cohort,
      width    = 12,
      height   = 6,
      dpi      = 300
    )
    message("-> Saved cohort overview plot to: ", cohort_fig_path)
  }
}
