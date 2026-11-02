# 1. Dependencies and configuration

#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(patchwork)
})

dir_input_plugin  <- here::here("results", "01_load_and_compute_entropy", "tables")
dir_input_alt     <- here::here("results", "03_compute_alternative_entropies", "tables")
dir_input_spanorm <- here::here("results", "06_compute_spanorm_entropy", "tables")

dir_results <- here::here("results", "09_entropy_standard_deviation")
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

all_csv_files <- list.files(dir_input_plugin, pattern = "_entropy\\.csv$", full.names = FALSE)
all_samples   <- gsub("_entropy\\.csv$", "", all_csv_files)

sample_numbers <- suppressWarnings(as.numeric(gsub("[^0-9]", "", all_samples)))
if (!any(is.na(sample_numbers))) {
  all_samples <- all_samples[order(sample_numbers)]
}

if (identical(SAMPLES_TO_RUN, "all")) {
  samples_to_process <- all_samples
} else {
  samples_to_process <- SAMPLES_TO_RUN
}

message("=== Starting Analysis: 09_entropy_standard_deviation ===")
message("Samples selected (", length(samples_to_process), "): ",
        paste(samples_to_process, collapse = ", "))

### 2. Processing loop: Compute spot-level entropy SD per sample

summary_list <- list()

for (sample_name in samples_to_process) {
  
  message("\nProcessing: ", sample_name)
  
  # A. Load tables
  file_plugin  <- file.path(dir_input_plugin,  paste0(sample_name, "_entropy.csv"))
  file_alt     <- file.path(dir_input_alt,     paste0(sample_name, "_alternative_entropies.csv"))
  file_spanorm <- file.path(dir_input_spanorm, paste0(sample_name, "_spanorm_entropy.csv"))
  
  if (!file.exists(file_plugin) || !file.exists(file_alt) || !file.exists(file_spanorm)) {
    warning("Missing input files for ", sample_name, ". Skipping...")
    next
  }
  
  df_plugin  <- read.csv(file_plugin,  stringsAsFactors = FALSE)
  df_alt     <- read.csv(file_alt,     stringsAsFactors = FALSE)
  df_spanorm <- read.csv(file_spanorm, stringsAsFactors = FALSE)
  
  # B. Compute standard deviation
  sd_plugin    <- sd(df_plugin$entropy[is.finite(df_plugin$entropy)])
  sd_chao_shen <- sd(df_alt$chao_shen_entropy[is.finite(df_alt$chao_shen_entropy)])
  sd_cwj       <- sd(df_alt$cwj_entropy[is.finite(df_alt$cwj_entropy)])
  sd_spanorm   <- sd(df_spanorm$spanorm_entropy[is.finite(df_spanorm$spanorm_entropy)])
  
  summary_list[[sample_name]] <- data.frame(
    sample    = sample_name,
    plugin    = sd_plugin,
    chao_shen = sd_chao_shen,
    cwj       = sd_cwj,
    spanorm   = sd_spanorm,
    stringsAsFactors = FALSE
  )
  
  # C. Memory clean up
  rm(df_plugin, df_alt, df_spanorm)
  
  message("Done with ", sample_name)
}

if (length(summary_list) == 0) {
  stop("No sample standard deviations could be computed. Please check input files.")
}

sd_df <- do.call(rbind, summary_list)
rownames(sd_df) <- NULL

### 3. Save sample-level SD tables

# Save tidy per-sample numeric table
csv_per_sample <- file.path(dir_tables, "entropy_sd_per_sample.csv")
message("\nSaving per-sample SD table to: ", csv_per_sample)
write.csv(sd_df, file = csv_per_sample, row.names = FALSE)

# Build ranked side-by-side table
ord_plugin    <- order(-sd_df$plugin,    sd_df$sample, na.last = TRUE)
ord_chao_shen <- order(-sd_df$chao_shen, sd_df$sample, na.last = TRUE)
ord_cwj       <- order(-sd_df$cwj,       sd_df$sample, na.last = TRUE)
ord_spanorm   <- order(-sd_df$spanorm,   sd_df$sample, na.last = TRUE)

ranked_table <- data.frame(
  "plugin"    = paste0(sd_df$sample[ord_plugin],    ": ", sd_df$plugin[ord_plugin]),
  "chao-shen" = paste0(sd_df$sample[ord_chao_shen], ": ", sd_df$chao_shen[ord_chao_shen]),
  "cwj"       = paste0(sd_df$sample[ord_cwj],       ": ", sd_df$cwj[ord_cwj]),
  "spanorm"   = paste0(sd_df$sample[ord_spanorm],   ": ", sd_df$spanorm[ord_spanorm]),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

csv_ranked <- file.path(dir_tables, "entropy_sd_ranked.csv")
message("Saving ranked SD table to: ", csv_ranked)
write.csv(ranked_table, file = csv_ranked, row.names = FALSE)

### 4. Overall variation: Standard deviation of the standard deviations

method_keys   <- c("plugin", "chao_shen", "cwj", "spanorm")
method_labels <- c(
  "plugin"    = "Plug-in (Raw)",
  "chao_shen" = "Chao-Shen",
  "cwj"       = "Chao-Wang-Jost (CWJ)",
  "spanorm"   = "SpaNorm"
)

summary_rows <- list()

for (m in method_keys) {
  vals <- sd_df[[m]]
  vals <- vals[is.finite(vals)]
  
  summary_rows[[m]] <- data.frame(
    method_id = m,
    method    = method_labels[[m]],
    sd_of_sds = if (length(vals) >= 2) sd(vals) else NA_real_,
    mean_sd   = mean(vals),
    median_sd = median(vals),
    min_sd    = min(vals),
    max_sd    = max(vals),
    n_samples = length(vals),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL

# Order by standard deviation of standard deviations descending
if (all(!is.na(summary_df$sd_of_sds))) {
  summary_df <- summary_df[order(-summary_df$sd_of_sds), ]
}

csv_summary <- file.path(dir_tables, "entropy_sd_summary.csv")
message("Saving overall SD summary table to: ", csv_summary)
write.csv(summary_df, file = csv_summary, row.names = FALSE)

message("\n=== Summary of Entropy Standard Deviations Across Methods ===")
print(summary_df[, c("method", "sd_of_sds", "mean_sd", "median_sd", "min_sd", "max_sd")], row.names = FALSE)

### 5. Visualization: Graphical comparison of entropy variation

if (nrow(sd_df) >= 2) {
  message("\nGenerating comparison figure...")
  
  method_colors <- c(
    "Plug-in (Raw)"        = "#2b5c8f",
    "Chao-Shen"            = "#d95f02",
    "Chao-Wang-Jost (CWJ)" = "#1b9e77",
    "SpaNorm"              = "#7570b3"
  )
  
  # Reshape to long format for ggplot2
  sd_long <- data.frame(
    sample = rep(sd_df$sample, length(method_keys)),
    method = rep(unname(method_labels[method_keys]), each = nrow(sd_df)),
    sd     = c(sd_df$plugin, sd_df$chao_shen, sd_df$cwj, sd_df$spanorm),
    stringsAsFactors = FALSE
  )
  sd_long$method <- factor(sd_long$method, levels = unname(method_labels[method_keys]))
  
  # Panel A: Distribution of spot-level SDs per sample across methods
  p1 <- ggplot(sd_long, aes(x = method, y = sd, fill = method)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5, show.legend = FALSE) +
    geom_jitter(width = 0.15, size = 2.2, alpha = 0.8, shape = 21, color = "black", show.legend = FALSE) +
    scale_fill_manual(values = method_colors) +
    labs(
      title = "Distribution of Sample Entropy SDs",
      x = "Entropy Method",
      y = "Sample Entropy SD"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 15, hjust = 1)
    )
  
  # Panel B: Bar plot of the Standard Deviation of Standard Deviations
  # Order methods by sd_of_sds descending matching summary_df
  summary_df_plot <- summary_df
  summary_df_plot$method <- factor(summary_df_plot$method, levels = summary_df_plot$method)
  
  p2 <- ggplot(summary_df_plot, aes(x = method, y = sd_of_sds, fill = method)) +
    geom_col(alpha = 0.75, width = 0.5, color = "black", show.legend = FALSE) +
    geom_text(
      aes(label = sprintf("%.3f", sd_of_sds)),
      vjust = -0.4,
      size = 3.5,
      fontface = "bold"
    ) +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Standard Deviation of Standard Deviations",
      x = "Entropy Method",
      y = "SD of Standard Deviations"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 15, hjust = 1)
    )
  
  combined_plot <- p1 + p2 +
    plot_layout(ncol = 2) +
    plot_annotation(
      title = "Entropy Variation Across Methods",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray30")
      )
    )
  
  fig_file <- file.path(dir_figures, "entropy_sd_comparison.png")
  message("Saving comparison figure to: ", fig_file)
  ggsave(fig_file, plot = combined_plot, width = 11, height = 5.2, dpi = 300)
} else {
  message("\nSkipping figure generation: at least 2 samples are required to compute SD of SDs.")
}

message("\nAnalysis complete.")

