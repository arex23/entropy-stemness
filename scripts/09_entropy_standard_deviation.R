#1. Dependencies and configuration

#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

library(here)

dir_input_plugin  <- here::here("results", "01_load_and_compute_entropy", "tables")
dir_input_alt     <- here::here("results", "03_compute_alternative_entropies", "tables")
dir_input_spanorm <- here::here("results", "06_compute_spanorm_entropy", "tables")

dir_results <- here::here("results", "09_entropy_standard_deviation")
dir_tables  <- file.path(dir_results, "tables")

dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)

all_csv_files <- list.files(dir_input_plugin, pattern = "_entropy\\.csv$", full.names = FALSE)
all_samples   <- gsub("_entropy\\.csv$", "", all_csv_files)

sample_numbers <- as.numeric(gsub("sample", "", all_samples))
all_samples    <- all_samples[order(sample_numbers)]

if (identical(SAMPLES_TO_RUN, "all")) {
  samples_to_process <- all_samples
} else {
  samples_to_process <- SAMPLES_TO_RUN
}

### 2. Processing loop

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
  
  # C. Memroy clean up
  rm(df_plugin, df_alt, df_spanorm)
  
  message("Done with ", sample_name)
}

# 3. Ranking and saving

sd_df <- data.frame()

for (df in summary_list) {
  sd_df <- rbind(sd_df, df)
}
rownames(sd_df) <- NULL

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

csv_file <- file.path(dir_tables, "entropy_sd_ranked.csv")
message("Saving ranked SD table to: ", csv_file)
write.csv(ranked_table, file = csv_file, row.names = FALSE)

message("Analysis complete.")
