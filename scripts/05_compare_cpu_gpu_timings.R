
suppressPackageStartupMessages({
  library(here)
})

dir_cpu_tables <- here::here("results", "05_spanorm_chunked_cpu", "tables")
dir_gpu_tables <- here::here("results", "05_spanorm_chunked_gpu", "tables")
dir_output_tables <- here::here("results", "05_spanorm_comparison", "tables")

dir.create(dir_output_tables, recursive = TRUE, showWarnings = FALSE)

load_timings <- function(tables_dir, backend_name) {
  csv_files <- sort(list.files(tables_dir, pattern = "^timings_.*\\.csv$", full.names = TRUE))
  if (length(csv_files) == 0) {
    stop("No timing CSV files found in: ", tables_dir)
  }

  timing_list <- lapply(csv_files, function(f) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    df[, c("sample", "elapsed_seconds")]
  })

  combined <- do.call(rbind, timing_list)
  # Keep the latest run for each sample if any sample was run multiple times
  latest <- combined[!duplicated(combined$sample, fromLast = TRUE), ]
  colnames(latest)[colnames(latest) == "elapsed_seconds"] <- paste0(backend_name, "_seconds")
  latest
}

cpu_df <- load_timings(dir_cpu_tables, "cpu")
gpu_df <- load_timings(dir_gpu_tables, "gpu")

# 4. Merge tables by sample
comparison <- merge(cpu_df, gpu_df, by = "sample")

# Order naturally by sample number (sample1, sample2, ..., sample22)
sample_numbers <- as.integer(sub("^sample", "", comparison$sample))
comparison <- comparison[order(sample_numbers), ]

# 5. Compute differences and metrics
# time_saved_seconds: positive value indicates GPU was faster than CPU
comparison$time_saved_seconds <- comparison$cpu_seconds - comparison$gpu_seconds
comparison$cpu_minutes <- round(comparison$cpu_seconds / 60, 2)
comparison$gpu_minutes <- round(comparison$gpu_seconds / 60, 2)
comparison$time_saved_minutes <- round(comparison$time_saved_seconds / 60, 2)
comparison$speedup <- round(comparison$cpu_seconds / comparison$gpu_seconds, 2)

# Round raw seconds to 1 decimal place for cleaner display
comparison$cpu_seconds <- round(comparison$cpu_seconds, 1)
comparison$gpu_seconds <- round(comparison$gpu_seconds, 1)
comparison$time_saved_seconds <- round(comparison$time_saved_seconds, 1)

# 6. Calculate total summary row
total_cpu_sec <- sum(comparison$cpu_seconds)
total_gpu_sec <- sum(comparison$gpu_seconds)
total_saved_sec <- total_cpu_sec - total_gpu_sec

total_cpu_min <- total_cpu_sec / 60
total_gpu_min <- total_gpu_sec / 60
total_saved_min <- total_saved_sec / 60
overall_speedup <- total_cpu_sec / total_gpu_sec

total_row <- data.frame(
  sample = "TOTAL",
  cpu_seconds = round(total_cpu_sec, 1),
  gpu_seconds = round(total_gpu_sec, 1),
  time_saved_seconds = round(total_saved_sec, 1),
  cpu_minutes = round(total_cpu_min, 2),
  gpu_minutes = round(total_gpu_min, 2),
  time_saved_minutes = round(total_saved_min, 2),
  speedup = round(overall_speedup, 2),
  stringsAsFactors = FALSE
)

final_table <- rbind(comparison, total_row)

# 7. Print to console
cat("\n=== SpaNorm Runtime Comparison: CPU vs GPU ===\n\n")
print(final_table, row.names = FALSE)

cat("\n=== Summary ===\n")
cat(sprintf("Total CPU runtime:    %.2f minutes (%.2f hours)\n", total_cpu_min, total_cpu_min / 60))
cat(sprintf("Total GPU runtime:    %.2f minutes (%.2f hours)\n", total_gpu_min, total_gpu_min / 60))
cat(sprintf("Total time saved:     %.2f minutes (%.2f hours)\n", total_saved_min, total_saved_min / 60))
cat(sprintf("Overall GPU speedup:  %.2fx\n\n", overall_speedup))

# 8. Save output CSV
output_csv <- file.path(dir_output_tables, "cpu_vs_gpu_timings.csv")
write.csv(final_table, file = output_csv, row.names = FALSE)
cat("Saved comparison table to:\n  ", output_csv, "\n\n")
