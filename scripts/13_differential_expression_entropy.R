#SAMPLES_TO_RUN <- "sample1"
SAMPLES_TO_RUN <- "all"

# Must match the GROUP_PERCENTAGE used in Script 12
GROUP_PERCENTAGE <- 15

MIN_PCT         <- 0.10
LOGFC_THRESHOLD <- 0.25
PADJ_THRESHOLD  <- 0.05
TOP_SHARED_PLOT <- 25

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(here)
})

if (!is.numeric(GROUP_PERCENTAGE) || length(GROUP_PERCENTAGE) != 1 ||
    GROUP_PERCENTAGE <= 0 || GROUP_PERCENTAGE >= 50) {
  stop("GROUP_PERCENTAGE must be a single number strictly between 0 and 50 (e.g., 20).")
}

pct_folder <- paste0("top_bottom_", GROUP_PERCENTAGE, "pct")

dir_input_objects <- here::here("results", "05_spanorm_normalization", "objects")
dir_input_groups  <- here::here("results", "12_classify_entropy_groups", pct_folder, "tables")

if (!dir.exists(dir_input_groups)) {
  stop("Input directory from Script 12 not found: ", dir_input_groups,
       "\nPlease run Script 12 with GROUP_PERCENTAGE = ", GROUP_PERCENTAGE, " first.")
}

dir_results <- here::here("results", "13_differential_expression_entropy", pct_folder)
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

all_group_files <- list.files(dir_input_groups, pattern = "_entropy_groups\\.csv$", full.names = FALSE)
all_samples     <- gsub("_entropy_groups\\.csv$", "", all_group_files)

sample_numbers <- suppressWarnings(as.numeric(gsub("[^0-9]", "", all_samples)))
if (!any(is.na(sample_numbers))) {
  all_samples <- all_samples[order(sample_numbers)]
}

if (identical(SAMPLES_TO_RUN, "all")) {
  samples_to_process <- all_samples
} else {
  missing_samples <- setdiff(SAMPLES_TO_RUN, all_samples)
  if (length(missing_samples) > 0) {
    stop("The following requested samples were not found in ", dir_input_groups, ": ",
         paste(missing_samples, collapse = ", "))
  }
  samples_to_process <- SAMPLES_TO_RUN
}

message("Quantile group folder: ", pct_folder)
message("DE parameters: min.pct = ", MIN_PCT,
        " | logfc.threshold = ", LOGFC_THRESHOLD,
        " | p_val_adj < ", PADJ_THRESHOLD)
message("Samples selected (", length(samples_to_process), "): ",
        paste(samples_to_process, collapse = ", "))

sample_summary_list <- list()
all_sig_de_list     <- list()


for (sample_name in samples_to_process) {
  message("\nProcessing: ", sample_name)
  
  rds_file   <- file.path(dir_input_objects, paste0(sample_name, "_seurat.rds"))
  group_file <- file.path(dir_input_groups,  paste0(sample_name, "_entropy_groups.csv"))
  
  if (!file.exists(rds_file) || !file.exists(group_file)) {
    warning("Missing Seurat object or group file for ", sample_name, ". Skipping...")
    next
  }
  
  seurat_obj <- readRDS(rds_file)
  groups_df  <- read.csv(group_file, stringsAsFactors = FALSE)
  
  spot_ids <- colnames(seurat_obj)
  matched_groups <- groups_df$entropy_group[match(spot_ids, groups_df$spot_id)]
  
  seurat_obj$entropy_group <- matched_groups
  Seurat::Idents(seurat_obj) <- "entropy_group"
  
  n_high <- sum(matched_groups == "High", na.rm = TRUE)
  n_low  <- sum(matched_groups == "Low",  na.rm = TRUE)
  
  if (n_high < 3 || n_low < 3) {
    warning("Too few High or Low spots in ", sample_name, ". Skipping...")
    rm(seurat_obj, groups_df)
    next
  }
  
  message("-> Running Wilcoxon DE test (High = ", n_high, " spots vs. Low = ", n_low, " spots)...")
  
  de_res <- Seurat::FindMarkers(
    object          = seurat_obj,
    ident.1         = "High",
    ident.2         = "Low",
    assay           = "Spatial",
    slot            = "data",
    test.use        = "wilcox",
    min.pct         = MIN_PCT,
    logfc.threshold = LOGFC_THRESHOLD,
    verbose         = FALSE
  )
  
  if (nrow(de_res) == 0) {
    warning("No genes passed the pre-filtering thresholds in ", sample_name, ".")
    rm(seurat_obj, groups_df)
    next
  }
  
  de_res$gene   <- rownames(de_res)
  de_res$sample <- sample_name
  rownames(de_res) <- NULL
  
  de_res$direction <- "Not_Significant"
  de_res$direction[de_res$p_val_adj < PADJ_THRESHOLD & de_res$avg_log2FC >=  LOGFC_THRESHOLD] <- "Up_in_High_Entropy"
  de_res$direction[de_res$p_val_adj < PADJ_THRESHOLD & de_res$avg_log2FC <= -LOGFC_THRESHOLD] <- "Up_in_Low_Entropy"
  
  col_order <- c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj", "direction", "sample")
  de_res    <- de_res[, intersect(col_order, colnames(de_res))]
  
  de_res <- de_res[order(de_res$p_val_adj, -abs(de_res$avg_log2FC)), ]
  
  n_up_high <- sum(de_res$direction == "Up_in_High_Entropy")
  n_up_low  <- sum(de_res$direction == "Up_in_Low_Entropy")
  
  message("-> Significant DE genes (p_val_adj < ", PADJ_THRESHOLD, "): ",
          n_up_high, " Up in High Entropy | ",
          n_up_low,  " Up in Low Entropy")
  
  out_de_csv <- file.path(dir_tables, paste0(sample_name, "_de_genes.csv"))
  write.csv(de_res, file = out_de_csv, row.names = FALSE)
  
  sample_summary_list[[sample_name]] <- data.frame(
    sample            = sample_name,
    group_percentage  = GROUP_PERCENTAGE,
    n_high_spots      = n_high,
    n_low_spots       = n_low,
    n_tested_genes    = nrow(de_res),
    n_up_high_entropy = n_up_high,
    n_up_low_entropy  = n_up_low,
    stringsAsFactors  = FALSE
  )
  
  sig_df <- de_res[de_res$direction != "Not_Significant", ]
  if (nrow(sig_df) > 0) {
    all_sig_de_list[[sample_name]] <- sig_df
  }
  
  plot_df <- de_res
  min_nonzero_p <- min(plot_df$p_val_adj[plot_df$p_val_adj > 0], na.rm = TRUE)
  if (!is.finite(min_nonzero_p)) min_nonzero_p <- 1e-300
  plot_df$neg_log10_padj <- -log10(pmax(plot_df$p_val_adj, min_nonzero_p * 0.1))
  
  top_high <- head(plot_df[plot_df$direction == "Up_in_High_Entropy", ], 10)
  top_low  <- head(plot_df[plot_df$direction == "Up_in_Low_Entropy",  ], 10)
  label_df <- rbind(top_high, top_low)
  
  volcano_colors <- c(
    "Up_in_High_Entropy" = "#d7191c",
    "Up_in_Low_Entropy"  = "#2c7bb6",
    "Not_Significant"    = "grey70"
  )
  
  p_volcano <- ggplot(plot_df, aes(x = avg_log2FC, y = neg_log10_padj, color = direction)) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_vline(xintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(PADJ_THRESHOLD), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = volcano_colors) +
    labs(
      title    = paste0(sample_name, ": High vs. Low Entropy Differential Expression"),
      subtitle = paste0("Top/Bottom ", GROUP_PERCENTAGE, "% SpaNorm Entropy | ",
                        "Up in High: ", n_up_high, " | Up in Low: ", n_up_low),
      x        = "Average Log2 Fold Change (High / Low Entropy)",
      y        = "-log10(Adjusted P-value)",
      color    = "Classification"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      legend.position = "top"
    )
  
  if (nrow(label_df) > 0) {
    p_volcano <- p_volcano +
      ggrepel::geom_text_repel(
        data          = label_df,
        aes(label     = gene),
        size          = 3.2,
        max.overlaps  = 20,
        box.padding   = 0.35,
        show.legend   = FALSE
      )
  }
  
  volcano_file <- file.path(dir_figures, paste0(sample_name, "_volcano.png"))
  ggsave(filename = volcano_file, plot = p_volcano, width = 8.5, height = 6.5, dpi = 300)
  
  rm(seurat_obj, groups_df, de_res, plot_df, p_volcano)
  gc(verbose = FALSE)
  
  message("-> Done with ", sample_name, "!")
}


if (length(sample_summary_list) > 0) {
  summary_df <- do.call(rbind, sample_summary_list)
  rownames(summary_df) <- NULL
  
  summary_csv <- file.path(dir_tables, "de_summary_per_sample.csv")
  write.csv(summary_df, file = summary_csv, row.names = FALSE)
  message("\nSaved per-sample DE summary to: ", summary_csv)
}

if (length(all_sig_de_list) > 0) {
  all_sig_df <- do.call(rbind, all_sig_de_list)
  rownames(all_sig_df) <- NULL
  n_total_samples <- length(sample_summary_list)
  
  summarize_shared_genes <- function(df_sub, total_samples) {
    if (nrow(df_sub) == 0) return(data.frame())
    
    genes_unique <- unique(df_sub$gene)
    rows_out     <- vector("list", length(genes_unique))
    
    for (i in seq_along(genes_unique)) {
      g     <- genes_unique[i]
      g_sub <- df_sub[df_sub$gene == g, ]
      
      rows_out[[i]] <- data.frame(
        gene          = g,
        n_samples     = nrow(g_sub),
        prop_samples  = nrow(g_sub) / total_samples,
        mean_log2FC   = mean(g_sub$avg_log2FC),
        median_log2FC = median(g_sub$avg_log2FC),
        mean_pct_1    = mean(g_sub$pct.1),
        mean_pct_2    = mean(g_sub$pct.2),
        samples       = paste(g_sub$sample, collapse = ", "),
        stringsAsFactors = FALSE
      )
    }
    
    res_df <- do.call(rbind, rows_out)
    res_df <- res_df[order(-res_df$n_samples, -abs(res_df$mean_log2FC)), ]
    rownames(res_df) <- NULL
    return(res_df)
  }
  
  high_sub <- all_sig_df[all_sig_df$direction == "Up_in_High_Entropy", ]
  low_sub  <- all_sig_df[all_sig_df$direction == "Up_in_Low_Entropy",  ]
  
  shared_high_df <- summarize_shared_genes(high_sub, n_total_samples)
  shared_low_df  <- summarize_shared_genes(low_sub,  n_total_samples)
  
  high_csv <- file.path(dir_tables, "shared_high_entropy_genes.csv")
  low_csv  <- file.path(dir_tables, "shared_low_entropy_genes.csv")
  
  write.csv(shared_high_df, file = high_csv, row.names = FALSE)
  write.csv(shared_low_df,  file = low_csv,  row.names = FALSE)
  
  message("Saved shared High-Entropy genes table to: ", high_csv)
  message("Saved shared Low-Entropy genes table to:  ", low_csv)
  
  
  if (nrow(shared_high_df) > 0 || nrow(shared_low_df) > 0) {
    
    p_high <- NULL
    p_low  <- NULL
    
    if (nrow(shared_high_df) > 0) {
      top_h <- head(shared_high_df, TOP_SHARED_PLOT)
      top_h$gene <- factor(top_h$gene, levels = rev(top_h$gene))
      
      p_high <- ggplot(top_h, aes(x = gene, y = n_samples, fill = mean_log2FC)) +
        geom_col(color = "black", width = 0.75) +
        coord_flip() +
        scale_fill_gradient(low = "#fcae91", high = "#cb181d") +
        scale_y_continuous(breaks = seq(0, n_total_samples, by = max(1, floor(n_total_samples / 5))),
                           limits = c(0, n_total_samples)) +
        labs(
          title    = paste0("Top Shared High-Entropy Genes"),
          x        = "Gene",
          y        = "Number of Samples Significant",
          fill     = "Mean log2FC"
        ) +
        theme_bw(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 12))
    }
    
    if (nrow(shared_low_df) > 0) {
      top_l <- head(shared_low_df, TOP_SHARED_PLOT)
      top_l$gene <- factor(top_l$gene, levels = rev(top_l$gene))
      
      p_low <- ggplot(top_l, aes(x = gene, y = n_samples, fill = abs(mean_log2FC))) +
        geom_col(color = "black", width = 0.75) +
        coord_flip() +
        scale_fill_gradient(low = "#bdd7e7", high = "#08519c") +
        scale_y_continuous(breaks = seq(0, n_total_samples, by = max(1, floor(n_total_samples / 5))),
                           limits = c(0, n_total_samples)) +
        labs(
          title    = paste0("Top Shared Low-Entropy Genes"),
          x        = "Gene",
          y        = "Number of Samples Significant",
          fill     = "Mean log2FC"
        ) +
        theme_bw(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 12))
    }
    
    if (!is.null(p_high) && !is.null(p_low)) {
      overview_plot <- (p_high | p_low) +
        patchwork::plot_annotation(
          title = paste0("Most Frequent DE Genes:",
                         GROUP_PERCENTAGE, "%"),
          theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))
        )
      overview_width <- 14
    } else if (!is.null(p_high)) {
      overview_plot  <- p_high
      overview_width <- 7.5
    } else {
      overview_plot  <- p_low
      overview_width <- 7.5
    }
    
    overview_fig <- file.path(dir_figures, "shared_de_genes_overview.png")
    ggsave(filename = overview_fig, plot = overview_plot, width = overview_width, height = 8, dpi = 300)
    message("Saved shared DE genes overview figure to: ", overview_fig)
  }
} else {
  message("\nNo significant differentially expressed genes found across the processed samples.")
}

