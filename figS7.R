setwd("/g/data/kr68/neysa/r_plotting/rcode_per_fig")

# ------------------------------------------------------------------------------
#                 Libraries
# ------------------------------------------------------------------------------
.libPaths(c("/g/data/kr68/andre/R_libs"))

library(data.table)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
#                 Configuration
# ------------------------------------------------------------------------------
base_output_dir <- "/g/data/kr68/andre/fshd_pipeline/subsampling_experiment"
sample_ids <- c("AS2603", "JOUB61166", "RJ1207", "GUAT0705")
sd_settings <- c("5sd", "10sd", "10sd_50draws")
combined_output_file <- file.path(base_output_dir, "combined_subsampling_summary.tsv")

# ------------------------------------------------------------------------------
#                 Limit of detection analysis
# ------------------------------------------------------------------------------
#pdf(file.path(base_output_dir, "Subsamplin_Experiment_Analysis.pdf"),height=8, width = 12)

# Initialize list for data.tables
all_summaries_list <- list()

for (sample_id in sample_ids) {
  for (sd_setting in sd_settings) {
    file_path <- file.path(
      base_output_dir,
      paste0(sample_id, "_subasampling_", sd_setting),
      paste0(sample_id, "_read_ids_and_lengths_draw_summary.tsv")
    )
    
    if (file.exists(file_path)) {
      current_dt <- fread(file_path, sep="\t")
      current_dt[, SampleID := sample_id]
      current_dt[, SD_Setting := sd_setting]
      all_summaries_list[[length(all_summaries_list) + 1]] <- current_dt
    } else {
      # This warning message is still useful for debugging if files are missing,
      # but will be silent if all files are found.
      # Consider keeping it commented out or logging to a separate file if needed.
      # cat(sprintf("Warning: File not found - %s. Skipping.\n", file_path)) 
    }
  }
}

if (length(all_summaries_list) > 0) {
  combined_summary_dt <- rbindlist(all_summaries_list, fill = TRUE)
} else {
  # This message is also useful if no files are processed at all.
  # Consider keeping it commented out or logging to a separate file.
  # cat("No summary files were found or processed. No combined output file created.\n")
}

plot_setting <- "10sd_50draws"

combined_summary_dt[, Multiplex_Factor := as.factor(Multiplex_Factor)]
combined_summary_dt[, Size_Gb := Size / 10^9]
hline_data <- unique(combined_summary_dt[, .(Base_Target_Size_Gb=Base_Target_Size/10^9)])

# ------------------------------------------------------------------------------
#                 Fig S7a
# ------------------------------------------------------------------------------
ggplot(combined_summary_dt[SD_Setting==plot_setting], aes(x = Multiplex_Factor, y = Size_Gb)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.25, size = 2) +
  scale_y_continuous(breaks = seq(0,15,5))+
  geom_hline(data = hline_data, aes(yintercept = Base_Target_Size_Gb), color = "black", linetype = "dashed") +
  facet_grid(. ~ SampleID) +
  labs(
    title = "Library Size by Multiplex Factor per Sample",
    x = "Multiplex Factor",
    y = "Library Size (Gb)"
  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold")
  )

# ------------------------------------------------------------------------------
#                 Fig S7b
# ------------------------------------------------------------------------------
ggplot(combined_summary_dt[SD_Setting==plot_setting], aes(x = Size_Gb, y = Complete_4qA_Reads, color = Multiplex_Factor)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.3, size = 3) +
  geom_smooth(method="lm", aes(group = 1)) +
  facet_grid(. ~ SampleID) +
  labs(
    title = "Detection of 4qA reads by Library Size per Sample",
    x = "Library Size (Gb)",
    y = "Number of Complete 4qA reads",
    color = "Multiplex Factor"
  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
