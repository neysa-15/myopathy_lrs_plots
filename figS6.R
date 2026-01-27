setwd("/path/to/working/directory")

# ------------------------------------------------------------------------------
#                 Libraries
# ------------------------------------------------------------------------------
.libPaths(c("/path/to/your/r_library"))

library(data.table)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
#                 Configuration
# ------------------------------------------------------------------------------
base_dir <- "/path/to/d4z4ling/results_directory"

fshd1 <- c("FSHD1 sample list")
fshd1_biallelic <- c("Biallelic FSHD1 sample list")
fshd2 <- c("FSHD2 sample list")
fshd1_n_2 <- c("FSHD1+2 sample list")

# Read age of onset
age_onset <- "/path/to/age_onset_group.tsv"
age_onset <- fread(age_onset, sep = "\t", header = TRUE)

# ------------------------------------------------------------------------------
#                 Data input and processing
# ------------------------------------------------------------------------------

# Retrieve relevant files
all_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

summary_files <- unlist(lapply(all_dirs, function(dir) {
  list.files(path = dir, pattern = "mapped_features_summary\\.tsv$", full.names = TRUE)
}))

combined_dt <- rbindlist(lapply(summary_files, function(f) {
  dt <- fread(f)
  sample_name <- sub("_mapped_features_summary\\.tsv$", "", basename(f))
  dt[, Sample := sample_name]
  return(dt)
}), fill = TRUE)


combined_dt[Sample %in% fshd1, Sample_Label := "FSHD1"]
combined_dt[Sample %in% fshd1_biallelic, Sample_Label := "Biallelic FSHD1"]
combined_dt[Sample %in% fshd2, Sample_Label := "FSHD2"]
combined_dt[Sample %in% fshd1_n_2, Sample_Label := "FSHD1+2"] 

# Define desired order of Sample_Label
label_order <- c("FSHD1", "FSHD2", "FSHD1+2", "Biallelic FSHD1")
combined_dt[, Sample_Label := factor(Sample_Label, levels = label_order)]

# Order Sample within each Sample_Label
combined_dt[, Sample := factor(Sample, levels = combined_dt[order(Sample_Label, Sample), unique(Sample)])]

# Filter out samples that's not in any group
combined_dt <- combined_dt[!is.na(Sample_Label)]

combined_df <- as.data.frame(combined_dt)

# ------------------------------------------------------------------------------
# Methylation levels vs Age of Onset
# ------------------------------------------------------------------------------

# group samples and calculate mean
methylation_summary <- combined_df %>%
  filter(ReadLabel == "Complete 4qA") %>%
  group_by(Sample, floor(MappedEstimatedCopies)) %>%
  filter(n() > 1) %>%
  summarize(average_pLAM_methylation = mean(pLAM_Methylation_Percentage, na.rm = TRUE))

fshd2_methylation_summary <- combined_df %>%
  filter(Sample_Label == "FSHD2") %>%
  filter(ReadLabel == "Partial distal 4qA") %>%
  group_by(Sample) %>%
  summarize(average_pLAM_methylation = mean(pLAM_Methylation_Percentage, na.rm = TRUE))

# Join tables
all_methylation_summary <- bind_rows(methylation_summary, fshd2_methylation_summary)

combined_data <- age_onset %>%
  left_join(all_methylation_summary, by = "Sample")

#####################
# No FSHD2
non_fshd2 <- fshd1_df <- combined_data[!(combined_data$Group == "FSHD2")]
combined_data <- non_fshd2
#####################

#####################
# If only fshd1
# fshd1_df <- combined_data[combined_data$Group == "FSHD1"]
# combined_data <- fshd1_df
#####################

# Drop NA
combined_data <- combined_data[!is.na(`Age at onset`)]

# Plot correlation of methylation vs age
label_colours <- c(
  "FSHD1" = "#ffe0e0",
  "FSHD2" = "#e0f7fa",
  "FSHD1+2"= "#efe0ff",
  "FSHD1-biallelic" = "#fff3cd"
)

# ------------------------------------------------------------------------------
#             Sup Figure 6a - Distal methylation vs Age at onset
# ------------------------------------------------------------------------------
# Correlation test
correlation_res <- cor.test(combined_data$average_pLAM_methylation, combined_data$`Age at onset`)
r_value <- round(correlation_res$estimate, 3)
p_value <- signif(correlation_res$p.value, 3)
methyl_corr_label <- paste0("r = ", r_value, ", p = ", p_value)

# Plot
methyl_corr_plot <- ggplot(combined_data, aes(x = average_pLAM_methylation, y = `Age at onset`)) +
  geom_point(aes(color = Group), size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se=FALSE) +
  annotate("text", x = Inf, y = -Inf, label = methyl_corr_label, hjust = 1.1, vjust = -1.1, size = 5) +
  labs(
    x = "Distal methylation (%)",
    y = "Age of onset"
  ) +
  theme_bw(base_size = 14) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank()
  )

# ------------------------------------------------------------------------------
#             Sup Figure 6c - D4Z4 copies vs Distal Methylation
# ------------------------------------------------------------------------------
# Correlation between average pLAM methylation and d4z4 copies
corr_res_methyl_copies <- cor.test(combined_data$`floor(MappedEstimatedCopies)`, combined_data$average_pLAM_methylation)
r_value_methyl_copies <- round(corr_res_methyl_copies$estimate, 3)
p_value_methyl_copies <- signif(corr_res_methyl_copies$p.value, 3)
methyl_copies_corr_label <- paste0("r = ", r_value_methyl_copies, ", p = ", p_value_methyl_copies)

# Plot
copies_methyl_corr_plot <- ggplot(combined_data, aes(x = `floor(MappedEstimatedCopies)`, y = average_pLAM_methylation)) +
  geom_point(aes(color = Group), size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se=FALSE) +
  annotate("text", x = Inf, y = -Inf, label = methyl_copies_corr_label, hjust = 1.1, vjust = -1.1, size = 5) +
  labs(
    x = "D4Z4 Copies",
    y = "Distal methylation (%)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank()
  )


# ------------------------------------------------------------------------------
#             Sup Figure 6b - D4Z4 copies vs Age at onset
# ------------------------------------------------------------------------------
d4z4copies_summary <- combined_df %>%
  filter(ReadLabel == "Complete 4qA") %>%
  group_by(Sample, floor(MappedEstimatedCopies)) %>%
  filter(n() > 1) %>%
  summarize(d4z4_copies = mean(floor(MappedEstimatedCopies), na.rm = TRUE))

combined_d4z4_data <- d4z4copies_summary %>%
  left_join(age_onset, by = "Sample")

# Drop NA
# combined_d4z4_data <- combined_d4z4_data[!is.na(`Age at onset`)]

# Correlation test
d4z4_correlation_res <- cor.test(combined_d4z4_data$d4z4_copies, combined_d4z4_data$`Age at onset`)
d4z4_r_value <- round(d4z4_correlation_res$estimate, 3)
d4z4_p_value <- signif(d4z4_correlation_res$p.value, 3)
d4z4_corr_label <- paste0("r = ", d4z4_r_value, ", p = ", d4z4_p_value)

# Plot
copies_corr_plot <- ggplot(combined_d4z4_data, aes(x = d4z4_copies, y = `Age at onset`)) +
  geom_point(aes(color = Group), size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se=FALSE) +
  annotate("text", x = Inf, y = -Inf, label = d4z4_corr_label, hjust = 1.1, vjust = -1.1, size = 5) +
  labs(
    x = "D4Z4 Copies",
    y = "Age of onset"
  ) +
  theme_bw(base_size = 14) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank()
  )

combined_correlation_plot <- methyl_corr_plot + copies_corr_plot + copies_methyl_corr_plot + plot_layout(ncol = 3, widths = c(1, 1, 1))
ggsave("correlation_plot.pdf", plot = combined_correlation_plot, width = 72, height = 24, units = "cm")
