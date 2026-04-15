setwd("/path/to/your/working/dir")

# ------------------------------------------------------------------------------
#                 Libraries
# ------------------------------------------------------------------------------
.libPaths(c("/path/to/your/R_libs"))

library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)

# ------------------------------------------------------------------------------
#                 Data loading and processing
# ------------------------------------------------------------------------------
cov_sum <- fread("/path/to/sample_coverage_summary.tsv", header = TRUE) 
setnames(cov_sum, "sample", "Sample")

# Load and clean sample key
sample_key <- fread("/path/to/sample_key.tsv", sep = "\t", header = FALSE)
names(sample_key) <- c("LRS_ID", "Sample")
sample_key[, LRS_ID := gsub("RS0*", "", LRS_ID)]

myopathy <- c(
  "list of samples ordered by however user like to see it on the figure"
)

# DATA FOR FIG S2A and S2C
# Merge annotations
cov_sum <- merge(cov_sum, sample_key, by = "Sample", all.x = TRUE)
cov_sum <- cov_sum[Sample %in% myopathy]

# Convert to Gbp
cov_sum[, total_bp_Gbp := `total_bp_(Mbp)` / 1000]
cov_sum[, on_target_Gbp := `on_target_bases_(Mbp)` / 1000]
cov_sum[, off_target_Gbp := pmax(total_bp_Gbp - on_target_Gbp, 0)]

# Set factor levels for consistent sample order
label_order <- c("FSHD1", "Undiagnosed FSHD1", "FSHD2", "Negative control")
cov_sum[, Sample_Label := factor(Sample, levels = myopathy)]
sample_order_dt <- unique(cov_sum[, .(Sample, Sample_Label)])
setorder(sample_order_dt, Sample_Label, Sample)
cov_sum[, Sample := factor(Sample, levels = sample_order_dt$Sample)]

# DATA FOR FIG S2B and SX - coverage per sample or per gene
# Path to myopathy panel
myop_panel <- fread("/path/to/myopathy_panel.bed", sep = "\t", header = FALSE)
names(myop_panel) <- c("chr", "start", "end", "gene")

cov_merged <- fread("/path/to/coverage_per_gene_per_sample.tsv", header = TRUE)

cov_merged <- merge(cov_merged, sample_key, by = "Sample", all.x = TRUE)

# DATA FOR FIG S2D
mt_cov_sum <- fread("/path/to/mt_dna_coverage.txt", header = FALSE)
setnames(mt_cov_sum, old = c("V1", "V2"), new = c("Sample", "Coverage"))

# ------------------------------------------------------------------------------
#        Sup Fig 2a - ON-TARGET and OFF-TARGET PLOTS FOR SIZE
# ------------------------------------------------------------------------------
# Set Sample_Label as factor with levels in myopathy order
cov_sum[, Sample_Label := factor(Sample, levels = myopathy)]

# Ensure sample order for plotting
sample_order_dt <- unique(cov_sum[, .(Sample, Sample_Label)])
setorder(sample_order_dt, Sample_Label, Sample)

cov_sum[, Sample := factor(Sample, levels = sample_order_dt$Sample)]

# Melt data for plotting
plot_data <- melt(
  cov_sum[, .(LRS_ID, Sample, Sample_Label, on_target_Gbp, off_target_Gbp)],
  id.vars = c("LRS_ID", "Sample", "Sample_Label"),
  variable.name = "category",
  value.name = "Gbp"
)
plot_data[category == "on_target_Gbp", category := "On-target"]
plot_data[category == "off_target_Gbp", category := "Off-target"]
plot_data$category <- factor(plot_data$category, levels = c("Off-target", "On-target"))

# Background fill colours
label_colours <- c(
  "FSHD1" = "#ffe0e0",
  "Undiagnosed FSHD1" = "#fff3cd",
  "FSHD2" = "#e0f7fa",
  "Negative control" = "#e0e0e0"
)

# # Set LRS_ID order based on diagnostic group
# cov_sum[, LRS_ID := factor(LRS_ID.y, levels = cov_sum[order(Sample_Label, Sample)]$LRS_ID.y)]
# plot_data[, LRS_ID := factor(LRS_ID, levels = levels(cov_sum$LRS_ID))]

# Reorder the rows of your key table to match the order in the 'myopathy' vector
ordered_key <- sample_key[match(myopathy, Sample), ]

# From this newly ordered table, extract the LRS_ID column.
# The !is.na() handles cases where a sample in 'myopathy' isn't in your key.
ordered_lrs_ids <- ordered_key[!is.na(LRS_ID)]$LRS_ID

# Apply the new LRS_ID order to your plot and summary data
cov_sum[, LRS_ID := factor(LRS_ID, levels = ordered_lrs_ids)]
plot_data[, LRS_ID := factor(LRS_ID, levels = ordered_lrs_ids)]

# Background shading data
bg_df <- unique(plot_data[, .(LRS_ID, Sample_Label)])
bg_df[, Sample_Num := as.numeric(LRS_ID)]
bg_df[, `:=`(
  xmin = Sample_Num - 0.5,
  xmax = Sample_Num + 0.5,
  ymin = -Inf,
  ymax = Inf,
  fill = label_colours[as.character(Sample_Label)]
)]

# Plot
ggplot(plot_data, aes(x = LRS_ID, y = Gbp)) +
  # draw background with hardcoded fill colour (not aes)
  # geom_rect(
  #   data = bg_df,
  #   inherit.aes = FALSE,
  #   aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Sample_Label),
  #   alpha = 0.8
  # )+
  # draw bars using category-based fill
  geom_bar(stat = "identity", aes(fill = category)) +
  scale_fill_manual(
    values = c(
      "On-target"        = "#3EB489",  # pastel green
      "Off-target"       = "#8C59A8",  # pastel purple
      "FSHD1"            = "#ffe0e0",
      "Undiagnosed FSHD1"= "#fff3cd",
      "FSHD2"            = "#e0f7fa",
      "Negative control" = "#e0e0e0"
    ),
    name = NULL,
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  labs(
    y = "Total bases (Gbp)",
    x = "Sample"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# ------------------------------------------------------------------------------
#      Sup Fig 2d - ON-TARGET and OFF-TARGET PLOTS FOR N50
# ------------------------------------------------------------------------------
# Melt data for plotting
plot_data <- melt(
  cov_sum[, .(LRS_ID, Sample, Sample_Label, n50_on_target, n50)],
  id.vars = c("LRS_ID", "Sample", "Sample_Label"),
  variable.name = "category",
  value.name = "n50"
)
plot_data[category == "n50_on_target", category := "On-target"]
plot_data[category == "n50", category := "Off-target"]
plot_data$category <- factor(plot_data$category, levels = c("Off-target", "On-target"))

# Background fill colours
label_colours <- c(
  "FSHD1" = "#ffe0e0",
  "Undiagnosed FSHD1" = "#fff3cd",
  "FSHD2" = "#e0f7fa",
  "Negative control" = "#e0e0e0"
)

# Set LRS_ID order based on diagnostic group
cov_sum[, LRS_ID := factor(LRS_ID, levels = cov_sum[order(Sample_Label, Sample)]$LRS_ID)]
plot_data[, LRS_ID := factor(LRS_ID, levels = levels(cov_sum$LRS_ID))]

# Background shading data
bg_df <- unique(plot_data[, .(LRS_ID, Sample_Label)])
bg_df[, Sample_Num := as.numeric(LRS_ID)]
bg_df[, `:=`(
  xmin = Sample_Num - 0.5,
  xmax = Sample_Num + 0.5,
  ymin = -Inf,
  ymax = Inf,
  fill = label_colours[as.character(Sample_Label)]
)]

# Plot
ggplot(plot_data, aes(x = LRS_ID, y = n50, color = category)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  scale_color_manual(
    values = c(
      "On-target"        = "#3EB489",  # pastel green
      "Off-target"       = "#8C59A8"   # pastel purple
    ),
    name = NULL,
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  labs(
    y = "n50",
    x = "Sample"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# ------------------------------------------------------------------------------
#      Sup Fig 2e - mt DNA coverage
# ------------------------------------------------------------------------------

mt_cov_sum <- merge(mt_cov_sum, sample_key, by = "Sample", all.x = TRUE)

# enforce order on Sample
mt_cov_sum$Sample <- factor(mt_cov_sum$Sample, levels = myopathy)

# get unique mapping of Sample → LRS_ID in myopathy order
id_order <- mt_cov_sum[!duplicated(mt_cov_sum$Sample), c("Sample", "LRS_ID")]
id_order <- id_order[order(factor(id_order$Sample, levels = myopathy)), ]

# set factor levels of LRS_ID based on that order
mt_cov_sum$LRS_ID <- factor(mt_cov_sum$LRS_ID, levels = id_order$LRS_ID)

mt_cov_sum_filter <- mt_cov_sum[!is.na(mt_cov_sum$Sample), ]

ggplot(mt_cov_sum_filter, aes(x = LRS_ID, y = Coverage)) +
  # geom_point(size = 4, position = position_dodge(width = 0.5), colour = "#ED872D") +
  geom_bar(stat = "identity", fill = "#ED872D") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    title = "Mitochondrial DNA Coverage distribution per sample",
    x = "Sample (LRS_ID)",
    y = "Coverage"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

mt_cov_df <- as.data.frame(mt_cov_sum_filter)
mean(mt_cov_df$Coverage)
sd(mt_cov_df$Coverage)
min(mt_cov_df$Coverage)
max(mt_cov_df$Coverage)


# ------------------------------------------------------------------------------
#      Sup Fig 2b - Coverage per target of each sample
# ------------------------------------------------------------------------------

# enforce order on Sample
cov_merged$Sample <- factor(cov_merged$Sample, levels = myopathy)

# get unique mapping of Sample → LRS_ID in myopathy order
id_order <- cov_merged[!duplicated(cov_merged$Sample), c("Sample", "LRS_ID")]
id_order <- id_order[order(factor(id_order$Sample, levels = myopathy)), ]

# set factor levels of LRS_ID based on that order
cov_merged$LRS_ID <- factor(cov_merged$LRS_ID, levels = id_order$LRS_ID)

cov_merged_filter <- cov_merged[!is.na(cov_merged$Sample), ]

# plot
ggplot(cov_merged_filter, aes(x = LRS_ID, y = Coverage)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    title = "Coverage distribution per sample",
    x = "Sample (LRS_ID)",
    y = "Coverage"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

# ------------------------------------------------------------------------------
#      Sup Fig 2c - Coverage per gene in all samples
# ------------------------------------------------------------------------------

gene_of_interest <- myop_panel[["gene"]]

cov_merged_filter <- cov_merged[!is.na(cov_merged$Sample), ]

# plot - all at one
plot_gene_cov_distribution <- function(cov_dt, plot_title) {
  ggplot(cov_dt, aes(x = Gene, y = Coverage)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(
      title = plot_title, #"Coverage distribution per gene"
      x = "Panel regions (gene)",
      y = "Coverage"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )
}

plot_gene_cov_distribution(cov_merged_filter, "Coverage distribution per gene")

#########################################
# Plot top 50 and bottom 50

# Filter top 50
sorted_cov_sum <- cov_merged_filter %>%
  group_by(Gene) %>%
  arrange(Coverage)

# plot_gene_cov_distribution(sorted_cov_sum, "Coverage distribution per gene sorted")

df_sorted_cov_sum <- as.data.frame(sorted_cov_sum)

# Group and median sort
df_median_myop_panel <- df_sorted_cov_sum %>%
  group_by(Gene) %>%
  summarise(median_gene_cov = median(Coverage)) %>%
  arrange(median_gene_cov)

##################################################################### ****
# Bottom 50 median cov
df_bottom50_median <- df_median_myop_panel %>%
  slice_head(n=40)
bottom_50_median_cov_gene <- df_bottom50_median[["Gene"]]

# Plot bottom 50 cov median
filtered_cov_sum_median_bottom50 <- df_sorted_cov_sum %>% 
  filter(Gene %in% bottom_50_median_cov_gene)
filtered_cov_sum_median_bottom50$Gene <- factor(filtered_cov_sum_median_bottom50$Gene, levels = bottom_50_median_cov_gene)
bottom50_plot <- plot_gene_cov_distribution(filtered_cov_sum_median_bottom50, "Coverage distribution per gene lowest 50 median")

##################################################################### ****
# Top 50 median cov
df_top50_median <- df_median_myop_panel %>%
  slice_tail(n=40)
top_50_median_cov_gene <- df_top50_median[["Gene"]]

# Plot bottom 50 cov median
filtered_cov_sum_median_top50 <- df_sorted_cov_sum %>% 
  filter(Gene %in% top_50_median_cov_gene)
filtered_cov_sum_median_top50$Gene <- factor(filtered_cov_sum_median_top50$Gene, levels = top_50_median_cov_gene)
top50_plot <- plot_gene_cov_distribution(filtered_cov_sum_median_top50, "Coverage distribution per gene top 50 median")


# # Ordering bottom 50 medians genes with bottom_50_median_cov_gene
# filtered_cov_sum_median_bottom50 <- filtered_cov_sum_median_bottom50 %>%
#   mutate(color_group = ifelse(Gene %in% bottom50_same_mean_median, "Both_mean_median", "Different"))
# filtered_cov_sum_median_bottom50$Gene <- factor(filtered_cov_sum_median_bottom50$Gene, levels = bottom_50_median_cov_gene)
# 
# 
# filtered_cov_sum_median_top50 <- filtered_cov_sum_median_top50 %>%
#   mutate(color_group = ifelse(Gene %in% top50_same_mean_median, "Both_mean_median", "Different"))
# filtered_cov_sum_median_top50$Gene <- factor(filtered_cov_sum_median_top50$Gene, levels = top_50_median_cov_gene)
# plot_cov_gene_highlight(filtered_cov_sum_median_top50, "Coverage distribution per gene top 50 median highlight diff to mean", custom_colors)

combined_plot <- bottom50_plot + top50_plot
combined_plot & scale_y_continuous(limits = c(0, 300))


