setwd("/g/data/kr68/neysa/r_plotting/rcode_per_fig")

# ------------------------------------------------------------------------------
#                 Libraries
# ------------------------------------------------------------------------------
.libPaths(c("/g/data/kr68/andre/R_libs"))

library(data.table)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
#                 Data loading and processing
# ------------------------------------------------------------------------------
cov_sum <- fread("/g/data/kr68/neysa/fshd_pipeline/coverage_analysis/summaries/sample_coverage_summary_complete_noreadlengths.tsv", header = TRUE)
setnames(cov_sum, "sample", "Sample")

# Load and clean sample key
sample_key <- fread("/g/data/kr68/puzzleapp/KISKUM_Myop/KISKUM_Myop.sample_key.tsv", sep = "\t", header = FALSE)
names(sample_key) <- c("LRS_ID", "Sample")
sample_key[, LRS_ID := gsub("RS0*", "", LRS_ID)]

myopathy <- c(
  "ZE2607", "RC1309", "BH0608", "ZD0608", "IB2806", "SA1110", "VQ2510", "BC2211", 
  "JOUB61166", "AS2603", "R230025", "JURA89", "KAHO2804", "JOBO3009", "RJ1207", 
  "GL2106", "R240177", "R240183", "R240059", "DL1104", "PN1206", "R220038", 
  "BP0703", "JZ2510", "AK2208", "R240186", "R240088", "CF2608", "GUAT0705", 
  "SAHI0207", "DOHO2501", "EW5762", "QOL0607", "ZB1207", "BA0908", "ZU1108", 
  "BF1708", "PS1509", "LU1110", "ZL0811", "FK1411", "ZL2011", "KN2211", "DC2702", 
  "BQ1303", "QQ0805", "BV2705", "WQ2407", "R250002", "R250028"
)

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
#      Sup Fig 2b - Coverage per target of each sample
# ------------------------------------------------------------------------------
# enforce order on Sample
cov_sum$Sample <- factor(cov_sum$Sample, levels = myopathy)

# get unique mapping of Sample → LRS_ID in myopathy order
id_order <- cov_sum[!duplicated(cov_sum$Sample), c("Sample", "LRS_ID")]
id_order <- id_order[order(factor(id_order$Sample, levels = myopathy)), ]

# set factor levels of LRS_ID based on that order
cov_sum$LRS_ID <- factor(cov_sum$LRS_ID, levels = id_order$LRS_ID)

cov_sum_filter <- cov_sum[!is.na(cov_sum$Sample), ]

# plot
ggplot(cov_sum_filter, aes(x = LRS_ID, y = Coverage)) +
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
#      Sup Fig 2c - ON-TARGET and OFF-TARGET PLOTS FOR N50
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