setwd("/g/data/kr68/neysa/r_plotting/rcode_per_fig")

# ------------------------------------------------------------------------------
#                 Libraries
# ------------------------------------------------------------------------------
.libPaths(c("/g/data/kr68/andre/R_libs"))

library(data.table)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
#                 Variables
# ------------------------------------------------------------------------------
# Directory of d4z4ling results
base_dir <- "/g/data/kr68/fshd/results_mapq0"

# Grouping of FSHD group
fshd1 <- c("JOUB61166","AS2603","R230025","JURA89","KAHO2804","JOBO3009","RJ1207","GL2106","R240177","R240183","CF2608")
fshd2 <- c("R240059","DL1104","R250109")
fshd1_n_2 <- c("R250119")
fshd1_borderline <- c("GUAT0705")
non_fshd <- c(
  "ZE2607", "RC1309", "BH0608", "ZD0608", "IB2806", "SA1110", "VQ2510", "BC2211", "PN1206", "R220038", "BP0703", "JZ2510", "AK2208", "R240186", "R250113",
  "R240088", "SAHI0207", "DOHO2501", "EW5762", "QOL0607", "ZB1207", "BA0908", "ZU1108", "BF1708", "PS1509", "LU1110", "ZL0811", "FK1411",
  "ZL2011", "KN2211", "DC2702", "BQ1303", "QQ0805", "BV2705", "WQ2407", "R250002", "R250028"
)

# Sample name conversion of different identifier
sample_key <- "/g/data/kr68/puzzleapp/KISKUM_Myop/KISKUM_Myop.sample_key.tsv"
sample_key <- fread(sample_key, sep = "\t", header = FALSE)
names(sample_key) <- c("LRS_ID","Sample")
sample_key[, LRS_ID := gsub("RS0*", "", LRS_ID)]
sample_key <- as.data.frame(sample_key)

# ------------------------------------------------------------------------------
#                 Colors
# ------------------------------------------------------------------------------
label_colours <- c(
  "FSHD1" = "#ffe0e0",
  "FSHD2" = "#e0f7fa",
  "FSHD1+2"= "#efe0ff",
  "Borderline FSHD1" = "#D6FAC3",
  "FSHD Negative" = "#e0e0e0"
)

# ------------------------------------------------------------------------------
#                 Data input and processing
# ------------------------------------------------------------------------------
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
# combined_dt[Sample %in% fshd1_undiagnosed, Sample_Label := "Undiagnosed FSHD1"]
combined_dt[Sample %in% fshd2, Sample_Label := "FSHD2"]
combined_dt[Sample %in% fshd1_n_2, Sample_Label := "FSHD1+2"] 
combined_dt[Sample %in% fshd1_borderline, Sample_Label := "Borderline FSHD1"]
combined_dt[Sample %in% non_fshd, Sample_Label := "FSHD Negative"]

# Define desired order of Sample_Label
label_order <- c("FSHD1", "FSHD1+2", "FSHD2", "Borderline FSHD1", "FSHD Negative")
combined_dt[, Sample_Label := factor(Sample_Label, levels = label_order)]

# Order Sample within each Sample_Label
combined_dt[, Sample := factor(Sample, levels = combined_dt[order(Sample_Label, Sample), unique(Sample)])]

# Filter out samples that's not in any group
combined_dt <- combined_dt[!is.na(Sample_Label)]

# Sample ordering
# filtered_dt <- combined_dt[Sample_Label %in% c("FSHD1", "FSHD2", "FSHD1+2", "Borderline FSHD1")]
# plot_dt <- filtered_dt[
#   grepl("4qA|4qB", Haplotype)
# ]
# 
# all_samples_ordered <- c(fshd1, fshd1_n_2, fshd2, fshd1_borderline)
# 
# # 2. Convert the 'Sample' column to a factor using this new master list
# plot_dt[, Sample := factor(Sample, levels = all_samples_ordered)]

# ------------------------------------------------------------------------------
#                 Figure 2b and Sup Figure 5c
# ------------------------------------------------------------------------------
plot_fshd_5mC_vs_d4z4copies_4qA <- function(sample_id, plot4q) {
  
  control_sample_dt <- plot_4q[MappedEstimatedCopies>=2 & Sample == sample_id]
  max_number_copies <- max(as.integer(control_sample_dt$MappedEstimatedCopies), na.rm = TRUE)
  
  ggplot(control_sample_dt, aes(x = pLAM_Methylation_Percentage, y = floor(MappedEstimatedCopies), color = ReadLabel)) +
    geom_point(size = 0.5, alpha = 0.6) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "black", size = 0.2) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "black", size = 0.2) +
    scale_x_continuous(limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, max_number_copies + 5)) +
    scale_color_manual(values = c("Complete 4qA" = "#1f77b4", "Partial distal 4qA" = "#ff7f0e")) +
    labs(x = "Distal D4Z4 5mC (%)", y = "Number of D4Z4 copies", color = "Read type") +
    theme_bw(base_size = 5) +
    theme(legend.position = "bottom",
          aspect.ratio = 1,
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
}

plot_4q <- combined_dt[
  pLAM_mapped == TRUE &
    grepl("^chr4:", GenomeCoords) &
    grepl("Complete 4qA|Partial distal 4qA", ReadLabel)
]

# Fig2 - FSHD1 Positive
sample_id <- "JOUB61166"
plot_fshd_5mC_vs_d4z4copies_4qA(sample_id, plot_4q)

# Resizing plot for manuscript
control_panel <- set_panel_size(
  p = plot_fshd_5mC_vs_d4z4copies_4qA(sample_id, plot_4q),
  width = unit(4, "cm"), # Example width
  height = unit(4, "cm")  # Example height
)
filename <- paste(sample_id, "_d4z4vsmethy.pdf", sep = "")
ggsave(filename, plot = control_panel)

# Fig2 - non-FSHD
sample_id <- "R250002"
plot_fshd_5mC_vs_d4z4copies_4qA(sample_id, plot_4q)

# Resizing plot for manuscript
control_panel <- set_panel_size(
  p = plot_fshd_5mC_vs_d4z4copies_4qA(sample_id, plot_4q),
  width = unit(4, "cm"), # Example width
  height = unit(4, "cm")  # Example height
)
filename <- paste(sample_id, "_d4z4vsmethy.pdf", sep = "")
ggsave(filename, plot = control_panel)

# # Fig2 - FSHD2
sample_id <- "DL1104"
####################
# Facet plot for FSHD2 and FSHD1+2
plot_fshd_5mC_vs_d4z4copies_facet <- function(sample_id) {
  ggplot(combined_dt[Sample==sample_id & !grepl("Unclassified",ReadLabel) & MappedEstimatedCopies>=2] , aes(x = pLAM_Methylation_Percentage, y = floor(MappedEstimatedCopies), color = ReadLabel)) +
    geom_point(size = 1.2, alpha = 0.6) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 50, linetype = "dashed", color = "black") +
    scale_x_continuous(limits = c(0, 100)) +
    # scale_y_continuous(limits = c(0, max_number_copies + 5)) +
    scale_color_manual(values = c("Complete 4qA" = "#1f77b4", "Partial distal 4qA" = "#ff7f0e", "Complete 10qA" = "#f3766d", "Partial distal 10qA" = "#25b151", "Partial distal 4qB" = "purple")) +
    labs(x = "Distal D4Z4 5mC (%)", y = "Number of D4Z4 copies", color = "Read type", title = sample_id ) +
    facet_wrap(~ Haplotype, drop = FALSE) +
    theme_bw(base_size = 5) +
    theme(legend.position = "bottom",
          aspect.ratio = 1,
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

# # Fig2 - FSHD2
sample_id <- "DL1104"
fshd2_facet_plot <- plot_fshd_5mC_vs_d4z4copies_facet(sample_id)

fshd2_panel <- set_panel_size(
  p = fshd2_facet_plot,
  width = unit(4, "cm"), # Example width
  height = unit(4, "cm")  # Example height
)
filename <- paste(sample_id, "_d4z4vsmethy.pdf", sep = "")
ggsave(filename, plot = fshd2_panel)

# Supplementary Fig5c - FSHD1+2
sample_id <- "R250119"
fshd1_2 <- plot_fshd_5mC_vs_d4z4copies_facet(sample_id)

fshd1_2_panel <- set_panel_size(
  p = fshd1_2,
  width = unit(4, "cm"), # Example width
  height = unit(4, "cm")  # Example height
)
filename <- paste(sample_id, "_d4z4vsmethy.pdf", sep = "")
ggsave(filename, plot = fshd1_2_panel)

# ------------------------------------------------------------------------------
#                 Figure 2c
# ------------------------------------------------------------------------------
filtered_dt <- combined_dt[Sample_Label %in% c("FSHD1", "FSHD2", "FSHD1+2", "Borderline FSHD1")]

# Take all 4q haplotype
plot_dt <- filtered_dt[
  grepl("4qA|4qB", Haplotype)
]

plot_dt[, Sample_Label := factor(Sample_Label, levels = label_order)]
sample_order_dt <- unique(plot_dt[, .(Sample, Sample_Label)])
setorder(sample_order_dt, Sample_Label, Sample)
ordered_samples <- sample_order_dt$Sample
plot_dt[, Sample := factor(Sample, levels = ordered_samples)]

setorder(plot_dt, Sample)

positive_samples <- unique(plot_dt$Sample)

filtered_key_dt <- sample_key[plot_dt$Sample %in% positive_samples, ]
sample_labels <- setNames(filtered_key_dt$LRS_ID, filtered_key_dt$Sample)

# Filter plot data:
filtered_plot_dt <- plot_dt[
  !(grepl("Partial distal", ReadLabel) & MappedEstimatedCopies < 3)
]

# change bg_df range
# Background shading
bg_df2 <- unique(filtered_plot_dt[, .(Sample, Sample_Label)])
bg_df2[, Sample := factor(Sample, levels = levels(filtered_plot_dt$Sample))]
bg_df2[, Sample_Num := as.numeric(Sample)]
bg_df2[, `:=`(
  xmin = Sample_Num - 0.5,
  xmax = Sample_Num + 0.5,
  fill = label_colours[Sample_Label]
)]


# 1f77b4  #ff7f0e  #9467bd
read_label_colors <- c(
  "Complete 4qA" = "#00429d", 
  "Partial distal 4qA" = "#ff7f0e",
  "Complete 4qB" = "#6a3d9a",
  "Partial distal 4qB" = "#255f38",
  "Complete 10qA" = "#e377c2",
  "Partial distal 10qA" = "#d62728",
  "Complete 10qB" = "#7f7f7f",
  "Partial distal 10qB" = "brown"
)

filtered_plot_dt <- filtered_plot_dt %>%
  mutate(Haplotype = case_when(
    grepl("4qA", ReadLabel) ~ "4qA",
    grepl("4qB", ReadLabel) ~ "4qB",
    grepl("10q", ReadLabel) ~ "10q", # Add other cases as needed
    TRUE ~ "Other" # A fallback for any other labels
  ))

plot_data_sorted <- filtered_plot_dt %>%
  arrange(grepl("Complete", ReadLabel))
# 

############
# d4z4 copies for postivive and borderline fshd cases
d4z4_non_neg_case <- ggplot(plot_data_sorted, aes(x = Sample, y = floor(MappedEstimatedCopies), color = ReadLabel, shape = Haplotype)) + # <-- ADD shape = Allele HERE
  geom_rect(
    data = bg_df2,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Sample_Label),
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(label_colours),
    breaks = c(names(label_colours)),
    name = NULL,
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "black") +
  scale_color_manual(
    values = read_label_colors,
    name = NULL
  ) +
  
  geom_jitter(alpha = 0.4, width = 0.25, height = 0.1, size = 1.5, stroke = 0.2) + # Added stroke for outlines
  
  scale_shape_manual(
    name = "Haplotype", # Add a title for the shape legend
    values = c("4qA" = 16, "4qB" = 17, "10q" = 15)
  ) +
  
  scale_x_discrete(labels = sample_labels) +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = c(0, 5, 10, 15, 20)
  ) +
  labs(x = "Sample", y = "Number of D4Z4 copies") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

# 2. Use set_panel_size() to create a new object with a fixed panel
#    You can set any dimensions you need here.
plot_with_fixed_panel <- set_panel_size(
  p = d4z4_non_neg_case,
  width = unit(10, "cm"), # Example width
  height = unit(5, "cm")  # Example height
)

# 3. Save the new, modified plot object
ggsave("d4z4_copies_plot_fixed.pdf", plot = plot_with_fixed_panel)

# ------------------------------------------------------------------------------
#                 Figure 2d - Limit of detection analysis
# ------------------------------------------------------------------------------

# Configuration
base_output_dir <- "/g/data/kr68/andre/fshd_pipeline/subsampling_experiment"
sample_ids <- c("AS2603", "JOUB61166", "RJ1207", "GUAT0705")
sd_settings <- c("5sd", "10sd", "10sd_50draws")
combined_output_file <- file.path(base_output_dir, "combined_subsampling_summary.tsv")

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

probability_dt_wide <- combined_summary_dt[SD_Setting == plot_setting, {
  .(
    Prob_at_least_1 = mean(Complete_4qA_Reads >= 1),
    Prob_at_least_2 = mean(Complete_4qA_Reads >= 2),
    Prob_at_least_3 = mean(Complete_4qA_Reads >= 3)
  )
}, by = .(SampleID, Multiplex_Factor)]

probability_dt_long <- melt(probability_dt_wide, 
                            id.vars = c("SampleID", "Multiplex_Factor"),
                            measure.vars = c("Prob_at_least_1", "Prob_at_least_2", "Prob_at_least_3"),
                            variable.name = "Threshold",
                            value.name = "Probability")

combined_summary_dt[, Multiplex_Factor := as.factor(Multiplex_Factor)] 

probability_dt_long[, Multiplex_Factor := as.factor(Multiplex_Factor)] 

probability_dt_long[, Threshold := fcase(
  Threshold == "Prob_at_least_1", ">= 1 Read",
  Threshold == "Prob_at_least_2", ">= 2 Reads",
  Threshold == "Prob_at_least_3", ">= 3 Reads"
)]
lod <- ggplot(probability_dt_long, aes(x = Multiplex_Factor, y = Probability, color = Threshold, group = Threshold)) +
  geom_line(linewidth = 0.15) + 
  geom_point(size = 1.0) + 
  facet_wrap(~ SampleID, nrow = 2) + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Probability of Detecting Complete 4qA Reads by Threshold",
    x = "Multiplex Factor",
    y = "Probability",
    color = "Detection Threshold"
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

# Determine the number of panels automatically
num_panels <- length(unique(probability_dt_long$SampleID))

# Calculate the required width for a single panel
single_panel_width_cm <- 8 / num_panels # 10 cm total width / 4 panels = 2.5 cm

# Apply set_panel_size() with the calculated dimensions
lod_fixed_panels <- set_panel_size(
  p = lod,
  width = unit(single_panel_width_cm, "cm"), # Set each panel to 2.5 cm wide
  height = unit(5, "cm")                     # Set each panel to 5 cm high
)

# 2x2
lod_fixed_panels <- set_panel_size(
  p = lod,
  width = unit(4, "cm"), # Set each panel to 2.5 cm wide
  height = unit(4, "cm")                     # Set each panel to 5 cm high
)

# Save the modified plot
# The final image will have a total panel width of 10 cm (4 * 2.5 cm)
ggsave("lod_plot.pdf", plot = lod_fixed_panels)

# ------------------------------------------------------------------------------
#                 Figure 2e
# ------------------------------------------------------------------------------

plot_methylation <- function(dt, region = "chr4", read_pattern = "4qA", 
                              read_labels = c("Complete", "Partial distal"), 
                              neg_control_label = "FSHD Negative") {
  
  label_colours <- c(
    "FSHD1" = "#ffe0e0",
    "FSHD1+2"= "#efe0ff",
    "FSHD2" = "#e0f7fa",
    "Borderline FSHD1" = "#D6FAC3",
    "FSHD Negative" = "#e0e0e0"
  )
  
  label_order <- c("FSHD1", "FSHD1+2", "FSHD2", "Borderline FSHD1", "FSHD Negative")
  
  # Filter data
  pattern <- paste(read_labels, read_pattern, sep = " ", collapse = "|")
  plot_dt <- dt[
    pLAM_mapped == TRUE &
      grepl(paste0("^", region, ":"), GenomeCoords) &
      grepl(pattern, ReadLabel)
  ]
  
  # Define ReadCategory
  plot_dt[, ReadCategory := fifelse(
    grepl(paste0("Complete ", read_pattern), ReadLabel),
    paste("Complete", read_pattern),
    paste("Partial distal", read_pattern)
  )]
  
  # Reorder sample factors
  plot_dt[, Sample_Label := factor(Sample_Label, levels = label_order)]
  sample_order_dt <- unique(plot_dt[, .(Sample, Sample_Label)])
  setorder(sample_order_dt, Sample_Label, Sample)
  plot_dt[, Sample := factor(Sample, levels = sample_order_dt$Sample)]
  
  # Background shading
  bg_df <- unique(plot_dt[, .(Sample, Sample_Label)])
  bg_df[, Sample := factor(Sample, levels = levels(plot_dt$Sample))]
  bg_df[, Sample_Num := as.numeric(Sample)]
  bg_df[, `:=`(
    xmin = Sample_Num - 0.5,
    xmax = Sample_Num + 0.5,
    fill = label_colours[Sample_Label]
  )]
  
  # Reference lines
  neg_ctrl_values <- plot_dt[
    Sample_Label == neg_control_label & 
      ReadCategory == paste("Partial distal", read_pattern),
    pLAM_Methylation_Percentage
  ]
  neg_ctrl_line <- list(
    x_start = min(bg_df[Sample_Label == neg_control_label]$Sample_Num) - 0.5,
    x_end   = max(bg_df[Sample_Label == neg_control_label]$Sample_Num) + 0.5,
    y       = median(neg_ctrl_values, na.rm = TRUE),
    mad     = mad(neg_ctrl_values, constant = 1, na.rm = TRUE)
  )
  
  fshd1_values <- plot_dt[
    Sample_Label == "FSHD1" & 
      ReadCategory == paste("Complete", read_pattern) &
      Sample != "GL2106",
    pLAM_Methylation_Percentage
  ]
  fshd1_line <- list(
    x_start = min(bg_df[Sample_Label == "FSHD1"]$Sample_Num) - 0.5,
    x_end   = max(bg_df[Sample_Label == "FSHD1"]$Sample_Num) + 0.5,
    y       = median(fshd1_values, na.rm = TRUE),
    mad     = mad(fshd1_values, constant = 1, na.rm = TRUE)
  )
  
  sample_labels <- setNames(sample_key$LRS_ID, sample_key$Sample)
  
  # Colour mapping
  custom_colours <- c(
    label_colours,
    setNames(c("#1f77b4", "#ff7f0e"),
             c(paste("Complete", read_pattern), paste("Partial distal", read_pattern)))
  )
  
  # Plot
  ggplot(plot_dt[MappedEstimatedCopies >= 1], aes(x = Sample, y = pLAM_Methylation_Percentage, fill = ReadCategory)) +
    geom_rect(
      data = bg_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Sample_Label),
      inherit.aes = FALSE, alpha = 0.8
    ) +
    geom_boxplot(position = position_dodge(preserve = "single"), alpha = 0.4, outlier.shape = NA, width = 0.5, linewidth = 0.2) +
    geom_segment(
      aes(x = neg_ctrl_line$x_start, xend = neg_ctrl_line$x_end,
          y = neg_ctrl_line$y, yend = neg_ctrl_line$y),
      inherit.aes = FALSE, linetype = "solid", colour = "#e65c00", size = 0.3
    ) +
    geom_segment(
      aes(x = fshd1_line$x_start, xend = fshd1_line$x_end,
          y = fshd1_line$y, yend = fshd1_line$y),
      inherit.aes = FALSE, linetype = "solid", colour = "#115e99", size = 0.3
    ) +
    geom_segment(
      aes(x = neg_ctrl_line$x_start, xend = neg_ctrl_line$x_end,
          y = neg_ctrl_line$y + 2 * neg_ctrl_line$mad,
          yend = neg_ctrl_line$y + 2 * neg_ctrl_line$mad),
      inherit.aes = FALSE, linetype = "dashed", colour = "#e65c00", size = 0.3
    ) +
    geom_segment(
      aes(x = neg_ctrl_line$x_start, xend = neg_ctrl_line$x_end,
          y = neg_ctrl_line$y - 2 * neg_ctrl_line$mad,
          yend = neg_ctrl_line$y - 2 * neg_ctrl_line$mad),
      inherit.aes = FALSE, linetype = "dashed", colour = "#e65c00", size = 0.3
    ) +
    geom_segment(
      aes(x = fshd1_line$x_start, xend = fshd1_line$x_end,
          y = fshd1_line$y + 2 * fshd1_line$mad,
          yend = fshd1_line$y + 2 * fshd1_line$mad),
      inherit.aes = FALSE, linetype = "dashed", colour = "#115e99", size = 0.3
    ) +
    geom_segment(
      aes(x = fshd1_line$x_start, xend = fshd1_line$x_end,
          y = fshd1_line$y - 2 * fshd1_line$mad,
          yend = fshd1_line$y - 2 * fshd1_line$mad),
      inherit.aes = FALSE, linetype = "dashed", colour = "#115e99", size = 0.3
    ) +
    scale_x_discrete(labels = sample_labels) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_fill_manual(
      values = custom_colours,
      breaks = c(paste("Complete", read_pattern), paste("Partial distal", read_pattern), names(label_colours)),
      name = NULL,
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    labs(x = "Sample", y = "Distal D4Z4 5mC (%)", fill = "Read Category") +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )
}

###############
# Only positive and borderline

# Filter the directory
# Combine all IDs
all_ids <- c(fshd1, fshd1_n_2, fshd2, fshd1_borderline)

# Get only subdirectories whose *basename* matches one of the IDs
all_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
selected_dirs <- all_dirs[basename(all_dirs) %in% all_ids]

selected_summary_files <- unlist(lapply(selected_dirs, function(dir) {
  list.files(path = dir, pattern = "mapped_features_summary\\.tsv$", full.names = TRUE)
}))

non_neg_dt <- rbindlist(lapply(selected_summary_files, function(f) {
  dt <- fread(f)
  sample_name <- sub("_mapped_features_summary\\.tsv$", "", basename(f))
  dt[, Sample := sample_name]
  return(dt)
}), fill = TRUE)

non_neg_dt[Sample %in% fshd1, Sample_Label := "FSHD1"]
non_neg_dt[Sample %in% fshd1_n_2, Sample_Label := "FSHD1+2"]
non_neg_dt[Sample %in% fshd1_borderline, Sample_Label := "Borderline FSHD1"]
non_neg_dt[Sample %in% fshd2, Sample_Label := "FSHD2"]

# Define desired order of Sample_Label
label_order <- c("FSHD1", "FSHD1+2", "FSHD2", "Borderline FSHD1")
non_neg_dt[, Sample_Label := factor(Sample_Label, levels = label_order)]

fshd_methylation_plot <- plot_methylation(non_neg_dt, region = "chr4", read_pattern = "4qA")
# ggsave("fshd_methylation_plot.pdf", plot = fshd_methylation_plot, width = 36, height = 24, units = "cm")

###############
# Negative controls

# other_dirs <- all_dirs[!basename(all_dirs) %in% all_ids]
other_dirs <- all_dirs[basename(all_dirs) %in% non_fshd]

neg_controls_summary_files <- unlist(lapply(other_dirs, function(dir) {
  list.files(path = dir, pattern = "mapped_features_summary\\.tsv$", full.names = TRUE)
}))

neg_dt <- rbindlist(lapply(neg_controls_summary_files, function(f) {
  dt <- fread(f)
  sample_name <- sub("_mapped_features_summary\\.tsv$", "", basename(f))
  dt[, Sample := sample_name]
  return(dt)
}), fill = TRUE)

neg_dt[Sample %in% non_fshd, Sample_Label := "FSHD Negative"]

fshd_methylation_plot_neg <- plot_methylation(neg_dt, region = "chr4", read_pattern = "4qA")
ggsave("fshd_methylation_plot_negative.pdf", plot = fshd_methylation_plot_neg, width = 36, height = 24, units = "cm")

# Merge and save
combined_plot <- fshd_methylation_plot + fshd_methylation_plot_neg + plot_layout(ncol = 2, widths = c(1, 2))
ggsave("combined_methylation_plot.pdf", plot = combined_plot, width = 72, height = 24, units = "cm")

