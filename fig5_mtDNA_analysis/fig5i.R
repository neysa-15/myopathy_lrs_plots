setwd("/path/to/analysis_dir")

.libPaths(c("/path/to/R_libs"))

library(data.table)
library(ggplot2)
library(patchwork)

sample_key <- "/path/to/sample_key.tsv"
sample_key <- fread(sample_key, sep = "\t", header = FALSE)
names(sample_key) <- c("LRS_ID","Sample")
sample_key[, LRS_ID_short := gsub("RS0*", "", LRS_ID)]
sample_key <- as.data.frame(sample_key)

# Add SRA sample from https://pmc.ncbi.nlm.nih.gov/articles/PMC10265159/ 
ont_muscle_row <- data.frame(
  LRS_ID = c("SRR23886856", "SRR23886840", "SRR23886855", "SRR23886851", "SRR23886844", "SRR23886845", "SRR23886841", "SRR23886842", "SRR23886843", "SRR23886857", "SRR23886858", "SRR23886859", "SRR23886860", "SRR23886861", "SRR23886862"),
  Sample = c("SRR23886856", "SRR23886840", "SRR23886855", "SRR23886851", "SRR23886844", "SRR23886845", "SRR23886841", "SRR23886842", "SRR23886843", "SRR23886857", "SRR23886858", "SRR23886859", "SRR23886860", "SRR23886861", "SRR23886862"),
  LRS_ID_short = c("SRR23886856", "SRR23886840", "SRR23886855", "SRR23886851", "SRR23886844", "SRR23886845", "SRR23886841", "SRR23886842", "SRR23886843", "SRR23886857", "SRR23886858", "SRR23886859", "SRR23886860", "SRR23886861", "SRR23886862"),
  stringsAsFactors = FALSE
)
sample_key <- rbind(sample_key, ont_muscle_row)

# read sample proportion file
sample_bigdel_prop_summary <- "/path/to/sample_proportion_summary_1kbp_deletion.tsv"

sample_bigdel_prop_summary <- fread(sample_bigdel_prop_summary, sep="\t")
sample_bigdel_prop_summary <- as.data.frame(sample_bigdel_prop_summary)
# Rename header
sample_bigdel_prop_summary <- sample_bigdel_prop_summary %>% rename(LRS_ID = Sample)

bigdel_summary_df <- inner_join(sample_key, sample_bigdel_prop_summary, by="LRS_ID")

# making proportion numeric
bigdel_summary_df <- bigdel_summary_df %>%
  mutate(
    # Remove '%' if it exists and convert to numeric
    Proportion_numeric = as.numeric(gsub("%", "", Proportion_with_big_deletion)),
    Proportion_numeric_over1kgp = as.numeric(gsub("%", "", Proportion_with_big_deletion_over1kbp))
  )

bigdel_summary_df <- bigdel_summary_df %>%
  mutate(Category = case_when(
    grepl("Sample_of_interest", Sample) ~ "patient_muscle",
    grepl("Sample_of_interest", Sample) ~ "patient_blood",
    grepl("SRR", Sample) ~ "sra_muscle",
    TRUE ~ "random_myop_blood"
  )) %>%
  arrange(Sample)

# map colors
cat_colors <- c("random_myop_blood" = "black", "patient_blood" = "red", "patient_muscle" = "blue", "sra_muscle"="darkgrey")

# Order by category
category_order <- c("random_myop_blood", "patient_blood", "patient_muscle", "sra_muscle")
bigdel_summary_df$Category <- factor(bigdel_summary_df$Category, levels = category_order)

bigdel_summary_df <- bigdel_summary_df %>%
  arrange(Category) %>%
  mutate(LRS_ID_short = factor(LRS_ID_short, levels = unique(LRS_ID_short)))

label_colors <- cat_colors[bigdel_summary_df$Category]

# Plot over 1 kgp
proportion_allover1kbp_plot <- ggplot(bigdel_summary_df, aes(x = LRS_ID_short, y = Proportion_numeric_over1kgp, fill = Category)) +
  geom_col() +
  # This links the bars to the colors and creates the legend
  scale_fill_manual(values = cat_colors) + 
  theme_bw(base_size = 14) +
  theme(
    # This colors the actual text on the X axis
    axis.text.x = element_text(angle = 45, hjust = 1, color = label_colors),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "Proportion of chrM reads with >1000bp deletion",
    x = "Sample",
    y = "Proportion (%)",
    fill = "Sample Group" # Changes the title of the legend
  )

print(proportion_allover1kbp_plot)

##
# FILTERED DF - only muscle samples
muscle_df <- bigdel_summary_df[bigdel_summary_df$Category %in% c("sra_muscle", "patient_muscle"), ]

# Add age if want to order by age
age_df <- fread("/path/to/sample_age.tsv", sep='\t')
age_df <- age_df[, c("Sample", "Age")]

merged_df <- left_join(muscle_df, age_df, by = c("LRS_ID" = "Sample"))

plot_data <- merged_df %>%
  arrange(Age) %>%
  mutate(LRS_ID_short = factor(LRS_ID_short, levels = unique(LRS_ID_short)))

ggplot(plot_data, aes(x = LRS_ID_short, y = Proportion_numeric_over1kgp, fill = Category)) +
  geom_col() +
  scale_fill_manual(values = cat_colors) + 
  theme_bw(base_size = 14) +
  theme(
    # color text on the X axis
    axis.text.x = element_text(angle = 45, hjust = 1, color = label_colors),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "Proportion of chrM reads with >1000bp deletion",
    x = "Sample",
    y = "Proportion (%)",
    fill = "Sample Group" # legend title
  )
