setwd("/g/data/kr68/neysa/r_plotting")

.libPaths(c("/g/data/kr68/andre/R_libs"))

library(data.table)
library(ggplot2)
library(tidyr)
library(patchwork)

myop_matrix <- fread("/g/data/kr68/neysa/r_plotting/demographics_matrix.tsv", header = TRUE)
setnames(myop_matrix, "ID", "Sample")

sample_key <- fread("/g/data/kr68/puzzleapp/KISKUM_Myop/KISKUM_Myop.sample_key.tsv", sep = "\t", header = FALSE)
names(sample_key) <- c("LRS_ID", "Sample")
sample_key[, LRS_ID := gsub("RS0*", "", LRS_ID)]

# Merge annotation
myop_matrix <- merge(myop_matrix, sample_key, by = "Sample", all.x = TRUE)

previously_unsolved <- c("L236", "L286", "L143", "L126", "L136", "L248", "L437", "L157", "L147", "L121", "L140", 
                         "L234", "L235", "L122", "L123", "L125", "L134", "L158", "L119", "L159", "L130", "L135", 
                         "L131", "L132", "L133", "L129", "L155", "L156", "L146", "L249", "L323")

# Convert to df
myop_df <- as.data.frame(myop_matrix)

# Add missing LRS_ID from sample_key
myop_df$LRS_ID <- as.character(myop_df$LRS_ID)
myop_df$LRS_ID[myop_df$Sample == "R220038"] <- "L286"

# Order
myop_df$LRS_ID <- factor(myop_df$LRS_ID, levels = previously_unsolved)
myop_df_ordered <- myop_df[order(myop_df$LRS_ID), ]

myop_df_ordered <- myop_df_ordered %>% 
  rename("Sex assigned at birth" = "Gender")

# Add matrix row for final diagnosis
myop_df_ordered$Diagnosis <- c(
  "LGMD", "OPDM", "OPDM", "OPDM", "OPDM", "OPDM", "OPDM", "OPDM", "FSHD1",
  "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""
)

# ------------------------------------------------------------------------------
# Plot age
# ------------------------------------------------------------------------------

# change format to long df
df_long <- myop_df_ordered %>%
  
  # Pivot the two age columns into one "Metric" column and one "Years" column
  pivot_longer(
    cols = c("Disease duration at time of recruitment (years)", "Age at onset (years)"), # The columns to stack
    names_to = "Metric",                         # New column for the names
    values_to = "Years"                          # New column for the values
  )

# View the new "long" data
print(df_long)

plot_age <- ggplot(df_long, aes(x = LRS_ID, y = Years, fill = Metric)) +
  geom_col(position = position_stack(reverse = TRUE)) +  # geom_col() is the same as geom_bar(stat="identity")
  
  scale_fill_manual(
    values = c(
      "Age at onset (years)" = "#3e7d81",
      "Disease duration at time of recruitment (years)" = "#DAB1DA"
    ),
    na.value = "white"   # Blanks (NA) are white
  ) +
  
  # Add nice labels
  labs(
    title = "Age at Onset and Disease Duration by Sample",
    x = "Sample ID",
    y = "Total Years",
    fill = "Age Component" # This renames the legend title
  ) +
  
  # A clean theme
  theme_minimal() +
  
  # Rotate x-axis labels if they overlap
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ------------------------------------------------------------------------------
# Plot matrix / heatmap
# ------------------------------------------------------------------------------

# Define metrics order
metric_order <- c(
  "Diagnosis", "Sex assigned at birth", "Targeted genel panel", "Exome-based neuromuscular panel",
  "Genome sequencing", "FSHD SB", "DM1 test", "DM2 test", "DMD MLPA", "Muscle biopsy"
)

# Pivot the data to "long" format for plotting
df_long_matrix <- myop_df_ordered %>%
  
  # Pivot all columns *except* LRS_ID
  pivot_longer(
    cols = -c(Sample, LRS_ID, `Age at onset (years)`, `Aget at time of recruitment (years)`, `Disease duration at time of recruitment (years)`),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  
  # 4. Apply the custom ordering
  mutate(
    # Set LRS_ID order (x-axis)
    LRS_ID = factor(LRS_ID, levels = previously_unsolved),
    # Set Metric order (y-axis, reversed for plotting)
    Metric = factor(Metric, levels = rev(metric_order))
  )

# 5. Create the plot
plot_matrix <- ggplot(df_long_matrix, aes(x = LRS_ID, y = Metric, fill = Value)) +
  
  # Add the tiles with a white border
  geom_tile(color = "black", linewidth = 0.2) +
  
  # Apply your custom color rules
  scale_fill_manual(
    values = c(
      "LGMD" = "#D83034",
      "OPDM" = "#213f7a",
      "FSHD1" = "#f173ac",
      "F" = "pink",      # F is pink
      "M" = "royalblue", # M is blue
      "Y" = "lightgreen"     # Y is yellow
    ),
    na.value = "white"   # Blanks (NA) are white
  ) +
  
  # Clean up the theme
  theme_minimal() +
  
  # Add labels
  labs(
    x = "LRS_ID",
    y = NULL,
    fill = "Result"
  ) +
  
  # Rotate x-axis labels
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  
  coord_fixed()


# ------------------------------------------------------------------------------
# Combine plot
# ------------------------------------------------------------------------------
plot_age_noxaxis <- plot_age + 
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) #+
  #theme(legend.position = "none")

plot_age_noxaxis / plot_matrix +
  plot_layout(guides = 'collect') +
  plot_layout(heights = c(1, 1))

