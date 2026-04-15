.libPaths(c("/path/to/R_libs", .libPaths()))
library(data.table)
library(ggplot2)
library(ggrepel)

# read
dt <- fread("/path/to/all_STR_loci.all_cohorts.tsv", na.strings=".")
setnames(dt, c("GENE","SAMPLE_ID","ALLELE_ID","COHORT"))
dt[,ALLELE_ID:=as.numeric(ALLELE_ID)]

cn <- fread("/path/to/all_STR_loci.all_cohorts.copy_numbers.tsv", na.strings=".")
cn <- cn[,.(V1,V2,V4,V6,V7)]
setnames(cn, c("GENE","ALLELE_ID","ALLELE_UNITS","MOTIF","COHORT"))

# key + merge
setkey(dt, GENE, ALLELE_ID, COHORT)
setkey(cn, GENE, ALLELE_ID, COHORT)
m <- merge(dt, cn, by=c("GENE","ALLELE_ID","COHORT"), all.x=TRUE)

# Reference bp sizes (per locus)
ref_dt <- data.table(
  GENE = c("ABCD3","NOTCH2NLC","RILPL1","GIPC1"),
  REF_BP = c(30, 41, 35, 41)
)

# Pathogenic ranges in repeat units
path_dt <- data.table(
  GENE = c("NOTCH2NLC","ABCD3","RILPL1","GIPC1"),
  LO = c(66, 118, 120, 73),
  HI = c(517, 694, 197, 164)
)

# Join + convert to allele bp + units (period=3)
m <- merge(m, ref_dt, by = "GENE", all.x = TRUE)
m[, SAMPLE_SHORT := sub("^LRS0*([0-9]+).*", "L\\1", SAMPLE_ID)]

# Join pathogenic ranges onto main table
m <- merge(m, path_dt, by = "GENE", all.x = TRUE)

# Flag alleles inside pathogenic range
m[, IN_PATH := !is.na(LO) & ALLELE_UNITS >= LO]

set.seed(1)  # reproducible jitter

# numeric x position per cohort
m[, x_num := as.numeric(factor(COHORT))]

# apply the same jitter to all layers
jitter_width <- 0.15
m[, x_jit := x_num + runif(.N, -jitter_width, jitter_width)]

# Plot in units
p <- ggplot(m[ALLELE_UNITS > 0],
            aes(x = x_jit, y = ALLELE_UNITS)) +
  
  geom_violin(
    aes(x = x_num, group = interaction(GENE, COHORT)),
    trim = FALSE, scale = "width", alpha = 0.6
  ) +
  
  geom_point(size = 0.8, alpha = 0.6) +
  
  facet_wrap(~ GENE, nrow = 1) +
  scale_x_continuous(
    breaks = unique(m$x_num),
    labels = levels(factor(m$COHORT))
  ) +
  scale_y_log10() +
  labs(x = "Cohort", y = "Allele size (repeat units, log10 scale)") +
  theme_bw()

p <- p +
  geom_text_repel(
    data = m[IN_PATH == TRUE],
    aes(label = SAMPLE_SHORT),
    size = 3,
    colour = "red",
    direction = "y",
    box.padding = 0.15,
    point.padding = 0.05,
    force = 0.5
  )

# Add pathogenic range bands per facet (optional but useful)
p <- p +
  geom_rect(
    data = path_dt,
    inherit.aes = FALSE,
    aes(xmin = -Inf, xmax = Inf, ymin = LO, ymax = HI),
    fill = "red",
    alpha = 0.08
  )

print(p)