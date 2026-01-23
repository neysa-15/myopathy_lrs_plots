setwd("/g/data/kr68/neysa/r_plotting/rcode_per_fig")

# ------------------------------------------------------------------------------
#                 Libraries
# ------------------------------------------------------------------------------
.libPaths(c("/g/data/kr68/andre/R_libs"))
# install.packages(c("data.table","ggplot2","patchwork","scales","zoo"))

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

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


# ------------------------------------------------------------------------------
#                 Figure 3b, h and k
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

source("fig2.R") # - taking function from fig2 - if not doing this then need to copy facet plot too
plot_4q <- combined_dt[
  pLAM_mapped == TRUE &
    grepl("^chr4:", GenomeCoords) &
    grepl("Complete 4qA|Partial distal 4qA", ReadLabel)
]

# Fig3c - FSHD1 Positive - Double contraction
sample_id <- "R240183"
plot_fshd_5mC_vs_d4z4copies_4qA(sample_id, plot_4q)

# Fig3k - Previously undiagnosed FSHD1 case
sample_id <- "CF2608"
plot_fshd_5mC_vs_d4z4copies_4qA(sample_id, plot_4q)

# Fig3h - Previously false positive FSHD1
sample_id <- "GL2106"
plot_fshd_5mC_vs_d4z4copies_facet(sample_id)


# ------------------------------------------------------------------------------
#                 Figure 3c and 3l - pileup plot
# ------------------------------------------------------------------------------
pileup_plot <- function(sample_id, result_dir, reference_feature = "pLAM", complete_only = TRUE) {

  # Set reference configuration
  if (reference_feature == "pLAM") {
    if (complete_only == TRUE) {
      ref_label_filter <- c("Complete 4qA")
    } else {
      ref_label_filter <- c("Partial distal 4qA", "Complete 4qA")
    }
    plot_label <- "start"
  } else if (reference_feature == "p13-E11") {
    if (complete_only == TRUE) {
      ref_label_filter <- c("Complete 4qA")
    } else {
      ref_label_filter <- c("Partial proximal Unclassified", "Complete 4qA")
    }
    plot_label <- "end"
  } else {
    stop("Unsupported reference feature")
  }
  
  # Load input files
  tsv <- fread(sprintf('%s/%s/%s_mapped_features_summary.tsv', result_dir, sample_id, sample_id), header = TRUE)
  fshd2 <- fread(sprintf('%s/%s/%s_all_features.bed', result_dir, sample_id, sample_id), header = FALSE)
  
  # Filter reads of interest
  filtered <- tsv[
    ReadLabel %in% ref_label_filter &
      duplex == FALSE &
      startsWith(GenomeCoords, "chr4")
  ]
  plus_strand_ids <- filtered[strand == "+", ReadID]
  minus_strand_ids <- filtered[strand == "-", ReadID]
  
  # Process plus strand
  fshd2_plus <- fshd2[V1 %in% plus_strand_ids]
  ref_plus <- fshd2_plus[V4 == reference_feature,
                         .(ref_coord = if (reference_feature == "pLAM") min(V2) else max(V3)),
                         by = V1]
  fshd2_plus <- merge(fshd2_plus, ref_plus, by = "V1")
  fshd2_plus[, `:=`(
    rel_start = V2 - ref_coord,
    rel_end   = V3 - ref_coord,
    strand = "+"
  )]
  
  # Process minus strand
  fshd2_minus <- fshd2[V1 %in% minus_strand_ids]
  read_lengths <- fshd2_minus[, .(read_start = min(V2), read_end = max(V3)), by = V1]
  fshd2_minus <- merge(fshd2_minus, read_lengths, by = "V1")
  ref_minus <- fshd2_minus[V4 == reference_feature,
                           .(ref_coord = if (reference_feature == "pLAM") max(V3) else min(V2)),
                           by = V1]
  fshd2_minus <- merge(fshd2_minus, ref_minus, by = "V1")
  fshd2_minus[, `:=`(
    rel_start = -(V3 - ref_coord),
    rel_end   = -(V2 - ref_coord)
  )]
  fshd2_minus[, `:=`(
    rel_start = pmin(rel_start, rel_end),
    rel_end   = pmax(rel_start, rel_end),
    strand = "-"
  )]
  
  # Combine
  combined <- rbind(fshd2_plus, fshd2_minus, fill = TRUE)
  
  # Relabel d4z4 units by proximity to reference point, then assign parity
  if (reference_feature == "pLAM") {
    combined[grepl("^d4z4_\\d+$", V4), d4z4_rank := frank(-rel_start, ties.method = "first"), by = V1]
  } else {
    combined[grepl("^d4z4_\\d+$", V4), d4z4_rank := frank(rel_start, ties.method = "first"), by = V1]
  }
  combined[!is.na(d4z4_rank), V4 := paste0("d4z4_", fifelse(d4z4_rank %% 2 == 0, "even", "odd"))]
  
  # Standardise labels
  combined[V4 == "4qA_marker_new", V4 := "4qA_marker1"]
  combined[V4 == "4qA_probe", V4 := "4qA_marker2"]
  
  # Order reads
  filtered[, label_priority := fifelse(ReadLabel == "Complete 4qA", 1L, 2L)]
  d4z4_counts <- combined[V4 %in% c("d4z4_odd", "d4z4_even"), .N, by = V1]
  d4z4_counts <- merge(d4z4_counts, filtered[, .(V1 = ReadID, label_priority, ReadLabel)], by = "V1", all.x = TRUE)
  ordered_ids <- d4z4_counts[order(-N, ReadLabel), V1]
  combined[, V1 := factor(V1, levels = ordered_ids)]
  
  # Compute read span
  read_lengths <- filtered[, .(V1 = ReadID, ReadLength)]
  ref_pos <- combined[V4 == reference_feature,
                      .(ref_coord = if (reference_feature == "pLAM" & strand[1] == "+") min(V2)
                        else if (reference_feature == "pLAM" & strand[1] == "-") max(V3)
                        else if (reference_feature == "p13-E11" & strand[1] == "+") max(V3)
                        else min(V2)),
                      by = .(V1, strand)]
  read_span <- data.table(V1 = ordered_ids, read_index = 1:length(ordered_ids))
  read_span <- merge(read_span, read_lengths)
  read_span <- merge(read_span, ref_pos, by = "V1")
  read_span[, `:=`(
    rel_start = ifelse(strand == "+", -ref_coord, ref_coord - ReadLength),
    rel_end   = ifelse(strand == "+", ReadLength - ref_coord, ref_coord)
  )]
  
  # Overlay/underlay
  overlay <- combined[V4 %in% c("pLAM", "SSLP", "p13-E11")]
  underlay <- combined[!V4 %in% c("pLAM", "SSLP", "p13-E11")]
  for (dt in list(combined, overlay, underlay)) {
    set(dt, j = "read_index", value = read_span[dt, on = "V1", read_index])
  }
  
  # Plot
  feature_colours <- c(
    "d4z4_chr4_proximal" = "#4D4D4D",
    "SSLP"               = "#FFD700",
    "p13-E11"            = "red",
    "d4z4_odd"           = "#B2E2B2",
    "d4z4_even"          = "#1F78B4",
    "pLAM"               = "magenta",
    "4qA_marker1"        = "orange",
    "4qA_marker2"        = "#00CED1"
  )
  
  ggplot() +
    geom_segment(data = read_span,
                 aes(x = rel_start, xend = rel_end, y = read_index, yend = read_index),
                 colour = "lightgray", size = 3) +
    geom_segment(data = underlay,
                 aes(x = rel_start, xend = rel_end, y = read_index, yend = read_index, colour = V4),
                 size = 3) +
    geom_segment(data = overlay,
                 aes(x = rel_start, xend = rel_end, y = read_index, yend = read_index, colour = V4),
                 size = 2.5) +
    scale_colour_manual(values = feature_colours) +
    theme_bw(base_size = 15) +
    labs(x = sprintf("Position relative to %s %s", reference_feature,plot_label), y = "Reads", colour = "Feature", title = sample_id) +
    theme(
      strip.text.y = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# 3c
pileup_plot("R240183", base_dir)
# 3l
pileup_plot("GUAT0705", base_dir)

# ------------------------------------------------------------------------------
#                 Figure 3e - FSHD2 assembly methylation of 4qA allele
# ------------------------------------------------------------------------------

plot_asm_meth <- function(sample_id, summary_tsv, annot_bed = NULL,
                          contig = NULL, smooth_n = 5, x_units = c("kb","Mb"),
                          reference_feature = c("pLAM","p13-E11")) {
  
  reference_feature <- match.arg(reference_feature)
  x_units <- match.arg(x_units)
  
  # 1. Load Methylation Data
  dt <- fread(summary_tsv)
  sel_contig <- if (is.null(contig)) dt$contig[1] else contig
  if (!sel_contig %in% dt$contig) stop("Contig not found: ", sel_contig)
  
  dt <- dt[contig == sel_contig]
  dt[, x := (start + end)/2]
  conv <- if (x_units == "kb") 1e3 else 1e6
  dt[, x_scaled := x/conv]
  
  # 2. Load and Prepare Annotations
  feat <- data.table()
  strand_val <- "+" # Default
  # if (!is.null(annot_bed) && file.exists(annot_bed)) {
  #   feat <- fread(annot_bed, header = FALSE)
  #   strand <- feat$V6[feat$V4 == "pLAM"][1]
  # } else {
  #   strand <- "+"
  # }
  
  if (!is.null(annot_bed) && file.exists(annot_bed)) {
    feat <- fread(annot_bed, header = FALSE)
    if (ncol(feat) >= 6) {
      setnames(feat, 1:6, c("contig","start","end","feature","score","strand"))
    } else if (ncol(feat) >= 4) {
      setnames(feat, 1:4, c("contig","start","end","feature"))
      feat[, strand := "+"]
    }
    
    feat <- feat[contig == sel_contig]
    
    if (nrow(feat) > 0) {
      # Harmonise labels
      feat[feature == "4qA_marker_new", feature := "4qA_marker1"]
      feat[feature == "4qA_probe",      feature := "4qA_marker2"]
      
      # Detect strand from pLAM if available
      if ("pLAM" %in% feat$feature) {
        strand_val <- feat[feature == "pLAM", strand][1]
      }
    }
  }
  
  # 3. Handle Coordinate Flipping (Keep axis ascending, flip data)
  xr <- range(dt$x_scaled, na.rm = TRUE)
  
  if (strand_val == "-") {
    # Transform methylation points: x' = (max + min) - x
    dt[, x_scaled := (xr[2] + xr[1]) - x_scaled]
    
    # Transform annotation coordinates
    if (nrow(feat) > 0) {
      feat[, `:=`(x0 = (xr[2] + xr[1]) - (end/conv), 
                  x1 = (xr[2] + xr[1]) - (start/conv))]
    }
  } else {
    if (nrow(feat) > 0) {
      feat[, `:=`(x0 = start/conv, x1 = end/conv)]
    }
  }
  
  # 4. Smoothing
  if (smooth_n > 1 && nrow(dt) >= smooth_n) {
    dt[, met_sm := frollmean(window_methy_avg, n = smooth_n, align = "center")]
  } else {
    dt[, met_sm := window_methy_avg]
  }
  
  # 5. D4Z4 Relabeling Logic
  if (nrow(feat) > 0 && any(feat$feature == reference_feature)) {
    ref_coord <- switch(
      reference_feature,
      "pLAM"   = min(feat[feature == "pLAM", start]),
      "p13-E11"= max(feat[feature == "p13-E11", end])
    )
    
    feat[, `:=`(rel_start = start - ref_coord, rel_end = end - ref_coord)]
    d4_idx <- grepl("^d4z4_\\d+$", feat$feature)
    
    if (any(d4_idx)) {
      if (reference_feature == "pLAM") {
        feat[d4_idx, d4_rank := frank(-pmin(rel_start, rel_end), ties.method = "first")]
      } else {
        feat[d4_idx, d4_rank := frank( pmin(rel_start, rel_end), ties.method = "first")]
      }
      feat[!is.na(d4_rank), feature := paste0("d4z4_", ifelse(d4_rank %% 2 == 0, "even", "odd"))]
    }
  }
  
  # 6. Plotting
  xlab <- paste0("Position on ", sel_contig, " (", x_units, ") (Strand: ", strand_val, ")")
  thr <- 50
  
  p_met <- ggplot(dt, aes(x = x_scaled)) +
    geom_hline(yintercept = thr, linetype = "dashed", colour = "grey70") +
    geom_line(aes(y = 100 * met_sm), colour = "grey60", linewidth = 0.35) +
    geom_line(aes(y = ifelse(100 * met_sm >= thr, 100 * met_sm, NA_real_)),
              colour = "firebrick", linewidth = 0.7, na.rm = TRUE) +
    geom_line(aes(y = ifelse(100 * met_sm <  thr, 100 * met_sm, NA_real_)),
              colour = "steelblue", linewidth = 0.7, na.rm = TRUE) +
    scale_x_continuous(limits = xr, expand = c(0, 0)) +
    labs(title = sample_id, x = xlab, y = "5mC (%)") +
    coord_cartesian(ylim = c(0, 100)) +
    theme_bw(12) +
    theme(panel.grid = element_blank())
  
  if (nrow(feat) > 0) {
    feature_colours <- c(
      "d4z4_chr4_proximal" = "#4D4D4D", "SSLP" = "#FFD700", "p13-E11" = "red",
      "d4z4_odd" = "#B2E2B2", "d4z4_even" = "#1F78B4", "pLAM" = "magenta",
      "4qA_marker1" = "orange", "4qA_marker2" = "#00CED1",
      "4qB_marker1" = "#A52A2A", "4qB_marker2" = "#6A5ACD"
    )
    
    overlay_names <- c("pLAM","SSLP","p13-E11")
    overlay  <- feat[feature %in% overlay_names]
    underlay <- feat[!feature %in% overlay_names]
    
    p_feat <- ggplot() +
      geom_segment(aes(x = xr[1], xend = xr[2], y = 0, yend = 0),
                   colour = "lightgray", linewidth = 12) +
      geom_segment(data = underlay,
                   aes(x = x0, xend = x1, y = 0, yend = 0, colour = feature),
                   linewidth = 12, lineend = "butt") +
      geom_segment(data = overlay,
                   aes(x = x0, xend = x1, y = 0, yend = 0, colour = feature),
                   linewidth = 12, lineend = "butt") +
      scale_x_continuous(limits = xr, expand = c(0, 0)) +
      scale_colour_manual(values = feature_colours, breaks = names(feature_colours)) +
      labs(x = xlab, y = NULL, colour = NULL) +
      theme_bw(12) +
      theme(panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "bottom")
    
    return(p_met / p_feat + plot_layout(heights = c(1, 0.25)))
  }
  
  p_met
}

sample_id <- "DL1104"
phase     <- "hap1"
base_dir  <- "/g/data/kr68/fshd/results_phased_assembly"

summary_tsv <- sprintf("%s/%s/%s_%s/%s_realigned_modfreq_cov_summary.tsv",
                       base_dir, sample_id, sample_id, phase, sample_id)
annot_bed   <- sprintf("%s/%s/%s_%s/%s_all_features.bed",
                       base_dir, sample_id, sample_id, phase, sample_id)

panel <- plot_asm_meth(
  sample_id   = sample_id,
  summary_tsv = summary_tsv,
  annot_bed   = annot_bed,
  contig      = "ptg000001l",   # or NULL to use the first in the file
  smooth_n    = 7,
  x_units     = "kb",
  reference_feature = "pLAM"
)

# strip any legends from the left panel explicitly (belt & braces)
panel_noleg <- panel +
  theme(legend.position = "none") +
  guides(colour = "none", fill = "none", linetype = "none",
         shape = "none", size = "none", alpha = "none")

panel_noleg

# ------------------------------------------------------------------------------
#                 Figure 3f - FSHD2 distal methylation levels
# ------------------------------------------------------------------------------
sample_id <- c("DL1104","R240059", "R250109")
fshd2_meth <- combined_dt[!is.na(Haplotype) & MappedEstimatedCopies>=2 & Sample %in% sample_id]

# Build the boxplot as its own object
ggplot(fshd2_meth, aes(x = Haplotype, y = pLAM_Methylation_Percentage)) +
  geom_boxplot() +
  facet_wrap(~Sample) +
  labs(x = "Haplotype", y = "Distal D4Z4 5mC (%)") +
  geom_hline(yintercept = 50, linetype = "dashed") +
  theme_bw(base_size = 25) +
  theme(                 # tall-ish to better match the left panel
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 18)
  )


# ------------------------------------------------------------------------------
#                 Figure 3i - False positive FSHD1 restriction sites
# ------------------------------------------------------------------------------
sample_id <- "GL2106"
sample_dt <- combined_dt[combined_dt$Sample == "GL2106"]
max_number_copies <- max(as.integer(sample_dt$MappedEstimatedCopies), na.rm = TRUE)

plot_fshd_read <- function(sample_id,
                           read_id,
                           reference_feature = c("pLAM","p13-E11"),
                           base_dir = "/g/data/kr68/fshd/results_mapq0",
                           track_layout = c("rows","overlay"),
                           use_file_counts = TRUE,
                           min_x = NULL, max_x = NULL) {
  reference_feature <- match.arg(reference_feature)
  track_layout <- match.arg(track_layout)
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
  })
  
  # ---------- Files ----------
  tsv_path  <- file.path(base_dir, sample_id, sprintf("%s_mapped_features_summary.tsv", sample_id))
  bed_path  <- file.path(base_dir, sample_id, sprintf("%s_all_features.bed", sample_id))
  xapi_path <- file.path(base_dir, sample_id, sprintf("%s_xapi_sites.bed", sample_id))
  blni_path <- file.path(base_dir, sample_id, sprintf("%s_blni_sites.bed", sample_id))
  
  # ---------- Load ----------
  tsv <- fread(tsv_path)
  bed <- fread(bed_path, header = FALSE)
  setnames(bed, c("V1","V2","V3","V4"), c("ReadID","Start","End","Feature"))
  
  # Meta and subset
  meta <- tsv[ReadID == read_id]
  if (nrow(meta) == 0L) stop("read_id not found in summary TSV for this sample")
  dt <- bed[ReadID == read_id]
  if (nrow(dt) == 0L) stop("No feature rows for this read in BED")
  if (!reference_feature %in% dt$Feature) {
    stop(sprintf("Reference feature '%s' not present on this read", reference_feature))
  }
  
  strand <- if ("strand" %in% names(meta)) meta$strand[1] else "+"
  
  # ---------- Anchor features so reference = 0 ----------
  if (reference_feature == "pLAM") {
    ref_coord <- if (strand == "+") min(dt[Feature == "pLAM", Start]) else max(dt[Feature == "pLAM", End])
    dt[, `:=`(
      rel_start = if (strand == "+") Start - ref_coord else -(End - ref_coord),
      rel_end   = if (strand == "+") End   - ref_coord else -(Start - ref_coord)
    )]
  } else {
    ref_coord <- if (strand == "+") max(dt[Feature == "p13-E11", End]) else min(dt[Feature == "p13-E11", Start])
    dt[, `:=`(
      rel_start = if (strand == "+") Start - ref_coord else -(End - ref_coord),
      rel_end   = if (strand == "+") End   - ref_coord else -(Start - ref_coord)
    )]
  }
  dt[, `:=`(rel_min = pmin(rel_start, rel_end), rel_max = pmax(rel_start, rel_end))]
  
  # ---------- Relabel D4Z4 alternating odd/even ----------
  d4z4_idx <- grepl("^d4z4_\\d+$", dt$Feature)
  if (any(d4z4_idx)) {
    if (reference_feature == "pLAM") {
      dt[d4z4_idx, d4z4_rank := frank(-rel_min, ties.method = "first")]
    } else {
      dt[d4z4_idx, d4z4_rank := frank(rel_min, ties.method = "first")]
    }
    dt[!is.na(d4z4_rank), Feature := paste0("d4z4_", ifelse(d4z4_rank %% 2 == 0, "even", "odd"))]
  }
  
  # ---------- Harmonise marker names ----------
  dt[Feature == "4qA_marker_new", Feature := "4qA_marker1"]
  dt[Feature == "4qA_probe",      Feature := "4qA_marker2"]
  dt[Feature == "4qB_probe",      Feature := "4qB_marker1"]
  dt[Feature == "4qA_down",       Feature := "4qB_marker2"]
  
  # ---------- Site BEDs ----------
  read_site_bed <- function(path, read_id) {
    if (!file.exists(path)) return(data.table())
    s <- fread(path, header = FALSE)
    if (ncol(s) < 3) stop("Site BED must have at least 3 cols: ReadID start end [strand]")
    s <- s[, 1:min(4, ncol(s)), with = FALSE]
    nm <- c("ReadID","Start","End","Strand")[seq_len(ncol(s))]
    setnames(s, nm)
    s <- s[ReadID == read_id]
    if (!nrow(s)) return(s)
    unique(s[, .(ReadID, Start, End)])
  }
  
  xapi_rows <- read_site_bed(xapi_path, read_id)
  blni_rows <- read_site_bed(blni_path, read_id)
  
  # Fallbacks from main BED if needed
  if (!nrow(xapi_rows) && any(grepl("XapI", dt$Feature))) {
    xapi_rows <- unique(dt[grepl("XapI", Feature), .(ReadID, Start, End)])
  }
  if (!nrow(blni_rows) && any(grepl("BlnI", dt$Feature))) {
    blni_rows <- unique(dt[grepl("BlnI", Feature), .(ReadID, Start, End)])
  }
  
  # Anchor sites
  anchor_to_ref <- function(s) {
    if (!nrow(s)) return(s[, `:=`(rel_min = numeric(), rel_max = numeric())])
    if (reference_feature == "pLAM") {
      s[, `:=`(
        rel_start = if (strand == "+") Start - ref_coord else -(End - ref_coord),
        rel_end   = if (strand == "+") End   - ref_coord else -(Start - ref_coord)
      )]
    } else {
      s[, `:=`(
        rel_start = if (strand == "+") Start - ref_coord else -(End - ref_coord),
        rel_end   = if (strand == "+") End   - ref_coord else -(Start - ref_coord)
      )]
    }
    s[, `:=`(rel_min = pmin(rel_start, rel_end), rel_max = pmax(rel_start, rel_end))]
  }
  
  xapi_rows <- anchor_to_ref(xapi_rows)
  blni_rows <- anchor_to_ref(blni_rows)
  
  xapi_pts <- if (nrow(xapi_rows)) xapi_rows[, .(x = (rel_min + rel_max) / 2)] else data.table(x = numeric(0))
  blni_pts <- if (nrow(blni_rows)) blni_rows[, .(x = (rel_min + rel_max) / 2)] else data.table(x = numeric(0))
  
  # ---------- Colours ----------
  feature_cols <- c(
    "d4z4_odd"           = "#B2E2B2",
    "d4z4_even"          = "#1F78B4",
    "SSLP"               = "#FFD700",
    "pLAM"               = "magenta",
    "p13-E11"            = "red",
    "4qA_marker1"        = "orange",
    "4qA_marker2"        = "#00CED1",
    "4qB_marker1"        = "#A52A2A",
    "4qB_marker2"        = "#6A5ACD",
    "d4z4_chr4_proximal" = "#4D4D4D"
  )
  site_cols <- c("XapI" = "#2b8cbe", "BlnI" = "#e34a33")
  
  # Remove site rows from coloured segments
  seg_dt <- dt
  if (nrow(xapi_rows)) seg_dt <- seg_dt[!(Start %in% xapi_rows$Start & End %in% xapi_rows$End)]
  if (nrow(blni_rows)) seg_dt <- seg_dt[!(Start %in% blni_rows$Start & End %in% blni_rows$End)]
  
  # Full-span (before windowing)
  read_span <- seg_dt[, .(x = min(rel_min, na.rm = TRUE), xend = max(rel_max, na.rm = TRUE))]
  
  # ---------- Optional x-window ----------
  xr_all <- c(read_span$x[1], read_span$xend[1])
  xr_req <- c(ifelse(is.null(min_x), xr_all[1], min_x),
              ifelse(is.null(max_x), xr_all[2], max_x))
  if (xr_req[1] >= xr_req[2]) stop("min_x must be < max_x")
  
  # Clip data to window (plotting speed; backbone handled per-branch)
  seg_dt   <- seg_dt  [rel_max >= xr_req[1] & rel_min <= xr_req[2]]
  if (nrow(xapi_pts)) xapi_pts <- xapi_pts[x >= xr_req[1] & x <= xr_req[2]]
  if (nrow(blni_pts)) blni_pts <- blni_pts[x >= xr_req[1] & x <= xr_req[2]]
  
  # Counts for labels
  file_x_cnt <- nrow(xapi_pts)
  file_b_cnt <- nrow(blni_pts)
  tsv_x_cnt  <- meta$XapI_Sensitive_Repeats
  tsv_b_cnt  <- meta$BlnI_Sensitive_Repeats
  show_x_cnt <- if (use_file_counts && file_x_cnt > 0) file_x_cnt else ifelse(length(tsv_x_cnt) && !is.na(tsv_x_cnt), tsv_x_cnt, NA)
  show_b_cnt <- if (use_file_counts && file_b_cnt > 0) file_b_cnt else ifelse(length(tsv_b_cnt) && !is.na(tsv_b_cnt), tsv_b_cnt, NA)
  
  # ---------- Plot ----------
  if (track_layout == "overlay") {
    site_dt <- rbindlist(list(
      if (nrow(xapi_pts)) data.table(x = xapi_pts$x, y0 = -0.5, y1 = 0.5, Site = "XapI") else NULL,
      if (nrow(blni_pts)) data.table(x = blni_pts$x, y0 = -0.5, y1 = 0.5, Site = "BlnI") else NULL
    ), use.names = TRUE, fill = TRUE)
    
    p <- ggplot() +
      geom_segment(
        data = read_span, aes(x = x, xend = xend, y = 0, yend = 0),
        colour = "grey80", linewidth = 4
      ) +
      geom_segment(
        data = seg_dt, aes(x = rel_min, xend = rel_max, y = 0, yend = 0, colour = Feature),
        linewidth = 4, lineend = "butt"
      ) +
      scale_colour_manual(values = feature_cols, na.translate = TRUE) +
      coord_cartesian(xlim = xr_req, ylim = c(-1, 1)) +
      labs(
        x = sprintf("Position relative to %s (%s = 0)", reference_feature,
                    ifelse(reference_feature == "pLAM", "start", "end")),
        y = NULL, colour = "Feature",
        title = sprintf("%s - read %s", sample_id, read_id)
      ) +
      theme_bw(base_size = 13) +
      theme(panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "bottom",
            aspect.ratio = 1/3) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4)
    
    if (nrow(site_dt)) {
      p <- p +
        geom_segment(
          data = site_dt,
          aes(x = x, xend = x, y = y0, yend = y1, colour = Site),
          inherit.aes = FALSE, linewidth = 0.9, lineend = "butt", alpha = 0.95
        ) +
        scale_discrete_manual(aesthetics = "colour",
                              values = c(feature_cols, site_cols),
                              breaks = names(feature_cols),
                              name = "Feature")
      if (is.finite(suppressWarnings(as.numeric(show_x_cnt))) ||
          is.finite(suppressWarnings(as.numeric(show_b_cnt)))) {
        lab <- sprintf("XapI: %s   BlnI: %s",
                       ifelse(!is.na(show_x_cnt), show_x_cnt, "NA"),
                       ifelse(!is.na(show_b_cnt), show_b_cnt, "NA"))
        p <- p + annotate("label", x = Inf, y = 0.9, hjust = 1.05, vjust = 0.5,
                          label = lab, size = 3, label.size = 0.2, fill = "white")
      }
    }
    
  } else {
    # rows layout
    feat_panel <- seg_dt[, .(Panel = "Features", rel_min, rel_max, Feature)]
    xapi_panel <- if (nrow(xapi_pts)) data.table(Panel = "XapI sites", x = xapi_pts$x) else data.table(Panel = "XapI sites", x = numeric(0))
    blni_panel <- if (nrow(blni_pts)) data.table(Panel = "BlnI sites", x = blni_pts$x) else data.table(Panel = "BlnI sites", x = numeric(0))
    
    sublab <- sprintf("XapI count: %s | BlnI count: %s",
                      ifelse(!is.na(show_x_cnt), show_x_cnt, "NA"),
                      ifelse(!is.na(show_b_cnt), show_b_cnt, "NA"))
    
    y0 <- -0.4; y1 <- 0.4
    
    # Facet order & backbone in window
    levs <- c("Features", "XapI sites", "BlnI sites")
    feat_panel[, Panel := factor(Panel, levels = levs)]
    xapi_panel[, Panel := factor(Panel, levels = levs)]
    blni_panel[, Panel := factor(Panel, levels = levs)]
    backbone_dt <- data.table(Panel = factor("Features", levels = levs),
                              x = xr_req[1], xend = xr_req[2], y = 0, yend = 0)
    
    p <- ggplot() +
      geom_segment(
        data = backbone_dt,
        aes(x = x, xend = xend, y = y, yend = y),
        colour = "grey80", linewidth = 4
      ) +
      geom_segment(
        data = feat_panel,
        aes(x = rel_min, xend = rel_max, y = 0, yend = 0, colour = Feature),
        linewidth = 4, lineend = "butt"
      ) +
      geom_segment(
        data = xapi_panel,
        aes(x = x, xend = x, y = y0, yend = y1),
        inherit.aes = FALSE, linewidth = 0.9, colour = "#2b8cbe", lineend = "butt"
      ) +
      geom_segment(
        data = blni_panel,
        aes(x = x, xend = x, y = y0, yend = y1),
        inherit.aes = FALSE, linewidth = 0.9, colour = "#e34a33", lineend = "butt"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
      scale_colour_manual(values = feature_cols, na.translate = TRUE) +
      scale_x_continuous(limits = xr_req, expand = c(0, 0)) +
      labs(
        x = sprintf("Position relative to %s (%s = 0)", reference_feature,
                    ifelse(reference_feature == "pLAM", "start", "end")),
        y = NULL, colour = "Feature",
        title = NULL
      ) +
      facet_grid(rows = vars(Panel), scales = "free_y", drop = FALSE) +
      theme_bw(base_size = 13) +
      theme(panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "none",
            strip.text.y = element_text(face = "bold"))
  }
  
  return(p)
}


# Separate rows for XapI/BlnI
read_example <- plot_fshd_read("GL2106", "260720b8-7c68-4003-bac2-bbad92bda96a",
                               reference_feature = "pLAM",
                               track_layout = "rows",
                               min_x = -42000, max_x = 10000)

# sample_id <- c("GL2106")
false_positive <- combined_dt[!is.na(Haplotype) & MappedEstimatedCopies>=2 & Sample %in% sample_id]

false_positive$Haplotype <- factor(false_positive$Haplotype,levels=c("4qA","4qB","10qA"))

false_positive_plot <- ggplot(false_positive, aes(x = pLAM_Methylation_Percentage, y = floor(MappedEstimatedCopies), color = ReadLabel)) +
  geom_point(size = 5, alpha = 0.6) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 50, linetype = "dashed", color = "black") +
  facet_wrap(~Sample)+
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, max_number_copies + 5)) +
  scale_color_manual(values = c(label_colours, "Complete 4qA" = "#1f77b4", "Partial distal 4qA" = "#ff7f0e","Complete 10qA" = "#1f77b4", "Partial distal 10qA" = "#ff7f0e","Partial distal 4qB" = "#ff7f0e")) +
  labs(x = "Distal D4Z4 5mC (%)", y = "Number of D4Z4 copies", color = "Read type") +
  facet_wrap(~ Haplotype, drop = FALSE) +
  theme_bw(base_size = 25) +
  theme(legend.position = "bottom",
        #aspect.ratio = 1,
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

no_strips <- theme(
  strip.background  = element_blank(),
  strip.background.x= element_blank(),
  strip.background.y= element_blank(),
  strip.text        = element_blank(),
  strip.text.x      = element_blank(),
  strip.text.y      = element_blank()
)

false_positive_plot <- false_positive_plot + no_strips
read_example        <- read_example        + no_strips
