#!/usr/bin/env Rscript
.libPaths(c("/path/to/R_libs", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Inputs
# ============================================================

list_file <- "/path/to/all_svlog_database.txt"
panel_overlap_path <- "/path/to/sniffles_overlapping_panel.tsv"

# ============================================================
# Parameters (edit these to change rarity definitions)
# ============================================================

GNOMAD_AF_THRESHOLD <- 0.01

# Your current logic uses "alt allele count" = HET + 2*HOM
ONT1000G_MAX_ALT_ALLELES_RARE <- 2     # ~1% for n=149 (your original choice)
INTERNAL_MAX_ALT_ALLELES_RARE <- 1     # ~1% for n=61  (your original choice)

# ============================================================
# Read list of svlog database files
# ============================================================

paths <- readLines(list_file, warn = FALSE)
paths <- trimws(paths)
paths <- paths[nzchar(paths)]
paths <- paths[!grepl("^#", paths)]

if (length(paths) == 0) {
  stop("No input paths found in: ", list_file)
}

get_family_id <- function(p) {
  bn <- basename(p)
  sub("\\..*$", "", bn)
}

get_family_from_lrsid <- function(p) {
  bn <- basename(p)
  bn <- sub("\\..*$", "", bn)   # remove extension
  sub("-.*$", "", bn)           # keep part before first dash
}

# ============================================================
# Read + bind all svlog database tables
# ============================================================

dt_list <- vector("list", length(paths))
names(dt_list) <- paths

for (i in seq_along(paths)) {
  p <- paths[i]
  
  if (!file.exists(p)) {
    warning("Missing file, skipping: ", p)
    next
  }
  
  fam <- get_family_id(p)
  
  d <- fread(p, sep = "\t", na.strings = c("", "NA", "NaN"))
  d[, family_id := fam]
  d[, source_file := p]
  
  dt_list[[i]] <- d
}

dt_list <- Filter(Negate(is.null), dt_list)

if (length(dt_list) == 0) {
  stop("No files were successfully read.")
}

merged <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# ============================================================
# Read + bind all svlog *static* tables
# ============================================================

static_list_file <- "/path/to/all_sv_annotations.txt"

static_paths <- readLines(static_list_file, warn = FALSE)
static_paths <- trimws(static_paths)
static_paths <- static_paths[nzchar(static_paths)]
static_paths <- static_paths[!grepl("^#", static_paths)]

if (length(static_paths) == 0) {
  stop("No input paths found in: ", static_list_file)
}

static_dt_list <- vector("list", length(static_paths))
names(static_dt_list) <- static_paths

for (i in seq_along(static_paths)) {
  p <- static_paths[i]
  
  if (!file.exists(p)) {
    warning("Missing static file, skipping: ", p)
    next
  }
  
  fam <- get_family_from_lrsid(p)
  
  d <- fread(p, sep = "\t", na.strings = c("", "NA", "NaN"))
  d[, family_id := fam]
  d[, source_file := p]
  
  static_dt_list[[i]] <- d
}

static_dt_list <- Filter(Negate(is.null), static_dt_list)

if (length(static_dt_list) == 0) {
  stop("No static files were successfully read.")
}

merged_static <- rbindlist(static_dt_list, use.names = TRUE, fill = TRUE)

# ============================================================
# Bring in panel overlap table and merge with svlog data
# ============================================================

panel_overlap <- fread(panel_overlap_path, sep = "\t", header = TRUE)
panel_overlap[, family_id := sub("-.*$", "", sample_id)]

# Make names consistent
setnames(panel_overlap, "variant_id", "ID")

# Remove positional columns from svlog table (panel_overlap provides CHROM/START/END)
merged[, c("CHROM", "START", "END") := NULL]

# Set keys for joining
setkey(merged, family_id, ID)
setkey(panel_overlap, family_id, ID)

merged_with_panel <- merge(
  merged,
  panel_overlap,
  by = c("family_id", "ID")
)

# Define a single base table of the panel-overlapping variants (one row per family+variant)
base_variants <- unique(
  merged_with_panel[, .(family_id, ID, panel, sample_id, CHROM, START, END)]
)

# Reciprocal overlap score used for ranking best gnomAD match
merged_with_panel[, recip_score := pmin(OVL_Q, OVL_T)]


wide <- dcast(
  unique(merged_with_panel[, .(family_id, ID, svlog_id, SRC, STRING)]),
  family_id + ID + svlog_id ~ SRC,
  value.var = "STRING",
  fun.aggregate = function(x) paste(unique(na.omit(x)), collapse = ",")
)

# ============================================================
# Pick "best" match per source (gnomAD/ONT1000G/Internal)
# ============================================================

best_gnomad <- merged_with_panel[
  SRC == "gnomAD",
  .SD[
    order(-recip_score, -OVL_T, -OVL_Q, DIST)
  ][1],
  by = .(family_id, ID)
]

best_gnomad[, gnomad_rare := !is.na(AF) & AF < GNOMAD_AF_THRESHOLD]

best_ont <- merged_with_panel[
  SRC == "ONT1000G",
  .SD[1],
  by = .(family_id, ID)
]

best_ont[, ont1000g_carriers := HET + HOM]
best_ont[, ont1000g_rare := (HET + 2*HOM) <= ONT1000G_MAX_ALT_ALLELES_RARE]

best_internal <- merged_with_panel[
  SRC == "InternalCohort",
  .SD[1],
  by = .(family_id, ID)
]

best_internal[, internal_carriers := HET + HOM]
best_internal[, internal_rare := (HET + 2*HOM) <= INTERNAL_MAX_ALT_ALLELES_RARE]

# ============================================================
# Build per-variant summary table
# ============================================================

summary_dt <- merge(
  base_variants,
  best_gnomad[, .(family_id, ID, gnomad_AF = AF, gnomad_rare)],
  by = c("family_id", "ID"),
  all.x = TRUE
)

summary_dt <- merge(
  summary_dt,
  best_ont[, .(family_id, ID, ont1000g_carriers, ont1000g_rare)],
  by = c("family_id", "ID"),
  all.x = TRUE
)

summary_dt <- merge(
  summary_dt,
  best_internal[, .(family_id, ID, internal_carriers, internal_rare)],
  by = c("family_id", "ID"),
  all.x = TRUE
)

summary_dt <- merge(
  summary_dt,
  merged_static[,.(ID,SVLOG_ID,GENE_ID,GENE_SYMBOL,VEP_CONSEQUENCE,family_id)],
  by = c("family_id", "ID"),
  all.x = TRUE
)

summary_dt[, gnomad_is_rare     := (gnomad_rare %in% TRUE) | is.na(gnomad_AF)]
summary_dt[, ont1000g_is_rare   := (ont1000g_rare %in% TRUE) | is.na(ont1000g_carriers)]
summary_dt[, internal_is_rare   := (internal_rare %in% TRUE) | is.na(internal_carriers)]

summary_dt[,panel:=NULL]

summary_dt <- unique(summary_dt)

# ============================================================
# Summarise per family
# ============================================================

sample_summary <- summary_dt[
  ,
  .(
    n_panel_sv = .N,
    
    # Keep original per-source metrics
    n_gnomad_rare   = sum(gnomad_is_rare, na.rm = TRUE),
    n_gnomad_absent = sum(is.na(gnomad_AF)),
    
    n_ont1000g_rare = sum(ont1000g_is_rare, na.rm = TRUE),
    n_internal_rare = sum(internal_is_rare, na.rm = TRUE),
    
    # NEW: rare across all sources (absence counts as rare)
    n_rare_all = sum(
      gnomad_is_rare & ont1000g_is_rare & internal_is_rare,
      na.rm = TRUE
    ),
    
    # Gene-level summaries
    n_gene_sv = sum(!is.na(GENE_SYMBOL) & nzchar(GENE_SYMBOL), na.rm = TRUE),
    
    n_gene_gnomad_rare = sum(
      (!is.na(GENE_SYMBOL) & nzchar(GENE_SYMBOL)) &
        (gnomad_is_rare),
      na.rm = TRUE
    ),
    
    n_gene_ont1000g_rare = sum(
      (!is.na(GENE_SYMBOL) & nzchar(GENE_SYMBOL)) &
        (ont1000g_is_rare),
      na.rm = TRUE
    ),
    
    n_gene_internal_rare = sum(
      (!is.na(GENE_SYMBOL) & nzchar(GENE_SYMBOL)) &
        (internal_is_rare),
      na.rm = TRUE
    ),
    
    # NEW: gene + rare in all sources (absence counts as rare)
    n_gene_rare_all = sum(
      (!is.na(GENE_SYMBOL) & nzchar(GENE_SYMBOL)) &
        gnomad_is_rare & ont1000g_is_rare & internal_is_rare,
      na.rm = TRUE
    )
  ),
  by = family_id
]

#sample_summary[, family_id := sub("^LRS0*", "L", family_id)]

setcolorder(sample_summary, c(
  "family_id",
  "n_panel_sv",
  "n_gene_sv",
  "n_gnomad_absent",
  "n_gnomad_rare",
  "n_ont1000g_rare",
  "n_internal_rare",
  "n_gene_gnomad_rare",
  "n_gene_ont1000g_rare",
  "n_gene_internal_rare",
  "n_rare_all",
  "n_gene_rare_all"
))

sample_summary[, family_id := sub("RS00", "", family_id)]

fwrite(sample_summary[family_id!="L152"], "/path/to/myopathy_svs_sample_summary.tsv", sep = "\t")