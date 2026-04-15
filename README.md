# Myopathy ONT figures

This repository contains the code to generate figures in [Targeted long-read sequencing enables comprehensive analysis of the genetic and epigenetic landscape of inherited myopathies](https://www.medrxiv.org/content/10.64898/2025.12.06.25340828v1).

Most figures are generated directly from output files produced by [d4z4ling](https://github.com/neysa-15/d4z4ling), our D4Z4 repeat annotation and methylation analysis tool. Where additional upstream processing was required (e.g. STR genotyping, subsampling experiments, mtDNA deletion analysis), the relevant scripts and input formats are described below. For full details of all methods, please refer to the Methods section of the paper.

---

## Repository structure

```
myopathy_lrs_plots/
├── fig2.R                    # FSHD D4Z4 copy number and methylation plots + limit of detection
├── fig3.R                    # FSHD case study plots (pileup, methylation, assembly)
├── fig4/
│   ├── fig4b_genotype_STR.sh # LongTR STR genotyping script
│   └── fig4c_promoter_methylation.sh  # Bedtools promoter methylation script
├── fig5_mtDNA_analysis/      # mtDNA deletion analysis scripts
│   ├── map_sra_data_to_ref.sh # script to map SRA data to hg38
│   ├── extract_deletion_cigar.sh # script to extract cigar and run big_del_proportion.py calculating big deletion proportion
│   └── fig5i.R               # script to create the proportion comparison plot shown on fig5i
├── figS1.R                   # Cohort demographics heatmap
├── figS2.R                   # Sequencing coverage summary
├── figS6.R                   # D4Z4 methylation/copies vs age of onset correlations
├── figS7.R                   # Subsampling experiment (limit of detection)
├── figS8.R                   # OPDM STR allele size distributions
└── sup_table_3.R             # Supplementary Table 3 (SV annotations)
```

---

## Dependencies

### R libraries

All R scripts require the following libraries. Install them into your R library path before running:

- `data.table`
- `ggplot2`
- `tidyr`
- `patchwork`
- `scales`
- `dplyr`
- `ggrepel`

### Command-line tools

- `LongTR` v1.0 — STR genotyping (fig4b)
- `bedtools` v2.31.0 — promoter methylation (fig4c)
- `TRF` (Tandem Repeat Finder) — repeat copy number from ALT sequences (figS8)
- `samtools`, `minimap2`, `sratoolkit`, `mosdepth`, `python3` — mtDNA analysis (fig5)

---

## Generic inputs (shared across figures)

### R library path
Set `.libPaths(c("/path/to/your/r_library"))` at the top of each script.

### d4z4ling results directory
Most figures read directly from d4z4ling output. The expected directory structure is:

```
/path/to/d4z4ling/results_directory/
├── sample1/
│   ├── sample1_mapped_features_summary.tsv
│   ├── sample1_all_features.bed
│   └── ...
├── sample2/
│   └── ...
```

### Sample key
A two-column TSV (`/path/to/sample_key.tsv`) used to map internal sample IDs to de-identified LRS IDs shown on plots:

| Column | Description |
|--------|-------------|
| `Sample` | Internal de-identified ID |
| `LRS_ID` | Long-read ID shown on plots |

---

## Figure-specific inputs and notes

### fig2.R — FSHD D4Z4 copy number, methylation, and limit of detection

Reads from the d4z4ling results directory. Key inputs:

- `/path/to/d4z4ling/results_directory` — d4z4ling output per sample
- `/path/to/subsampling_experiment/input_output` — subsampling experiment output directory (for fig2d limit-of-detection analysis)
- `/path/to/sample_key.tsv` — sample ID mapping

The subsampling experiment draws reads at varying multiplex factors and estimates the probability of detecting complete 4qA reads at different thresholds. See Methods for details.

---

### fig3.R — FSHD case studies

Reads from d4z4ling output. Also sources `fig2.R` for shared plotting functions, so `fig2.R` must be in the same working directory.

For the FSHD2 in-cis duplication case (fig3e), additional inputs are required from a phased hifiasm assembly re-annotated with d4z4ling:

- `/path/to/phased/assembly/results` — with the structure:
```
/{sample_id}/{sample_id}_{hap}/
    {sample_id}_realigned_modfreq_cov_summary.tsv
    {sample_id}_all_features.bed
```

Set `fshd2_incis_duplication_sample` and `asm_phase` variables at the top of the script accordingly.

---

### fig4/ — STR genotyping and promoter methylation

#### fig4b — OPDM STR genotyping (LongTR)

STR genotyping was performed with LongTR v1.0 on a pre-defined set of OPDM-associated loci (*ABCD3*, *GIPC1*, *NOTCH2NLC*, *RILPL1*). Run:

```bash
./fig4/fig4b_genotype_STR.sh
```

See Methods for LongTR parameters and the STR locus definitions (Table S1).

#### fig4c — Promoter methylation (bedtools)

```bash
./fig4/fig4c_promoter_methylation.sh
```

The resulting TSV can be visualised using Prism, Excel, or a similar tool.

---

### fig5_mtDNA_analysis/ — mtDNA deletion analysis (fig5i)

#### SRA data download

```bash
prefetch ${SRA_ID}
./fig5_mtDNA_analysis/map_sra_data_to_ref.sh ${SRA_ID}
```

Dependencies: `samtools`, `minimap2`, `sratoolkit`.

#### Run analysis

Create the output table header:

```bash
output_table="/path/to/analysis_dir/sample_proportion_summary_1kbp_deletion.tsv"
echo -e "Sample\tTotal_reads\tTotal_reads_over_1kbp\tReads_with_big_deletion\tchrM_coverage_mean\tmyop_panel_cov_mean\tratio_mtcov_to_ontargetcov\tProportion_with_big_deletion\tProportion_with_big_deletion_over1kbp" > "$output_table"
```

Then per sample:

```bash
sample=your_sample_name
HAPLOTAGGED_BAM=/path/to/sample_bam_file
root_dir=/path/to/analysis_dir
REGION_BED=/path/to/chrM_region_hg38.bed
output_table="/path/to/analysis_dir/sample_proportion_summary_1kbp_deletion.tsv"

./fig5_mtDNA_analysis/extract_deletion_cigar.sh "$sample" "$HAPLOTAGGED_BAM" "$root_dir" "$REGION_BED" "$output_table"
```

Dependencies: `python3` (pandas), `samtools`, `mosdepth`.

The resulting `$output_table` is used by `fig5i.R` to compare the proportion of mtDNA deletion-containing reads across samples.

---

### figS1.R — Cohort demographics

Requires `/path/to/demographics_matrix.tsv` with the following columns:

- `ID`, `Sex assigned at birth`, `Age at onset (years)`, `Age at time of recruitment (years)`, `Disease duration at time of recruitment (years)`
- Previous testing columns (Y/blank): `Exome-based neuromuscular panel`, `Genome sequencing`, `FSHD SB`, `DM1 test`, `DM2 test`, `DMD MLPA`, `Muscle biopsy`

---

### figS2.R — Sequencing coverage summary

Requires:

- `/path/to/sample_coverage_summary.tsv` — columns: `sample`, `total_bp_(Mbp)`, `on_target_bases_(Mbp)`, `n50`, `n50_on_target`
- `/path/to/myopathy_panel.bed` — columns: `chr`, `start`, `end`, `gene`
- `/path/to/coverage_per_gene_per_sample.tsv` — columns: `Sample`, `Gene`, `Coverage`
- `/path/to/mt_dna_coverage.txt` — columns: `Sample`, `Coverage`

---

### figS6.R — D4Z4 methylation and copy number vs age of onset

Reads from d4z4ling output. Also requires `/path/to/age_onset_group.tsv` with columns: `LRS_ID`, `Sample`, `Group` (FSHD1/FSHD2/FSHD1-biallelic/FSHD1+2), `Age at onset`.

Pearson correlation is computed between distal D4Z4 methylation, copy number, and age of onset. FSHD2 samples are excluded from the methylation vs age correlation by default (see script comments).

---

### figS7.R — Subsampling experiment (limit of detection)

Reads subsampling output from `/path/to/subsampling_experiment/output`. Expected file structure per sample and SD setting:

```
{sample_id}_subasampling_{sd_setting}/
    {sample_id}_read_ids_and_lengths_draw_summary.tsv
```

Three SD settings are supported: `5sd`, `10sd`, `10sd_50draws`. Plots show library size distribution (figS7a) and number of complete 4qA reads recovered vs library size (figS7b).

---

### figS8.R — OPDM STR allele size distributions

This figure compares allele sizes (in repeat units, log10 scale) between the study cohort (MyopCohort) and the 1000 Genomes long-read cohort (ont1000G) for four OPDM-associated STR loci.

**MyopCohort:** LongTR v1.0 was run per sample on each OPDM locus. ALT sequences were extracted from the resulting VCFs and repeat copy number was estimated using Tandem Repeat Finder (TRF; parameters: `2 7 7 80 10 50 500`, `-ngs` mode). Only TRF hits matching the expected repeat motif for each locus were retained. See `run_trf_from_alt_tsv.v2.sh` for the processing script.

**ont1000G:** Allele sizes were derived from a pre-existing TRGT-genotyped VCF (`1KGP.500.adotto.trgt.vcf.gz`).

Two input TSVs are required:

`/path/to/all_STR_loci.all_cohorts.tsv`:

| Column | Description |
|--------|-------------|
| `GENE` | Locus name |
| `SAMPLE_ID` | Sample identifier |
| `ALLELE_ID` | Allele index |
| `COHORT` | MyopCohort or ont1000G |

`/path/to/all_STR_loci.all_cohorts.copy_numbers.tsv`:

| Column | Description |
|--------|-------------|
| `GENE` | Locus name |
| `ALLELE_ID` | Allele index |
| `COPIES` | Repeat copy number (from TRF or TRGT) |
| `COPIES NORMALISED` | Normalised copy number |
| `MOTIF` | Repeat motif matched |
| `MOTIF_CANONICAL` | Canonical motif form |
| `COHORT` | MyopCohort or ont1000G |

Pathogenic expansion thresholds (in repeat units) are hardcoded in the script and shown as pink shading:

| Locus | Lower threshold | Upper threshold |
|-------|----------------|----------------|
| ABCD3 | 118 | 694 |
| NOTCH2NLC | 66 | 517 |
| RILPL1 | 120 | 197 |
| GIPC1 | 73 | 164 |

---

### sup_table_3.R — Supplementary Table 3 (SV annotations)

This script aggregates panel-overlapping SV calls across all samples, annotates them for rarity against gnomAD, the 1KGP long-read cohort (ONT1000G), and the internal cohort, and produces a per-sample summary TSV. SV prioritisation criteria and matching logic (reciprocal overlap thresholds, insertion distance criteria) are described in full in the Methods section of the paper.

Rarity thresholds (adjustable at the top of the script):

| Source | Threshold |
|--------|-----------|
| gnomAD | AF < 0.01 |
| ONT1000G | ≤ 2 alt alleles (~1% for n=149) |
| Internal cohort | ≤ 1 alt allele (~1% for n=61) |

Absence from a database is treated as rare. The output (`myopathy_svs_sample_summary.tsv`) reports per-sample counts of panel SVs, gene-overlapping SVs, and rare SVs per source, including a combined `n_rare_all` column for variants rare across all three sources. This table underlies Table S3 in the paper.

Three input files are required:

**`/path/to/all_svlog_database.txt`** — plain text list of paths, one per line, to per-sample `database_annotation.tsv` files from the svlog database.

**`/path/to/all_sv_annotations.txt`** — plain text list of paths to per-sample `*.hg38.svs.annotated.tsv` files containing VEP consequence and gene annotations.

**`/path/to/sniffles_overlapping_panel.tsv`** — panel-overlapping Sniffles SV calls with columns: `sample_id`, `CHROM`, `START`, `END`, `variant_id`, `panel`.

---
