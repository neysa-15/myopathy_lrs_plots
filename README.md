# Myopathy ONT figures
This repository contains the code to generate figures in [Targeted long-read sequencing enables comprehensive analysis of the genetic and epigenetic landscape of inherited myopathies](https://www.medrxiv.org/content/10.64898/2025.12.06.25340828v1)

## Project abstract
Inherited myopathies are a group of disorders with diverse and complex genetic aetiologies. The causative genetic variants vary widely in type, size and sequence context, encompassing small sequence variants and large structural variants in protein-coding genes, mitochondrial variants, repeat expansions, and more complex events, such as the chromosome 4 D4Z4 macrosatellite contraction and hypomethylation that causes facioscapulohumeral muscular dystrophy (FSHD). This poses a challenge for analysis with next generation sequencing and other molecular methods. This is further compounded by phenotypic variability between patients and phenotypic overlap of different forms of genetic myopathy. To accelerate myopathy research and improve diagnosis we have developed a targeted long-read sequencing assay and an integrated bioinformatics analysis framework that captures the full suite of genes, variants and epigenetic signatures currently implicated in all forms of inherited myopathy. Applying this to a cohort of 53 myopathy patients with and without previous genetic diagnoses, we demonstrate the analytical validity of our approach, as well as its improved accuracy and resolution compared to existing methods – especially for FSHD. Our LRS assay identified an array of new information about the genetic and epigenetic landscape of inherited myopathies and provided new diagnoses for 29% of patients who had remained undiagnosed following clinical genetic testing. In our cohort, FSHD and oculopharyngodistal myopathies were the most common new diagnoses that were missed or mis-diagnosed by standard clinical genetic testing. Our new method constitutes a single streamlined assay for comprehensive genetic and epigenetic characterisation of inherited myopathies.

 ## Inputs

 ### Path to R library
 `/path/to/your/r_library` containing these libraries:
 * `data.table`
 * `ggplot2`
 * `tidyr`
 * `patchwork`
 * `scales`
 * `dplyr`

 ### Sample key for most figures
 `/path/to/sample_key.tsv` containing two different de-identified name:
 * `Sample`: de-identified ID
 * `LRS_ID`: long read id (this is what we show on the plots)

### fig S1
`/path/to/demographics_matrix.tsv` containing:
* `ID`: Sample's de-identified ID
* `Sex assigned at birth`: Male/female
* `Age at onset (years)`
* `Age at time of recruiment (years)`
* `Disease duration at time of recruitment (years)`
* `Previous testing`: Any list of previous testing which was done, put `Y` if test is done and `` if not, test done in our cases
   * `Exome-based neuromuscular panel`
   * `Genome sequencing`
   * `FSHD SB`
   * `DM1 test`
   * `DM2 test`
   * `DMD MLPA`
   * `Muscle biopsy`

 ### fig S2
`/path/to/sample_coverage_summary.tsv` need to have at least:
* `sample`: sample name
* `total_bp_(Mbp)`: Total bases obtained from sequencing
* `on_target_bases_(Mbp)`: Total bases obtained from sequencing for the target panel
* `n50`: Sequencing N50 
* `n50_on_target`: Sequencing N50 for target panel 

`/path/to/myopathy_panel.bed` need to have at least 4 columns:
* `chr`: chromosome
* `start`: chromosome start coordinates
* `end`: chromosome end coordinates
* `gene`: Gene of interest

`/path/to/coverage_per_gene_per_sample.tsv` containing:
* `Sample`: sample name
* `Gene`: Gene of interest
* `Coverage`: which is the coverage calculated per gene per sample

`/path/to/mt_dna_coverage.tx` containing:
* `Sample`: sample name
* `Coverage`: mtDNA coverage of sample
