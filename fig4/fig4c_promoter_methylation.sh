#!/bin/bash

## bedtools versions used v2.31.0

## name of sample input files
PREFIX=sample1
METH=${PREFIX}.hg38.minimod.modfreqs_combined.bed

## list of promoter sites tested here are as follows. These are stored in BED files name ${GENE}_promoter_hg38.bed
# chr1:94417576-94419501_abcd3_promoter
# chr19:14495136-14496685_gipc1_promoter
# chr1:149389995-149392084_notch2nlc_promoter
# chr12:123531630-123536213_rilpl1_promoter

## extract relevant fields for methylation analysis from minimod modfreqs BED file 
cat ${METH} | awk '$9 == "*" && $5 >= 5' > ${PREFIX}_tmp.bed

## intersect minimod modfreqs BED file with each promoter region to obtain methylation frequency of each CpG site within each promoter region in a convenient TSV format
for GENE in 'abcd3' 'gipc1' 'notch2nlc' 'rilpl1'
do
SITE=${GENE}_promoter_hg38.bed
SITE_FREQS=${PREFIX}.${GENE}_promoter_modfreqs_combined.tsv
bedtools intersect -a ${SITE} -b ${PREFIX}_tmp.bed -wb | cut -f 11 > ${SITE_FREQS}
done