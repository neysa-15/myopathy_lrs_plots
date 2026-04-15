#!/bin/bash

## LongTR version used: v1.0

## genome reference file
GENOME=/path/to/hg38.reference.fa

## STR site for genotyping (in this example NOTCH2NLC)
REGION=/path/to/notch2nlc_str_hg38.bed

## name of sample input files
PREFIX=sample1
INBAM=${PREFIX}.hg38.minimap2.whatshap.sorted.haplotagged.bam # path to your mapped bam file

## name of sample output files
OUTVCF=${PREFIX}_longtr_notch2nlc.vcf.gz
OUTTSV=${PREFIX}_longtr_notch2nlc.tsv

## run LongTR to genotpye a given STR site for a given sample (in this case NOTCH2NLC genotyping for sample a sample called 'sample1')
LongTR --bams ${INBAM} \
	--fasta ${GENOME} \
	--regions ${REGION} \
	--tr-vcf ${OUTVCF} \
	--bam-samps ${PREFIX} --bam-libs ${PREFIX} \
	--min-reads 3 \
	--max-tr-len 10000 \
	--phased-bam

## extract the two allele sequences for genotpyed STR site from the LongTR VCF; take the reverse complement and print to a TSV file.
zcat ${OUTVCF} | grep -v "^#"| cut -f5 | sed 's/,/\t/g' | awk -v prefix=${PREFIX} '{print prefix"\t"$1"\t"$2}' > ${OUTTSV}

