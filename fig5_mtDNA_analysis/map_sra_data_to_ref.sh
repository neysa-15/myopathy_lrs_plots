#!/bin/bash
#PBS -P kr68
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l storage=gdata/kr68+gdata/if89+scratch/kr68
#PBS -l jobfs=10GB
#PBS -l wd

# module unload sratoolkit
# module unload minimap2
# module unload samtools

# module load sratoolkit/3.1.1
# module load minimap2/2.28
# module load samtools/1.22

SRA=$1 # e.g. SRR23886843

fastq-dump --stdout $SRA | gzip > muscle/${SRA}/${SRA}.fastq.gz

# Input
FASTQ=muscle/${SRA}/${SRA}.fastq.gz
OUTPUT_BAM=muscle/${SRA}/${SRA}_chrM_mapped.bam

# reference genome
REF=/path/to/hg38.reference.fa
CHR_M_BED=chrM_region_hg38.bed

# Map and filter chrM reads
minimap2 -ax map-ont $REF $FASTQ -t ${PBS_NCPUS:-8}| \
    samtools view -L ${CHR_M_BED} -Sb |\
    samtools sort -o $OUTPUT_BAM
samtools index $OUTPUT_BAM
