#!/bin/bash
#PBS -P ox63
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l storage=gdata/kr68+gdata/if89+scratch/kr68+gdata/ox63+scratch/ox63
#PBS -l jobfs=10GB
#PBS -l wd

# export MODULEPATH=$MODULEPATH:/g/data/if89/apps/modulefiles/

# module unload samtools
# module unload python3
# module unload mosdepth

# module load samtools/1.22
# module load python3/3.12.1
# module load mosdepth/0.3.9

sample=$1
HAPLOTAGGED_BAM=$2
root_dir=$3
REGION_BED=$4
output_table=$5

OUTDIR="${root_dir}/result_cigar/${sample}"

# Extract the whole chrM
chrM_BAM="${OUTDIR}/chrM.bam"
# If no bam then run samtools
if [ ! -f "$chrM_BAM" ]; then
    echo "Extracting chrM reads to ${chrM_BAM}"
    samtools view -L "$REGION_BED" -Sb "$HAPLOTAGGED_BAM" |\
        samtools sort -o "$chrM_BAM"
    samtools index "$chrM_BAM"
fi

echo "Calculating coverage depth for chrM ${chrM_BAM}"
cov_depth=$(samtools coverage -r chrM "$chrM_BAM" | grep "chrM" | cut -f7)

# coverage to the rest of myop panel region
mosdepth_summary="${OUTDIR}/on_target_cov/depth.mosdepth.summary.txt"
if [ -f "$mosdepth_summary" ]; then
    echo "On-target coverage summary found: ${mosdepth_summary}"
else
    echo "Calculating on-target coverage for ${HAPLOTAGGED_BAM}"
    MYOP_PANEL_BED=/path/to/myopathy_panel.hg38.bed
    mkdir -p ${OUTDIR}/on_target_cov
    mosdepth ${OUTDIR}/on_target_cov/depth \
        $HAPLOTAGGED_BAM \
        -b $MYOP_PANEL_BED
    mosdepth_summary="${OUTDIR}/on_target_cov/depth.mosdepth.summary.txt"
fi

# grab cigar strings
echo "Get CIGAR strings for ${chrM_BAM}"
CIGAR="${OUTDIR}/cigarstrings.txt"
samtools view "$chrM_BAM" | cut -f6 > $CIGAR

# grep "[0-9]*D" $CIGAR
echo "Calculating proportion of big deletions for sample ${sample}"
bigdel_pyscript="${root_dir}/bigdel_proportion.py"
python3 "$bigdel_pyscript" "$sample" "$CIGAR" "$cov_depth" "$mosdepth_summary" >> "$output_table"

