import pandas as pd
import os
import re
import sys

DEL_SIZE = 1000

sample = sys.argv[1]
cigar_file = sys.argv[2]
cov_depth = sys.argv[3]
mosdepth_summary_file = sys.argv[4]

# -------------------------------------------------------------------
# COUNT coverage
def get_region_cov(filename, pattern="total_region"):
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            for line in f:
                if re.search(pattern, line):
                    summary = line.strip().split('\t')
                    return {
                        "region": summary[0],
                        "length": int(summary[1]),
                        "bases": int(summary[2]),
                        "mean": float(summary[3]),
                        "min": float(summary[4]),
                        "max": float(summary[5])
                    }

    # If pattern not found, return None or zeros
    return {
        "region": pattern,
        "length": 0,
        "bases": 0,
        "mean": 0,
        "min": 0,
        "max": 0
    }

ontarget_region_cov = get_region_cov(mosdepth_summary_file)
ontarget_mean_cov = ontarget_region_cov["mean"]

ratio = float(cov_depth) / ontarget_mean_cov if ontarget_mean_cov > 0 else 0
# print(f"Ratio of mt coverage to on-target coverage: {ratio:.2f}")
# -------------------------------------------------------------------

count_big_del = 0
reads_over_1kbp = 0
with open(cigar_file, "r+") as file:
    content = file.readlines()
    for line in content:
        #print(line.strip('\n'))

        matches = re.findall(r'(\d+)([MIDNSHP=X])', line)
        length = sum(int(length) for length, cigar_symbol in matches)

        if length >= DEL_SIZE:
        
            deletion_list = re.findall(r'(\d+)D', line)
            
            big_deletion = [x for x in deletion_list if int(x) > DEL_SIZE]
            if big_deletion:
                count_big_del += 1

            reads_over_1kbp += 1
            
num_reads = len(content)
proportion = count_big_del / num_reads if num_reads > 0 else 0
proportion_denom_over1kbp = count_big_del / reads_over_1kbp if reads_over_1kbp > 0 else 0
print(f"{sample}\t{num_reads}\t{reads_over_1kbp}\t{count_big_del}\t{cov_depth}\t{ontarget_mean_cov}\t{ratio:.2f}\t{proportion * 100:.2f}%\t{proportion_denom_over1kbp * 100:.2f}%")
