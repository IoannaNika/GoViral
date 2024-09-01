import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches
from utils.utils import check_and_update_if_haplotype_exists

def adjust_relative_abundance(merged_file_path):
    
    df = pd.read_csv(merged_file_path, sep = "\t", header = 0)

    unique_grs = df["region"].unique()

    for gr in unique_grs:
        df_gr = df[df["region"] == gr]
        total_abundance = df_gr["rel_abundance"].sum()
        

        if total_abundance != 1:
            df.loc[df["region"] == gr, "rel_abundance"] = df_gr["rel_abundance"] / total_abundance
    
    df.to_csv(merged_file_path, sep='\t', index=False, header=True)

    return

def main():
    parser = argparse.ArgumentParser(description="Merge the standard output files from different seeds")
    parser.add_argument('--directory', dest = 'directory', required=True, type=str, help="Directory for input and output files")
    parser.add_argument('--gr_start', dest = 'gr_start', required=True, type=int, help="Genomic region start")
    parser.add_argument('--coverage', dest = 'coverage', default=100, required=False, type=int, help="Number of reads to consider in each subsample")
    parser.add_argument('--seed_limit', dest = 'seed_limit', default = 10, required=True, type=int, help="Number of subsamples to consider")
    parser.add_argument('--ab_threshold', dest = 'ab_threshold', required=True, type=float, help="Abundance threshold for filtering out low abundance sequences")
    args = parser.parse_args()

    seeds = range(1,args.seed_limit + 1)

    merged_file_path = os.path.join(args.directory, "merged_standard_output.tsv")

    if not os.path.exists(merged_file_path):

        with open(merged_file_path, 'w') as file:
            file.write('haplotype_id\tregion\trel_abundance\tsequence\n')  
            file.close()

    merged_results_tsv = pd.read_csv(merged_file_path,  sep='\t', header=0)

    for seed in seeds:
        consensus_file = os.path.join(args.directory, "seed_" + str(seed), "subsampled_reads_" + str(args.coverage) + "_" + str(args.gr_start),"standard_output.tsv")
        consensus_tsv = pd.read_csv(consensus_file, sep='\t', header=0)

        # write the sequences to the new tsv
        for index, row in consensus_tsv.iterrows():
            haplotype_id = row['haplotype_id']
            region = row["region"]
            sequence = row["sequence"]
            rel_ab = row['rel_abundance']
            
            if check_and_update_if_haplotype_exists(merged_file_path, region, sequence, rel_ab):
                continue
                
            with open(merged_file_path, 'a') as file:
                file.write(f"{haplotype_id}\t{region}\t{rel_ab}\t{sequence}\n")
                file.close()
    
    # adjust relative abundances to ensure that they add up to 1
    adjust_relative_abundance(merged_file_path)

    merged_results_tsv = pd.read_csv(merged_file_path, sep='\t', header = 0)

    filtered_df = merged_results_tsv[merged_results_tsv['rel_abundance'] >= args.ab_threshold]

    # Save the filtered DataFrame to a new TSV file
    filtered_df.to_csv(merged_file_path, sep='\t', index=False, header=True)

    adjust_relative_abundance(merged_file_path)

if __name__ == '__main__':
    main()
