

import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO
import json
from utils.read_simulation import create_template, simulate_hifi_reads, simulate_ONT_HQ_reads
from utils.utils import remove_directory_contents


def main():
    parser = argparse.ArgumentParser(description="Simulates pacbio hifi amplicon reads")
    parser.add_argument('--data_dir', dest = 'data_dir', required=True, type=str, help="path to data directory")
    parser.add_argument('--input_fasta', dest = 'input_fasta', required=True, type=str, help="")
    parser.add_argument('--strategy', dest = 'strategy',required=True, type=str, help="ONT or pacbio-hifi")   
    parser.add_argument('--ids_and_ab', type=str, dest = 'ids_and_ab', required=True, help='A dictionary of ids and their relative abundances (%)')
    parser.add_argument('--coverage', dest = 'coverage', type=int, required=True, help='The coverage to simulate')
    parser.add_argument('--outdir', dest = 'outdir', type=str, required=True, help='The output directory')
    parser.add_argument('--primers', dest = 'primers', type=str, required=True, help ="")
    parser.add_argument('--cores', dest = 'cores', default=8, required=False, type=int, help="cpu cores for the simulation")

    args = parser.parse_args()

    data_dir = args.data_dir
    cores = args.cores
    strategy = args.strategy
    
    outdir = os.path.join(args.outdir, "output")
    remove_directory_contents(outdir)

    with open(args.ids_and_ab) as json_file:
        ids_and_abundances = json.load(json_file)
        
    coverage = args.coverage
    primers = args.primers

    # open file to write the ids and their abundances
    ids_abs_file = os.path.join(outdir, "ids_and_abundances.tsv")
    
    ids_abs_file = open(ids_abs_file, "w")
    
    for s_id in ids_and_abundances.keys():
        abundance = ids_and_abundances[s_id]
        ab_fraction = abundance
        n_reads = int(ab_fraction * coverage)
        ids_abs_file.write(s_id + "\t" + str(n_reads) + "\n")

        seq_path = os.path.join(data_dir, s_id + ".fasta")
        
        template_file = os.path.join(data_dir, s_id + ".template")

        final_template = create_template(args.input_fasta, seq_path, s_id, args.primers, n_reads)

        # if template file already exists, remove it
        if os.path.exists(template_file):
            os.remove(template_file)

        with open(template_file, "w") as out_file:
            out_file.write(final_template)
            out_file.close()

        if strategy == "pacbio-hifi":
            simulate_hifi_reads(data_dir, s_id, cores)
            
        if strategy == "ONT": 
            simulate_ONT_HQ_reads(data_dir, s_id, cores)

    ids_abs_file.close()

if __name__ == "__main__":
    sys.exit(main())

    

    

