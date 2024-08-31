
import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO
from utils.read_simulation import simulate_ONT_HQ_reads, simulate_hifi_reads

def main():
    parser = argparse.ArgumentParser(description="Simulates pacbio hifi or ONT amplicon reads")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="path to data directory")
    parser.add_argument('--cores', dest = 'cores', default=8, required=False, type=int, help="cpu cores for the simulation")
    parser.add_argument('--strategy', dest = 'strategy', default="pacbio-hifi", required=True, type=str, help="ONT or pacbio-hifi")    
    args = parser.parse_args()

    cores = args.cores
    strategy = args.strategy
    directory = args.dir

    # get all fasta files in directory
    files = os.listdir(directory)
    # take only fasta files
    files = [file for file in files if file.endswith(".template")]
    
    for file in files:
        identifier = file.split(".template")[0]
        if strategy == "pacbio-hifi":
            simulate_hifi_reads(directory, identifier, cores)
        if strategy == "ONT": 
            simulate_ONT_HQ_reads(directory, identifier, cores)


    print("Done")



if __name__ == "__main__":
    sys.exit(main())