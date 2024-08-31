

import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO
import subprocess
from utils.read_simulation import create_template

def main():
    parser = argparse.ArgumentParser(description="Creates amplicon fragments given  bed file from primalScheme")
    parser.add_argument('--input_fasta', dest = 'input_fasta', required=True, type=str, help="")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="path to data directory")
    parser.add_argument('--primers', dest = 'primers', required=True, type=str, help="path to bed file with primer locations")
    parser.add_argument('--cnt_templates', dest = 'cnt_templates',  required=False, default=10, type=int, help="number of templates to create per sequence")
    args = parser.parse_args()

    # get directories in data directory
    data_dir = args.dir
  
    files = os.listdir(data_dir)
    # take only fasta files
    files = [file for file in files if file.endswith(".fasta")]

    for file in files:

        file_path = os.path.join(data_dir, file)

        identifier = file.split(".fasta")[0]

        if identifier == "msa_alignment":
            continue
    
        seq_path = os.path.join(data_dir, identifier + ".fasta")

        template = create_template(args.input_fasta, seq_path, identifier, args.primers, args.cnt_templates)
        # write template to file
        template_file = os.path.join(data_dir, identifier + ".template")
        # if template file already exists, remove it
        if os.path.exists(template_file):
            os.remove(template_file)

        with open(template_file, "w") as out_file:
            out_file.write(template)

    
if __name__ == "__main__":
    sys.exit(main())