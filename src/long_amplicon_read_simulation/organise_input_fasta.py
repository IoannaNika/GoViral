
import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--input_fasta', dest = 'input_fasta', required=True, type=str, help="")
    args = parser.parse_args()

    # read sequences.fasta file
    fasta_sequences = SeqIO.parse(open(args.input_fasta),'fasta')

    directory = "/".join(args.input_fasta.split("/")[:-1])
    # create output directory in data directory
    out_dir = os.path.join(directory, "processed")

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # total number of sequences considered
    total_seqs = 0

    for fasta in fasta_sequences:

        strain, sequence = fasta.id, str(fasta.seq)
        identifier = strain.split(" ")[0].strip()
        
        total_seqs += 1

        # create output file name
        output_file = os.path.join(out_dir, identifier  + ".fasta")
        # if file exists, delete and rewrite
        if os.path.exists(output_file):
            os.remove(output_file)
            
        # write to file
        with open(output_file, "w") as out_file:
            out_file.write(">" + strain + "\n" + sequence + "\n")

    print("Total number of sequences considered: " + str(total_seqs))

if __name__ == "__main__":
    sys.exit(main())





