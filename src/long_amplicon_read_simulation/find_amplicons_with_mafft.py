import argparse
import sys
import os
from Bio import SeqIO
import pandas as pd


def mafft_msa(fasta_seqs, directory): 
    os.system(f"mafft --auto --quiet --thread 4 {fasta_seqs} > {directory}/msa_alignment.fasta")
    return 

def find_amplicons(directory, start, end, record_id):

    msa_alignment = "{}/msa_alignment.fasta".format(directory)
    fasta_sequences = SeqIO.parse(open(msa_alignment),'fasta')
    
    for record in fasta_sequences:
        if record.id == record_id:
            if start >= len(record.seq) or end >= len(record.seq):
                raise ValueError("Start or end is out of range, skipping")
            amplicon = record.seq[start:end]
            return amplicon
    
    raise ValueError("Record not found")
      

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--input_fasta', dest = 'input_fasta', required=True, type=str, help="")
    parser.add_argument('--directory', dest = 'directory', required=True, type=str, help="")
    parser.add_argument('--start', dest = 'start', required=True, type=int, help="start of amplicon")
    parser.add_argument('--end', dest = 'end', required=True, type=int, help="end of amplicon")
    parser.add_argument('--record_id', dest = 'record_id', required=True, type=str, help="record id of sequence to be used for amplicon detection")
    args = parser.parse_args()

    input_fasta = args.input_fasta
    directory = args.directory
    start = args.start
    end = args.end
    record_id = args.record_id

    mafft_msa(input_fasta, directory)

    amplicon = find_amplicons(directory, start, end, record_id)
    # remove any gaps
    amplicon = amplicon.replace("-", "").strip().upper()

    print(amplicon)


if __name__ == "__main__":
    sys.exit(main())