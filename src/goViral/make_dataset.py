
import pandas as pd
import os
import sys
import argparse
from typing import List, Tuple
from utils.utils import get_genomic_regions, map_to_correct_region
import numpy as np
from multiprocessing import Process
from Bio import SeqIO

def get_bams_and_beds(ref_seq_file_path:str, input_fastq_file_path:str) -> str:
    """
    Maps the sequencing reads to their corresponding genomic region and strand

    Args:
    ref_seq_file_path: Path to the reference sequence file
    input_fastq_file_path: Path to the input file containing the fastq reads

    Returns:
    bed_file: Path to the bed file
    """
 
    output_file = input_fastq_file_path[:-6]

    os.system("minimap2 -ax map-pb {} {} > {}.sam".format(ref_seq_file_path, input_fastq_file_path, output_file))
    os.system("samtools view -bS {}.sam > {}.bam".format(output_file, output_file))
    os.system("samtools sort {}.bam -o {}.sorted.bam".format(output_file, output_file))
    os.system("samtools index {}.sorted.bam".format(output_file))
    os.system("bedtools bamtobed -i {}.sorted.bam > {}.bed".format(output_file, output_file))

    bed_file = output_file + ".bed"

    return bed_file

def main():
    parser = argparse.ArgumentParser(description="Maps sequencing reads to their corresponding genomic region and strand")
    parser.add_argument('--fastq', dest = 'fastq', required=True, type=str, help="Input file containing the fastq reads")
    parser.add_argument('--primers', dest = 'primers', required=True, type=str, help="File containing primer positions. To be used for mapping the amplicon reads to the correct genomic regions")
    parser.add_argument('--ref_seq', dest = 'ref_seq', required=True, type=str, help="Reference sequence. To be used for mapping the amplicon reads to the correct genomic regions")
    parser.add_argument('--out', dest = 'out', required=True, type=str, help="Output file path (must be a TSV file; ends in .tsv)")
    args = parser.parse_args()

    fastq_file = args.fastq
    ref_seq = args.ref_seq
    outfile = args.out
    primers = args.primers

    # if outfile exists, delete it
    if os.path.exists(outfile):
        os.remove(outfile)

    outdir = "/".join(outfile.split("/")[:-1])
    os.system("mkdir -p {}".format(outdir))

    # create it again
    open(outfile, 'a').close()

    bed_file = get_bams_and_beds(ref_seq, fastq_file)

    genomic_regions = get_genomic_regions(primers)

    # Read in the bed file
    bed = pd.read_csv(bed_file, sep="\t", header=None)
    bed.columns = ["ref_id", "start", "end", "id", "score", "strand"]

    reads = {}

    for record in SeqIO.parse(fastq_file, "fastq"):
        reads[record.id] = str(record.seq)

    for read_name in reads.keys():
        row = bed[bed['id'] == read_name]
        # check if row is empty
        if row.empty:
            continue

        start = int(row['start'].values[0])
        end = int(row['end'].values[0])
        strand = row['strand'].values[0]
            
        # map to correct region
        start, end = map_to_correct_region(start, genomic_regions)

        # write to tsv file
        with open(outfile, "a") as f:
            f.write(read_name + "\t" + reads[read_name] + "\t" + str(start) + "\t" + str(end) + "\t" + strand +"\n")
            f.close()
    
               
if __name__ == "__main__":
    sys.exit(main())