import pandas as pd
import os
import sys
import argparse
from utils.file_parsing import parse_fastq_ONT

def parse_maf_file_ONT(maf_file):
    maf_dict = {}
    with open(maf_file, "r") as f:
        lines = f.readlines()

    cnt = 0
    while cnt < len(lines):
        if lines[cnt].startswith("a"):
            next_line = lines[cnt+1].strip().split(" ")
            genomic_region = next_line[1].split(":")[1]
            strand = next_line[1][0]
            next_line = lines[cnt+2].strip().split(" ")
            sim_read_id = next_line[1].strip()
            if genomic_region not in maf_dict.keys():
                maf_dict[genomic_region] = set()
                maf_dict[genomic_region].add((sim_read_id, strand))
            else:
                maf_dict[genomic_region].add((sim_read_id, strand))

            cnt += 4
    
    return maf_dict


def find_identifier(file):
    # remove the .fasta
    identifier = file[:-6]
    return identifier

def make_output_file(output):
    if os.path.exists(output):
        os.remove(output)
    
    if  os.path.exists("/".join(output.split("/")[:-1])) != True:
        os.mkdir("/".join(output.split("/")[:-1]))
    # create new file
    open(output, "x").close()  

    with open(output, "w") as f:
        f.write("read_id\tgenomic_region\tlabel\n")
        f.close()
    
    return

def parse_maf_file(maf_file):
    maf_dict = {}
    with open(maf_file, "r") as f:
        lines = f.readlines()

    cnt = 0
    while cnt < len(lines):
        if lines[cnt].startswith("a"):
            next_line = lines[cnt+1].strip().split(" ")
            genomic_region = next_line[1].split(":")[1]
            strand = next_line[1].split("_")[0]
            next_line = lines[cnt+2].strip().split(" ")
            sim_read_id = "/".join(next_line[1].split("/")[:-1])
            if genomic_region not in maf_dict.keys():
                maf_dict[genomic_region] = set()
                maf_dict[genomic_region].add((sim_read_id, strand))
            else:
                maf_dict[genomic_region].add((sim_read_id, strand))

            cnt += 4
            
    return maf_dict


def parse_fastq(fastq_file):

    reads = {}

    with open(fastq_file, "r") as f:
        lines = f.readlines()

    lines = [x.strip() for x in lines]

    l = 0
    while(l!= len(lines)):
        # @S/1/ccs -> lines[l].split(" ")[0]
        # S/1/ccs -> lines[l].split(" ")[0][1:]
        # [S, 1, css] lines[l].split(" ")[0][1:].split("/")
        # [S, 1] -> lines[l].split(" ")[0][1:].split("/")[:-1]
        # S/1 -> "/".join(lines[l].split(" ")[0][1:].split("/")[:-1])
        read_name = "/".join(lines[l].split(" ")[0][1:].split("/")[:-1])
        l+=1
        read = lines[l]
        l+=2
        # quality_score = lines[l]
        l+=1
        reads[read_name] = read

    return reads

def write_read_entry_to_file(read, genomic_region, label, file):
    # open file and write read to file
    with open(file, "a") as tsvfile:
        tsvfile.write(read + "\t" + genomic_region + "\t" + label + "\n")
    return 

def write_sequence_to_file(full_read_id, sequence, outdir_read_dir):
    # check if file exists
    if os.path.exists(outdir_read_dir + "/" + full_read_id + ".fasta"):
        os.remove(outdir_read_dir + "/" + full_read_id + ".fasta")

    if not os.path.exists(outdir_read_dir + "/" + full_read_id + ".fasta"):
        with open(outdir_read_dir + "/" + full_read_id + ".fasta", "w") as f:
            f.write(">" + full_read_id + "\n")
            f.write(sequence + "\n")
    return

def write_all_reads(reads, maf_info, label, file, outdir_read_dir, ont):
    for genomic_region in maf_info:
        read_ids_and_strands = maf_info[genomic_region]
        print(read_ids_and_strands, genomic_region, label)

        for (read_id, strand) in read_ids_and_strands:
            if ont == "True": 
                full_read_id = label + "_" + read_id
            else: 
                full_read_id = label + "_" + read_id.split("/")[0] + read_id.split("/")[1]
            if read_id in reads: 
                sequence = reads[read_id]
                if strand == "-":
                    sequence = sequence[::-1].translate(str.maketrans("ATGC", "TACG"))
                write_read_entry_to_file(full_read_id, genomic_region, label, file)
                write_sequence_to_file(full_read_id, sequence, outdir_read_dir)


def main():
    parser = argparse.ArgumentParser(description="List all reads in a directory")
    parser.add_argument('--data_dir', dest = 'data_dir', required=True, type=str, help="data directory to sample from")
    parser.add_argument('--out', dest = 'out_dir', required=True, type=str, help="output directory to write samples to")
    parser.add_argument('--ont', dest = 'ont', action="store_true", help="ONT reads or not")

    args = parser.parse_args()

    data_dir = args.data_dir
    out_dir = args.out_dir

    out_dir_tsv_file = out_dir + "/all_reads.tsv"
    make_output_file(out_dir_tsv_file)

    # make a /reads directory in the out_dir
    if not os.path.exists(out_dir + "/reads"):
        os.makedirs(out_dir + "/reads")

    # get all fastq files in data directory
    files = os.listdir(data_dir)
    # take only fastq files
    files = [file for file in files if file.endswith(".fastq")]

    for file in files:
        file_path = data_dir + "/" + file
        identifier = find_identifier(file)
        if args.ont:
            reads = parse_fastq_ONT(file_path)
        else:
            reads = parse_fastq(file_path)
        
        if args.ont:
            maf_info = parse_maf_file_ONT(data_dir + "/" + identifier + ".maf")
        else: 
            maf_info = parse_maf_file(data_dir + "/" + identifier + ".maf")

        write_all_reads(reads, maf_info, identifier, out_dir_tsv_file, out_dir + "/reads", args.ont)

if __name__ == "__main__":
    sys.exit(main())