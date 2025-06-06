import argparse
import os
from typing import Tuple
import pandas as pd
from Bio import SeqIO
from utils.utils import cut_amplicon


def read_fasta_file(fasta_file: str) -> list[Tuple[str, str]]:
    """
    Function to read a fasta file and return a list of tuples with the sequence id and sequence
    
    Args:
        fasta_file: str, path to the fasta file
    
    Returns:
        list of tuples with the sequence id and sequence
    """
    seqs = []
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            seqs.append((record.id, record.seq))
    return seqs


def process_cliquesnv_output(out_file_dir: str, ref_seq:str):
    """
    Function to process the output of cliquesnv and standardize it

    Args:
        out_file_dir: str, path to the output directory
        out_file_dir: str, reference sequence

    Returns:
        None
    """
    # locate .fasta file
    sample = out_file_dir.split("/")[-2].strip()
    region = out_file_dir.split("/")[-1].strip()
 
    fasta_file_path = os.path.join(out_file_dir,  "region_" + region + ".fasta")
    
    # if file does not exist, return
    if not os.path.exists(fasta_file_path):
        print("File does not exist", fasta_file_path)
        return

    # read fasta file
    seqs = read_fasta_file(fasta_file_path)
    # figure out if it is per region or whole genome

    region = out_file_dir.split("/")[-1]


    # make output_file within the same directory named standard_output.tsv
    output_file = os.path.join(out_file_dir, "standard_output.tsv")

    with open(output_file, "w") as f:
        f.write("haplotype_id\tregion\trel_abundance\tsequence\n")
        for sample_id, seq in seqs:
            # >haplotype_id_fr_rel_abundance
            haplotype_id = sample_id.split("_")[0].strip()
            rel_abundance = float(sample_id.split("_")[2])
            seq = seq.upper().strip()
            start_region = int(region.split("_")[0])
            end_region = int(region.split("_")[1])
            try:
                reconstructed_amplicon = cut_amplicon(seq, ref_seq, start_region, end_region)
            except: 
                print("Could not reconstruct the amplicon")
                continue
            print("New haplotype with id: ", haplotype_id)
            f.write("{}\t{}\t{}\t{}\n".format(haplotype_id, region, rel_abundance, reconstructed_amplicon))

        f.close()

    return

def process_haplodmf_output(out_file_dir: str):
    """
    Function to process the output of haplodmf and standardize it

    Args:
        out_file_dir: str, path to the output directory

    Returns:
        None
    """

    fasta_file = os.path.join(out_file_dir, "haplodmf_haplotypes.fasta")

    # if file does not exist, return
    if not os.path.exists(fasta_file):
        print("File does not exist", fasta_file)
        return

    seqs = read_fasta_file(fasta_file)
    region = out_file_dir.split("/")[-1]
    
    # make output_file within the same directory named standard_output.tsv
    output_file = os.path.join(out_file_dir, "standard_output.tsv")
    
    f = open(output_file, "w")
    f.write("haplotype_id\tregion\trel_abundance\tsequence\n")

    for sample_id, seq in seqs:
        # >haplotype_0_length_1072_abundance_1_number_of_reads_29_depth_28.86007462686567
        haplotype_id = sample_id.split("_")[1]
        rel_abundance = float(sample_id.split("_")[5])
        seq = seq.upper().strip()
        f.write("{}\t{}\t{}\t{}\n".format(haplotype_id, region, rel_abundance, seq))
    
    f.close()

    return

def process_rvhaplo_output(out_file_dir: str):
    """

    Args:
        out_file_dir: str, path to the output directory
    
    Returns:
        None
    """
    
    fasta_file = os.path.join(out_file_dir, "rvhaplo_haplotypes.fasta")

    # if file does not exist, return
    if not os.path.exists(fasta_file):
        print("File does not exist", fasta_file)
        return

    seqs = read_fasta_file(fasta_file)

    region = out_file_dir.split("/")[-1]

    output_file = os.path.join(out_file_dir, "standard_output.tsv")
    f = open(output_file, "w")
    f.write("haplotype_id\tregion\trel_abundance\tsequence\n")

    for sample_id, seq in seqs:
        # >haplotype_0_length_25801_abundance_1_number_of_reads_4040_depth_169.81744586900103
        haplotype_id = sample_id.split("_")[1]
        rel_abundance = float(sample_id.split("_")[5])
        seq = seq.upper().strip()
        f.write("{}\t{}\t{}\t{}\n".format(haplotype_id, region, rel_abundance, seq))
    
    f.close()
    
    return 

def main():
    parser = argparse.ArgumentParser(description='Read tools output and standardize it')
    parser.add_argument('--results_dir', dest='results', type=str, required=True, help='Results directory')
    parser.add_argument('--hrt', dest='hrt', type=str, required=True, help='Haplotype reconstruction tool')
    parser.add_argument('--ec', dest='ec', type=str, required=True, help='Error correction tool')
    parser.add_argument('--sample', dest='sample', type=str, required=True, help='LUMC sample name')
    parser.add_argument('--region', dest='region', type=str, required=True, help='Genomic region')
    parser.add_argument('--ref_seq', dest='ref_seq', type=str, required=True, help='Reference sequence')
    parser.add_argument('--virus', dest='virus', type=str, required=True, help='Virus')
    parser.add_argument('--coverage', dest='coverage', type=int, required=True, help='Coverage')

    args = parser.parse_args()

    # folder structure is as follows:
    # results_dir
    #   - tool_1
    #      - virus_name
    #         - coverage
    #           - regions
    #               - ec_tool
    #                   - sample
    #                       - region
    #                           - output_file


    results_dir = args.results
    hrt = args.hrt
    ec = args.ec
    region = args.region
    sample = args.sample
    ref_seq = args.ref_seq
    ref_seq = str(SeqIO.read(ref_seq, "fasta").seq.upper())
    virus = args.virus
    coverage = str(args.coverage)
    # get the output directory for the tool
    outptut_dir = os.path.join(results_dir, hrt, virus, coverage, "regions", ec, sample, region)

    if hrt == "cliquesnv":
        process_cliquesnv_output(outptut_dir, ref_seq)
    elif hrt == "haplodmf":
        process_haplodmf_output(outptut_dir)        
    elif hrt == "rvhaplo":
        process_rvhaplo_output(outptut_dir)
    else:
        raise ValueError("Invalid haplotype reconstruction tool")



if __name__ == '__main__':
    main()