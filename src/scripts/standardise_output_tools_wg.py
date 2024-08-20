import argparse
import os
from typing import Tuple
import pandas as pd
from Bio import SeqIO

def cut_amplicon(seq: str, ref_seq: str, start: int, end: int) -> str:
    """
    Function to cut the amplicon from the sequence based on the reference sequence

    Args:
        seq: str, sequence
        ref_seq: str, reference sequence
        start: int, start position of the amplicon
        end: int, end position of the amplicon

    Returns:
        str, the reconstructed amplicon sequence
    """ 
    # write the reference sequence and the sequence to a temporary file
    with open("temp.fasta", "w") as f:
        f.write(">ref\n{}\n".format(ref_seq))
        f.write(">seq\n{}\n".format(seq))
        f.close()
    
    # align the sequences using mafft
    os.system("mafft --quiet temp.fasta > temp_aligned.fasta")

    # read the aligned sequences
    aligned_seqs = read_fasta_file("temp_aligned.fasta")

    # get the amplicon sequence 
    # ensure that the reference sequence is the first sequence in the aligned sequences
    if aligned_seqs[0][0] == "seq":
        amplicon_seq = aligned_seqs[0][1][start:end]
    else:
        amplicon_seq = aligned_seqs[1][1][start:end]

    amplicon_seq = str(amplicon_seq).replace("-", "").strip().strip("N")
    
    # clean up
    os.system("rm temp.fasta temp_aligned.fasta")

    return amplicon_seq

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


def process_cliquesnv_output(out_file_dir: str, genomic_regions: list[Tuple[int, int]], ref_seq: str):
    """
    Function to process the output of cliquesnv and standardize it

    Args:
        out_file_dir: str, path to the output directory

    Returns:
        None
    """
    # locate .fasta file
    fasta_file = [f for f in os.listdir(out_file_dir) if f.endswith(".fasta")][0]
    fasta_file_path = os.path.join(out_file_dir, fasta_file)    
    # read fasta file
    seqs = read_fasta_file(fasta_file_path)

    # make output_file within the same directory named standard_output.tsv
    output_file = os.path.join(out_file_dir, "standard_output.tsv")

    with open(output_file, "w") as f:
        f.write("haplotype_id\tregion\trel_abundance\tsequence\n")
        for sample_id, seq in seqs:
            # >haplotype_id_fr_rel_abundance
            haplotype_id = sample_id.split("_")[0][1:] 
            rel_abundance = float(sample_id.split("_")[1])
            # this is a full length sequence
            seq = seq.upper().strip()

            for region in genomic_regions:
                start, end = region
                region_str = "{}_{}".format(start, end)
                # cut the amplicon
                reconstructed_amplicon = cut_amplicon(seq, ref_seq, start, end)
                f.write("{}\t{}\t{}\t{}\n".format(haplotype_id, region_str, rel_abundance, reconstructed_amplicon))

        f.close()
    return

def process_haplodmf_output(out_file_dir: str, genomic_regions: list[Tuple[int, int]], ref_seq: str):
    """
    Function to process the output of haplodmf and standardize it

    Args:
        out_file_dir: str, path to the output directory

    Returns:
        None
    """

    fasta_file = os.path.join(out_file_dir, "haplodmf_haplotypes.fasta")
    seqs = read_fasta_file(fasta_file)

    # make output_file within the same directory named standard_output.tsv
    output_file = os.path.join(out_file_dir, "standard_output.tsv")
    f = open(output_file, "w")
    f.write("haplotype_id\tregion\trel_abundance\tsequence\n")

    for sample_id, seq in seqs:
        # >haplotype_0_length_1072_abundance_1_number_of_reads_29_depth_28.86007462686567
        haplotype_id = sample_id.split("_")[1]
        rel_abundance = float(sample_id.split("_")[5])
        # this is a full length sequence
        seq = seq.upper().strip()

        for region in genomic_regions:
            start, end = region
            region_str = "{}_{}".format(start, end)
            # cut the amplicon
            reconstructed_amplicon = cut_amplicon(seq, ref_seq, start, end)
            f.write("{}\t{}\t{}\t{}\n".format(haplotype_id, region_str, rel_abundance, reconstructed_amplicon))

    f.close()

    return

def process_rvhaplo_output(out_file_dir: str, genomic_regions: list[Tuple[int, int]], ref_seq: str):
    """

    Args:
        out_file_dir: str, path to the output directory
    
    Returns:
        None
    """
    
    fasta_file = os.path.join(out_file_dir, "rvhaplo_haplotypes.fasta")

    seqs = read_fasta_file(fasta_file)

    output_file = os.path.join(out_file_dir, "standard_output.tsv")
    f = open(output_file, "w")
    f.write("haplotype_id\tregion\trel_abundance\tsequence\n")

    for sample_id, seq in seqs:
        # >haplotype_0_length_25801_abundance_1_number_of_reads_4040_depth_169.81744586900103
        haplotype_id = sample_id.split("_")[1]
        rel_abundance = float(sample_id.split("_")[5])
        # this is a full length sequence
        seq = seq.upper().strip()

        for region in genomic_regions:
            start, end = region
            region_str = "{}_{}".format(start, end)
            # cut the amplicon
            reconstructed_amplicon = cut_amplicon(seq, ref_seq, start, end)
            f.write("{}\t{}\t{}\t{}\n".format(haplotype_id, region_str, rel_abundance, reconstructed_amplicon))

    f.close()
    
    return 

def main():
    parser = argparse.ArgumentParser(description='Read tools output and standardize it')
    parser.add_argument('--results_dir', dest='results', type=str, required=True, help='Results directory')
    parser.add_argument('--hrt', dest='hrt', type=str, required=True, help='Haplotype reconstruction tool')
    parser.add_argument('--ec', dest='ec', type=str, required=True, help='Error correction tool')
    parser.add_argument('--sample', dest='sample', type=str, required=True, help='LUMC sample name')
    parser.add_argument('--ref_seq', dest='ref_seq', type=str, required=True, help='Reference sequence')
    args = parser.parse_args()

    # folder structure is as follows:
    # results_dir
    #   - tool_1
    #       - ec_tool_1
    #           - per_region
    #               - sample_1
    #                   - region_1
    #                       - output_file
    #           - whole_genome
    #               - sample_1
    #                   - output_file

    results_dir = args.results
    hrt = args.hrt
    ec = args.ec
    sample = args.sample
    ref_seq = args.ref_seq

    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                    (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                    (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                    (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]


    # get the output directory for the tool
    outptut_dir = os.path.join(results_dir, hrt, ec, "whole_genome", sample)

    if hrt == "cliquesnv":
        process_cliquesnv_output(outptut_dir, genomic_regions, ref_seq)
    elif hrt == "haplodmf":
        process_haplodmf_output(outptut_dir, genomic_regions, ref_seq)    
    elif hrt == "rvhaplo":
        process_rvhaplo_output(outptut_dir, genomic_regions, ref_seq)
    else:
        raise ValueError("Invalid haplotype reconstruction tool")


if __name__ == '__main__':
    main()