import argparse
import os
from typing import Tuple, List
import pandas as pd
from Bio import SeqIO

def compare_sequences(seq1: str, seq2: str) -> Tuple[bool, str]:
    """
    Function to compare two sequences and return True if they are the same, False otherwise

    Args:
        seq1: str, sequence 1
        seq2: str, sequence 2

    Returns:
        bool, True if the sequences are the same, False otherwise
        str, the corrected sequence if the sequences are the same, None otherwise
    """

    # align the sequences with MAFFT
    # write the sequences to a temporary file
    with open("temp.fasta", "w") as f:
        f.write(">seq1\n{}\n".format(seq1))
        f.write(">seq2\n{}\n".format(seq2))
        f.close()
    
    # align the sequences using mafft
    os.system("mafft --quiet temp.fasta > temp_aligned.fasta")

    # read the aligned sequences
    aligned_seqs = read_fasta_file("temp_aligned.fasta")

    # get the aligned sequences
    seq1_aligned = str(aligned_seqs[0][1]).upper()
    seq2_aligned = str(aligned_seqs[1][1]).upper()

    # remove temporary files
    os.system("rm temp.fasta temp_aligned.fasta")

    # replace the surrounding - with N
    curr_char = 0
    while curr_char < len(seq1_aligned) and seq1_aligned[curr_char] == "-" :
        seq1_aligned = seq1_aligned[:curr_char] + "N" + seq1_aligned[curr_char+1:]
        curr_char += 1
    
    curr_char = len(seq1_aligned) - 1
    while  curr_char >= 0 and seq1_aligned[curr_char] == "-" :
        seq1_aligned = seq1_aligned[:curr_char] + "N" + seq1_aligned[curr_char+1:]
        curr_char -= 1
    
    curr_char = 0
    while curr_char < len(seq2_aligned) and seq2_aligned[curr_char] == "-":
        seq2_aligned = seq2_aligned[:curr_char] + "N" + seq2_aligned[curr_char+1:]
        curr_char += 1

    curr_char = len(seq2_aligned) - 1
    while  curr_char >= 0 and seq2_aligned[curr_char] == "-":
        seq2_aligned = seq2_aligned[:curr_char] + "N" + seq2_aligned[curr_char+1:]
        curr_char -= 1
    
    corrected_seq = ""

    # check if they only differ by Ns, if one is the substring of the other, or if they are the same
    for i in range(len(seq1_aligned)):
        if seq1_aligned[i] == seq2_aligned[i]:
            corrected_seq += seq1_aligned[i]
            continue
        if seq1_aligned[i] != seq2_aligned[i] and (seq1_aligned[i] != "N" and seq2_aligned[i] != "N"):
            return False, None
            
        elif seq1_aligned[i] == "N" and seq2_aligned[i] != "N":
            corrected_seq += seq2_aligned[i]

        elif seq1_aligned[i] != "N" and seq2_aligned[i] == "N":
            corrected_seq += seq1_aligned[i]
    
    return True, corrected_seq

def check_and_update_if_haplotype_exists(standard_output_path: str, region: str, sequence: str, rel_ab: float) -> bool:
    """
    Function to check if a haplotype with the same sequence exists in the standard output file
    If it does, update the sequence and the relative abundance and return True else return False

    Args:
        standard_output_path: str, path to the standard output file
        region: str, genomic region
        sequence: str, sequence
        rel_ab: float, relative abundance
    
    Returns:
        bool, True if the haplotype exists and was updated, False otherwise
    """
    try:
        haplotypes = pd.read_csv(standard_output_path, sep="\t", header=0)
    except: 
        # file is empty
        return False

    if haplotypes.shape[0] == 0:
        return False

    for idx, row in haplotypes.iterrows():

        same_bool, corrected_seq = compare_sequences(row["sequence"], sequence)
        
        if row["region"] == region and same_bool:
            # update the row with the corrected sequence
            haplotypes.at[idx, "sequence"] = corrected_seq
            # make sure rel_abudance is a float
            haplotypes.at[idx, "rel_abundance"] = float(haplotypes.at[idx, "rel_abundance"])
            haplotypes.at[idx, "rel_abundance"] += rel_ab
            haplotypes.to_csv(standard_output_path, sep="\t", index=False, header=True)
            return True
        
    return False
   

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

def read_fasta_file(fasta_file: str) -> List[Tuple[str, str]]:
    """
    Function to read a fasta file and return a list of tuples with the sequence id and sequence
    
    Args:
        fasta_file: str, path to the fasta file
    
    Returns:
        list of tuples with the sequence id and sequence
    """

    if not os.path.exists(fasta_file):
        print("File does not exist", fasta_file)
        return []

    seqs = []

    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            seqs.append((record.id, str(record.seq).upper()))
    
    return seqs


def process_cliquesnv_output(out_file_dir: str, genomic_regions: List[Tuple[int, int]], ref_seq: str):
    """
    Function to process the output of cliquesnv and standardize it

    Args:
        out_file_dir: str, path to the output directory

    Returns:
        None
    """

    sample = out_file_dir.split("/")[-1].strip()
    # locate .fasta file
    fasta_file_path = os.path.join(out_file_dir, sample + ".fasta")
    
    if not os.path.exists(fasta_file_path):
        print("File does not exist", fasta_file_path)
        return
    
    # read fasta file
    seqs = read_fasta_file(fasta_file_path)

    if len(seqs) == 0:
        print("No sequences found", fasta_file_path)
        return

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
                
                # check if for the same genomic region, a haplotype with the same sequence has already been written
                # if so, skip writing it
                # if not, write it
                if not check_and_update_if_haplotype_exists(output_file, region_str, reconstructed_amplicon, rel_abundance) and len(reconstructed_amplicon) > 600:
                    f.write("{}\t{}\t{}\t{}\n".format(haplotype_id, region_str, rel_abundance, reconstructed_amplicon))

        f.close()

    return

def process_haplodmf_output(out_file_dir: str, genomic_regions: List[Tuple[int, int]], ref_seq: str):
    """
    Function to process the output of haplodmf and standardize it

    Args:
        out_file_dir: str, path to the output directory

    Returns:
        None
    """

    fasta_file = os.path.join(out_file_dir, "haplodmf_haplotypes.fasta")

    if not os.path.exists(fasta_file):
        print("File does not exist", fasta_file)
        return
    
    seqs = read_fasta_file(fasta_file)

    if len(seqs) == 0:
        print("No sequences found", fasta_file_path)
        return

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

            if not check_and_update_if_haplotype_exists(output_file, region_str, reconstructed_amplicon, rel_abundance) and len(reconstructed_amplicon) > 600:
                f.write("{}\t{}\t{}\t{}\n".format(haplotype_id, region_str, rel_abundance, reconstructed_amplicon))

    f.close()

    return

def process_rvhaplo_output(out_file_dir: str, genomic_regions: List[Tuple[int, int]], ref_seq: str):
    """
    Args:
        out_file_dir: str, path to the output directory
    """
    
    fasta_file = os.path.join(out_file_dir, "rvhaplo_haplotypes.fasta")

    if not os.path.exists(fasta_file):
        print("File does not exist", fasta_file)
        return
        
    seqs = read_fasta_file(fasta_file)

    if len(seqs) == 0:
        print("No sequences found", fasta_file)
        return

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
            
            if not check_and_update_if_haplotype_exists(output_file, region_str, reconstructed_amplicon, rel_abundance) and len(reconstructed_amplicon) > 600:
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

    # read the genomic sequence of the reference
    ref_seq = str(SeqIO.read(ref_seq, "fasta").seq.upper())

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