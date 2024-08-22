import os 
from typing import List, Tuple
from Bio import SeqIO
import pandas as pd


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