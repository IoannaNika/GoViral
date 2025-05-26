import argparse
import os
from typing import Tuple, List
import pandas as pd
from Bio import SeqIO
from utils.utils import cut_amplicon, read_fasta_file, check_and_update_if_haplotype_exists, get_genomic_regions

def process_cliquesnv_output(out_file_dir: str, genomic_regions: List[Tuple[int, int]], ref_seq: str, ec_tool: str):
    """
    Function to process the output of cliquesnv and standardize it

    Args:
        out_file_dir: str, path to the output directory

    Returns:
        None
    """

    sample = out_file_dir.split("/")[-1].strip()
    # locate .fasta file
    fasta_file_path = os.path.join(out_file_dir, "reads_{}.fasta".format(ec_tool))
    
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
            haplotype_id = sample_id.split("_")[0].strip()
            rel_abundance = float(sample_id.split("_")[2])
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
                if not check_and_update_if_haplotype_exists(output_file, region_str, reconstructed_amplicon, rel_abundance):
                    print("Haplotype id: ", haplotype_id)
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
        print("No sequences found", fasta_file)
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

            if not check_and_update_if_haplotype_exists(output_file, region_str, reconstructed_amplicon, rel_abundance):
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
            
            if not check_and_update_if_haplotype_exists(output_file, region_str, reconstructed_amplicon, rel_abundance):
                f.write("{}\t{}\t{}\t{}\n".format(haplotype_id, region_str, rel_abundance, reconstructed_amplicon))

    f.close()
    
    return 

def main():
    parser = argparse.ArgumentParser(description='Read tools output and standardize it')
    parser.add_argument('--results_dir', dest='results', type=str, required=True, help='Results directory')
    parser.add_argument('--hrt', dest='hrt', type=str, required=True, help='Haplotype reconstruction tool')
    parser.add_argument('--ec', dest='ec', type=str, required=True, help='Error correction tool')
    parser.add_argument('--sample', dest='sample', type=str, required=True, help='Sample name')
    parser.add_argument('--ref_seq', dest='ref_seq', type=str, required=True, help='Reference sequence')
    parser.add_argument('--coverage', dest='coverage', type=str, required=True, help='Coverage')
    parser.add_argument('--virus', dest='virus', type=str, required=True, help='Virus name')
    parser.add_argument('--primers', dest='primers', type=str, required=False, default = "None", help='Primers file, if not given assumes sars-cov-2 amplicons (LUMC)')


    args = parser.parse_args()

    # folder structure is as follows:
    # results_dir
    #   - tool_1
    #       - virus
    #           - covreage
    #               - per_region
    #                   - ec_tool
    #                       - sample
    #                           - region
    #                               - output_file
    #               - whole_genome
    #                   - ec_tool
    #                       -sample
    #                           - region
    #                               - output_file

    results_dir = args.results
    hrt = args.hrt
    coverage = args.coverage
    ec = args.ec
    sample = args.sample
    ref_seq = args.ref_seq
    virus = args.virus
    primers = args.primers

    # read the genomic sequence of the reference
    ref_seq = str(SeqIO.read(ref_seq, "fasta").seq.upper())

    if (primers == "None"): 
        genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                    (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                    (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                        (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                    (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                        (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                        (26766, 27872), (27808, 28985), (28699, 29768)]
    else: 
        genomic_regions = get_genomic_regions(primers, mode = 2)

    # get the output directory for the tool
    outptut_dir = os.path.join(results_dir, hrt, virus, coverage, "whole_genome", ec, sample)

    if hrt == "cliquesnv":
        process_cliquesnv_output(outptut_dir, genomic_regions, ref_seq, ec)
    elif hrt == "haplodmf":
        process_haplodmf_output(outptut_dir, genomic_regions, ref_seq)    
    elif hrt == "rvhaplo":
        process_rvhaplo_output(outptut_dir, genomic_regions, ref_seq)
    else:
        raise ValueError("Invalid haplotype reconstruction tool")


if __name__ == '__main__':
    main()