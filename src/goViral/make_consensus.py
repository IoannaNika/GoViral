import argparse
import sys
import pandas as pd
import igraph as ig
import random
import os
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
from utils.utils import check_and_update_if_haplotype_exists
from typing import IO, Any
import uuid

def make_output_file(output_file_name: str) -> IO[Any]:
    """
    Create the output file and write the header

    Args:
    output_file_name: The name of the output file

    Returns:
    file: The file object
    """    
    if os.path.exists(output_file_name):
        os.remove(output_file_name)
    
    # open file to write the consensus
    file = open(output_file_name, "w")
    file.write("haplotype_id" + "\t" + "region" + "\t" + "rel_abundance" + "\t" + "sequence" + "\n")
    file.close()
    return file

def mafft_alignment(sequences: list, unique_id: str) -> None:
    """
    Perform multiple sequence alignment using mafft

    Args:
    sequences: List of sequences to be aligned
    unique_id: unique_id to be used for temp files

    Returns:
    None
    """
    # clear any previous temp files
    if os.path.exists("temp_input_{}.fasta".format(unique_id)):
        os.remove("temp_input_{}.fasta".format(unique_id))
    if os.path.exists("temp_output_{}.fasta".format(unique_id)):
        os.remove("temp_output_{}.fasta".format(unique_id))

    with open("temp_input_{}.fasta".format(unique_id), "w") as file:
        for i, seq in enumerate(sequences):
            file.write(">" + str(i) + "\n")
            file.write(seq + "\n")
    # run mafft
    os.system("mafft --auto --quiet --thread 4 temp_input_{}.fasta > temp_output_{}.fasta".format(unique_id, unique_id))
    return


def create_consensus(community, genomic_region, results):
    # get consensus sequence
    sequences = []
    sequence_ids = []
    # get rows that correspond to the community and genomic region
    community_results = results[results["Community"] == community]
    community_results = community_results[community_results["Genomic_regions"] == genomic_region]

    unique_id = str(uuid.uuid4()) 
    # get sequences 
    for i in range(len(community_results)):
        sequences.append(community_results.iloc[i]["Sequence"])
        sequence_ids.append(community_results.iloc[i]["Sequence_id"])

    print("Number of sequences: ", len(sequences))
    print("Mafft alignment to be performed...")
    mafft_alignment(sequences, unique_id)
    print("Mafft alignment done.")
    alignment = AlignIO.read("temp_output_{}.fasta".format(unique_id), "fasta")  
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.gap_consensus(ambiguous='N', threshold=0.5)
    # to upper case
    consensus = consensus.upper()
    # # get string representation of the consensus
    consensus = str(consensus).replace("-", "")

    os.system("rm temp_output_{}.fasta temp_input_{}.fasta".format(unique_id, unique_id))
    return consensus

def get_results_per_community_and_genomic_region(community, genomic_region, results):
    results_per_community = results[results["Community"] == int(community)]
    results_per_community_and_gr = results_per_community[results_per_community["Genomic_regions"] == genomic_region]
    return results_per_community_and_gr

def main():
    parser = argparse.ArgumentParser(description="Create consensus")
    parser.add_argument('--communities', dest = 'results', required=True, type=str, help="tsv file with communities")
    parser.add_argument('--output', dest = 'output', required=True, type=str, help="output file")
    args = parser.parse_args()

    results = pd.read_csv(args.results, sep='\t', header=0)
    output = args.output

    file = make_output_file(output)

    # get genomic regions
    genomic_regions = results["Genomic_regions"].unique()

    # get communities per genomic region
    for genomic_region in genomic_regions:
        communities = results[results["Genomic_regions"] == genomic_region]
        communities = communities["Community"].unique()

        for community in communities:
            consensus = create_consensus(community, genomic_region, results)
            community = str(community)
            results_per_community_and_gr = get_results_per_community_and_genomic_region(community, genomic_region, results)

            number_of_sequences = len(results_per_community_and_gr)

            relative_abundance = number_of_sequences / len(results[results["Genomic_regions"] == genomic_region])

            if check_and_update_if_haplotype_exists(output, genomic_region, consensus, relative_abundance):
                continue
                                    
            if number_of_sequences >= 3: 
                with open(output, "a") as file:
                    file.write(community + "\t" + genomic_region + "\t" + str(relative_abundance) + "\t" + consensus + "\n")
        
    
    file.close()

if __name__ == "__main__":
    sys.exit(main())