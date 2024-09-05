import editdistance
from typing import Tuple
import os
from utils.utils import read_fasta_file

def calculate_average_edit_distance(edit_distance_from_closest_consensus_per_region:dict) -> float:
    """
    Calculates the average edit distance across all genomic regions and true haplotypes

    Args:
    edit_distance_from_closest_consensus_per_region: dict, edit distance from the closest consensus per region dictionary

    Returns:
    average_edit_distance: float, average edit distance
    """
    average_edit_distance_wuhan = 0
    average_edit_distance_omicron = 0

    for region in edit_distance_from_closest_consensus_per_region:
        average_edit_distance_wuhan += sum(edit_distance_from_closest_consensus_per_region[region]['Wuhan'])
        average_edit_distance_omicron += sum(edit_distance_from_closest_consensus_per_region[region]['Omicron'])
    
    average_edit_distance_wuhan = average_edit_distance_wuhan / len(edit_distance_from_closest_consensus_per_region.keys())
    average_edit_distance_omicron = average_edit_distance_omicron / len(edit_distance_from_closest_consensus_per_region.keys())

    average_edit_distance = round((average_edit_distance_wuhan + average_edit_distance_omicron) / 2, 2)
    
    return average_edit_distance


def calculate_average_number_of_haplotypes(number_of_haplotypes_per_region:dict) -> float:
    """
    Calculates the average number of haplotypes across all genomic regions and true haplotypes

    Args:
    number_of_haplotypes_per_region: dict, number of haplotypes per region dictionary

    Returns:
    average_number_of_haplotypes: float, average number of haplotypes
    """
    average_number_of_haplotypes = 0

    for region in number_of_haplotypes_per_region:
        average_number_of_haplotypes += number_of_haplotypes_per_region[region]['Wuhan'] + number_of_haplotypes_per_region[region]['Omicron']

    average_number_of_haplotypes = round(average_number_of_haplotypes / len(number_of_haplotypes_per_region.keys()), 2)

    return average_number_of_haplotypes


def calculate_recall(number_of_haplotypes_per_region:dict, true_number_haps_per_sample_region:dict, sample_name:str) -> Tuple[float, float]:
    """
    Calculates the recall for Wuhan and Omicron haplotypes. We define as recall the ratio of the number of true reconstructed haplotypes to the number of true haplotypes.
    How many times did we find the expected haplotypes?

    Args:
    number_of_haplotypes_per_region: dict, number of haplotypes per region dictionary
    true_number_haps_per_sample_region: dict, number of haplotypes for each sample per region
    sample_name: str, sample name

    Returns:
    recall_wuhan: float, recall for Wuhan haplotypes
    recall_omicron: float, recall for Omicron haplotypes
    """

    recall_wuhan = []
    recall_omicron = []

    for region in number_of_haplotypes_per_region:
        if true_number_haps_per_sample_region[sample_name][region]['Wuhan'] == 0:
            recall_wuhan_region = 1
        else:
            recall_wuhan_region = min (number_of_haplotypes_per_region[region]['Wuhan'] / true_number_haps_per_sample_region[sample_name][region]['Wuhan'], 1)
        if true_number_haps_per_sample_region[sample_name][region]['Omicron'] == 0:
            recall_omicron_region = 1
        else:
            recall_omicron_region = min (number_of_haplotypes_per_region[region]['Omicron'] / true_number_haps_per_sample_region[sample_name][region]['Omicron'], 1)

        recall_wuhan.append(recall_wuhan_region)
        recall_omicron.append(recall_omicron_region)

    recall_wuhan = round(sum(recall_wuhan) / len(recall_wuhan), 2)
    recall_omicron = round(sum(recall_omicron) / len(recall_omicron), 2)

    return recall_wuhan, recall_omicron


def calculate_duplication_ratio(number_of_haplotypes_per_region:dict, true_number_haps_per_sample_region:dict, sample_name:str) -> float:
    """
    Calculates the duplication ratio. The ratio of the number of reconstructed haplotypes to the number of true haplotypes.

    Args:
    number_of_haplotypes_per_region: dict, number of haplotypes per region dictionary
    true_number_haps_per_sample_region: dict, number of haplotypes for each sample per region
    sample_name: str, sample name

    Returns:
    duplication_ratio: float, duplication ratio
    """
    
    duplication_ratios = []
    
    for region in number_of_haplotypes_per_region:
        duplication_ratio = number_of_haplotypes_per_region[region]['Wuhan'] + number_of_haplotypes_per_region[region]['Omicron']
        duplication_ratio = duplication_ratio / sum(true_number_haps_per_sample_region[sample_name][region].values())
        duplication_ratios.append(duplication_ratio)

    duplication_ratio = round(sum(duplication_ratios) / len(duplication_ratios), 2)

    return duplication_ratio

def calculate_relative_absolute_abundance_error(abundance_per_region:dict, true_abundances_per_sample_per_region:dict, sample_name:str) -> float:
    """
    Calculates the average relative absolute abundance error

    Args:
    abundance_per_region: dict, abundance per region dictionary
    true_abundances_per_sample_per_region: dict, true abundances for each sample per region

    Returns:
    rel_ab_error: float, average relative absolute abundance error
    """

    rel_ab_error_wuhan = []
    rel_ab_error_omicron = []

    for region in abundance_per_region:
        rel_ab_error_wuhan.append(abs(abundance_per_region[region]['Wuhan'] - true_abundances_per_sample_per_region[sample_name][region]['Wuhan']) / (true_abundances_per_sample_per_region[sample_name][region]['Wuhan'] + 1e-6))
        rel_ab_error_omicron.append(abs(abundance_per_region[region]['Omicron'] - true_abundances_per_sample_per_region[sample_name][region]['Omicron']) / (true_abundances_per_sample_per_region[sample_name][region]['Omicron'] + 1e-6))

    rel_ab_error_wuhan = sum(rel_ab_error_wuhan) / len(abundance_per_region.keys())  
    rel_ab_error_omicron = sum(rel_ab_error_omicron) / len(abundance_per_region.keys())
       
    rel_ab_error = round((rel_ab_error_wuhan + rel_ab_error_omicron) / 2, 2)


    return rel_ab_error

def closest_haplotype(seq:str, rel_ab:float, omicron_amplicon:str, wuhan_amplicon:str, true_abundances_per_sample_per_region:dict, region:str, abundance_per_region:dict, sample_name:str) -> Tuple[str, int]:
    """
    Returns the closest haplotype to the input sequence

    Args:
    seq: str, input sequence
    rel_ab: float, relative abundance of the input sequence
    omicron_amplicon: str, omicron consensus sequence
    wuhan_amplicon: str, wuhan consensus sequence
    true_abundances_per_sample_per_region: dict, true abundances for each sample per region
    region: str, genomic region
    abundance_per_region: dict, current abundance per region
    sample_name: str, sample name

    Returns:
    haplotype: str, closest haplotype
    distance: int, edit distance between the input sequence and the closest haplotype
    """

    # compare the sequence with the wuhan and omicron consensus sequences
    wuhan_distance =  edit_distance_on_overlap(wuhan_amplicon, seq) #editdistance.eval(wuhan_amplicon, seq)
    omicron_distance = edit_distance_on_overlap(omicron_amplicon, seq)   #.eval(omicron_amplicon, seq)

    
    if wuhan_distance < omicron_distance:
        return 'Wuhan', wuhan_distance
    elif wuhan_distance > omicron_distance:
        return 'Omicron', omicron_distance
    else:
        updated_ab_wuhan = abundance_per_region[region]['Wuhan'] + rel_ab
        updated_ab_omicron = abundance_per_region[region]['Omicron'] + rel_ab
        try: 
            rel_ab_wuhan = abs(updated_ab_wuhan - true_abundances_per_sample_per_region[sample_name][region]['Wuhan']) / (true_abundances_per_sample_per_region[sample_name][region]['Wuhan'] + 1e-6)
        except: 
            print("KeyError: ", true_abundances_per_sample_per_region)
            raise KeyError

        rel_ab_omicron = abs(updated_ab_omicron - true_abundances_per_sample_per_region[sample_name][region]['Omicron']) / (true_abundances_per_sample_per_region[sample_name][region]['Omicron'] + 1e-6)

        if rel_ab_wuhan < rel_ab_omicron:
            return 'Wuhan', wuhan_distance
        else:
            return 'Omicron', omicron_distance


def edit_distance_on_overlap(seq1:str, seq2:str) -> Tuple[str, str]:
    # align the two sequences using mafft

    with open("temp.fasta", "w") as f:
        f.write(">seq1\n{}\n".format(seq1))
        f.write(">seq2\n{}\n".format(seq2))
        f.close()
    
    # align the sequences using mafft
    os.system("mafft --quiet temp.fasta > temp_aligned.fasta")

    # read the aligned sequences
    aligned_seqs = read_fasta_file("temp_aligned.fasta")

    # clean up
    os.system("rm temp.fasta temp_aligned.fasta")

    aligned_seq1 =""
    aligned_seq2 =""

    if aligned_seqs[0][0] == "seq1":
        aligned_seq1 = aligned_seqs[0][1]
        aligned_seq2 = aligned_seqs[1][1]
    else:
        aligned_seq1 = aligned_seqs[1][1]
        aligned_seq2 = aligned_seqs[0][1]
    
    # get the overlapping sequence
    start = 0

    while aligned_seq1[start] == "-" or aligned_seq2[start] == "-":
        start += 1
    
    end = -1

    while aligned_seq1[end] == "-" or aligned_seq2[end] == "-":
        end -= 1

    overlap_seq1 = aligned_seq1[start:end].replace("-", "").strip().strip("N")
    overlap_seq2 = aligned_seq2[start:end].replace("-", "").strip().strip("N")

    edit_distance = editdistance.eval(overlap_seq1, overlap_seq2)

    return edit_distance