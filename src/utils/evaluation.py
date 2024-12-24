import editdistance
from typing import Tuple
import os
from utils.utils import read_fasta_file
import uuid
import numpy as np

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

    denom_wuhan = 0
    denom_omicron = 0

    for region in edit_distance_from_closest_consensus_per_region:
        average_edit_distance_wuhan += sum(edit_distance_from_closest_consensus_per_region[region]['Wuhan'])
        denom_wuhan += len(edit_distance_from_closest_consensus_per_region[region]['Wuhan'])
        average_edit_distance_omicron += sum(edit_distance_from_closest_consensus_per_region[region]['Omicron'])
        denom_omicron += len(edit_distance_from_closest_consensus_per_region[region]['Omicron'])
    
    average_edit_distance_wuhan = average_edit_distance_wuhan / denom_wuhan #len(edit_distance_from_closest_consensus_per_region.keys())
    average_edit_distance_omicron = average_edit_distance_omicron / denom_omicron #len(edit_distance_from_closest_consensus_per_region.keys())

    average_edit_distance = round((average_edit_distance_wuhan + average_edit_distance_omicron) / 2, 3)
    
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

    average_number_of_haplotypes = round(average_number_of_haplotypes / len(number_of_haplotypes_per_region.keys()), 3)

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
            recall_wuhan_region = np.nan
        else:
            recall_wuhan_region = min (number_of_haplotypes_per_region[region]['Wuhan'] / true_number_haps_per_sample_region[sample_name][region]['Wuhan'], 1)
        if true_number_haps_per_sample_region[sample_name][region]['Omicron'] == 0:
            recall_omicron_region = np.nan
        else:
            recall_omicron_region = min (number_of_haplotypes_per_region[region]['Omicron'] / true_number_haps_per_sample_region[sample_name][region]['Omicron'], 1)
            
        recall_wuhan.append(recall_wuhan_region)
        recall_omicron.append(recall_omicron_region)
    
    recall_wuhan = [x for x in recall_wuhan if str(x) != 'nan']
    recall_omicron = [x for x in recall_omicron if str(x) != 'nan']

    if (len(recall_wuhan)> 0):
        recall_wuhan = round(sum(recall_wuhan) / len(recall_wuhan), 3)
    else:
        recall_wuhan = np.nan

    if (len(recall_omicron)> 0):
        recall_omicron = round(sum(recall_omicron) / len(recall_omicron), 3)
    else:
        recall_omicron = np.nan

    return recall_wuhan, recall_omicron


def calculate_precision(number_of_assigned_haplotypes_per_region:dict, number_of_all_haplotypes_per_region:dict) -> Tuple[float, float]:
    """
    Calculates the precision for Wuhan and Omicron haplotypes. How many of the reconstructed haplotypes assigned to Wuhan and Omicron are exact matches to the true haplotypes?

    Args:
    number_of_assigned_haplotypes_per_region: dict, number of haplotypes assigned to Wuhan or Omicron per region (not discarded)
    number_of_all_haplotypes_per_region: dict, number of all reconstructed haplotypes per region

    Returns:
    precision_wuhan: float, precision for Wuhan haplotypes
    precision_omicron: float, precision for Omicron haplotypes
    """

    precision_wuhan = []
    precision_omicron = []

    for region in number_of_assigned_haplotypes_per_region:
        if number_of_all_haplotypes_per_region[region]['Wuhan'] == 0:
            precision_wuhan_region = np.nan
        else:
            precision_wuhan_region = number_of_assigned_haplotypes_per_region[region]['Wuhan'] / number_of_all_haplotypes_per_region[region]['Wuhan']
        if number_of_all_haplotypes_per_region[region]['Omicron'] == 0:
            precision_omicron_region = np.nan
        else:
            precision_omicron_region = number_of_assigned_haplotypes_per_region[region]['Omicron'] / number_of_all_haplotypes_per_region[region]['Omicron']
        
        precision_wuhan.append(precision_wuhan_region)
        precision_omicron.append(precision_omicron_region)
    
    precision_wuhan = [x for x in precision_wuhan if str(x) != 'nan']
    precision_omicron = [x for x in precision_omicron if str(x) != 'nan']
    
    if (len(precision_wuhan)> 0):
        precision_wuhan = round(sum(precision_wuhan) / len(precision_wuhan), 3)
    else:
        precision_wuhan = np.nan

    if (len(precision_omicron)> 0):
        precision_omicron = round(sum(precision_omicron) / len(precision_omicron), 3)
    else:
        precision_omicron = np.nan

    return precision_wuhan, precision_omicron


def calculate_f1_score(precision:float, recall:float) -> float:
    """
    Calculates the F1 score. The harmonic mean of precision and recall.

    Args:
    precision: float, precision
    recall: float, recall

    Returns:
    f1_score: float, F1 score
    """

    f1_score = round(2 * (precision * recall) / (precision + recall + 1e-6), 3)

    return f1_score
    


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

    duplication_ratio = round(sum(duplication_ratios) / len(duplication_ratios), 3)

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
       
    rel_ab_error = round((rel_ab_error_wuhan + rel_ab_error_omicron) / 2, 3)


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
    wuhan_distance =  edit_distance_on_overlap(wuhan_amplicon, seq)
    omicron_distance = edit_distance_on_overlap(omicron_amplicon, seq)

    if wuhan_distance < omicron_distance:
        return 'Wuhan', wuhan_distance
    elif wuhan_distance > omicron_distance:
        return 'Omicron', omicron_distance
    else:
        updated_ab_wuhan = abundance_per_region[region]['Wuhan'] + rel_ab
        updated_ab_omicron = abundance_per_region[region]['Omicron'] + rel_ab
        try: 
            rel_ab_err_wuhan = abs(updated_ab_wuhan - true_abundances_per_sample_per_region[sample_name][region]['Wuhan']) / (true_abundances_per_sample_per_region[sample_name][region]['Wuhan'] + 1e-6)
        except: 
            print("KeyError: ", true_abundances_per_sample_per_region)
            raise KeyError

        rel_ab_err_omicron = abs(updated_ab_omicron - true_abundances_per_sample_per_region[sample_name][region]['Omicron']) / (true_abundances_per_sample_per_region[sample_name][region]['Omicron'] + 1e-6)

        if rel_ab_err_wuhan < rel_ab_err_omicron:
            return 'Wuhan', wuhan_distance
        else:
            return 'Omicron', omicron_distance


def find_overlapping_region(seq1:str, seq2:str) -> Tuple[str, str]:
    """
    Finds the overlapping region between two sequences

    Args:
    seq1: str, input sequence 1
    seq2: str, input sequence 2

    Returns:
    overlap_seq1: str, overlapping region of sequence 1
    overlap_seq2: str, overlapping region of sequence 2
    """
    unique_id = str(uuid.uuid4())

    # align the two sequences using mafft
    with open("temp_{}.fasta".format(unique_id), "w") as f:
        f.write(">seq1\n{}\n".format(seq1))
        f.write(">seq2\n{}\n".format(seq2))
        f.close()
    
    # align the sequences using mafft
    os.system("mafft --quiet temp_{}.fasta > temp_aligned_{}.fasta".format(unique_id, unique_id))

    # read the aligned sequences
    aligned_seqs = read_fasta_file("temp_aligned_{}.fasta".format(unique_id))

    # clean up
    os.system("rm temp_{}.fasta temp_aligned_{}.fasta".format(unique_id, unique_id))

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

    return overlap_seq1, overlap_seq2


def edit_distance_on_overlap(seq1:str, seq2:str) -> Tuple[str, str]:
    """
    Calculates the edit distance between two sequences on the overlapping region

    Args:
    seq1: str, input sequence 1
    seq2: str, input sequence 2

    Returns:
    edit_distance: int, edit distance between the two sequences
    """
    overlap_seq1, overlap_seq2 = find_overlapping_region(seq1, seq2)

    edit_distance = editdistance.eval(overlap_seq1, overlap_seq2)

    return edit_distance


def normalised_edit_distance_on_overlap(seq1:str, seq2:str) -> int:
    """
    Calculates the normalized edit distance between two sequences on the overlapping region

    Args:
    seq1: str, input sequence 1
    seq2: str, input sequence 2

    Returns:
    normalised_edit_distance: float, normalized edit distance between the two sequences
    """

    overlap_seq1, overlap_seq2 = find_overlapping_region(seq1, seq2)

    edit_distance = editdistance.eval(overlap_seq1, overlap_seq2)

    normalized_edit_distance = round(edit_distance / len(overlap_seq1), 3)

    return normalized_edit_distance
    
def percent_identity_on_overlap(seq1:str, seq2:str) -> int:
 """
 Calculates the percent identity on the overlapping region between the two sequences

 Args:
 seq1: str, input sequence 1
 seq2: str, input sequence 2

 Returns:
 percent identity: float, percent identity between the two sequences
 """

 overlap_seq1, overlap_seq2 = find_overlapping_region(seq1, seq2)

 edit_distance = editdistance.eval(overlap_seq1, overlap_seq2)

 percent_identity = round(((len(overlap_seq1) - edit_distance) / len(overlap_seq1)), 3)

 return percent_identity

