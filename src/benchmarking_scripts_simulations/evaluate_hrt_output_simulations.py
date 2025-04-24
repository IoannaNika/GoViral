import argparse
import pandas as pd
import editdistance
import os
from typing import List, Tuple
import json
from utils.evaluation import percent_identity_on_overlap, calculate_f1_score
from utils.utils import get_genomic_regions
import numpy as np

def calculate_precision(number_of_exact_haplotypes_per_region: dict, number_of_assigned_haplotypes_per_region: dict):
    """
    Calculates the precision for multiple haplotypes in each region. Precision is defined as the ratio of 
    the number of true reconstructed haplotypes to the total number of assigned haplotypes. It is calculated 
    for each haplotype in each region.

    Args:
    number_of_exact_haplotypes_per_region: dict, number of exact haplotypes per region dictionary (true haplotypes)
    number_of_assigned_haplotypes_per_region: dict, number of assigned haplotypes per region dictionary

    Returns:
    precision_dict: dict, where keys are haplotype IDs and values are precision scores for that haplotype
    """
    precision_dict = {}
    
    # Loop over all regions
    for region in number_of_assigned_haplotypes_per_region.keys():

        for hap_id in number_of_assigned_haplotypes_per_region[region].keys():
           
            n_exact_haps = number_of_exact_haplotypes_per_region[region][hap_id]
            n_assigned_haps = number_of_assigned_haplotypes_per_region[region][hap_id]

            if n_assigned_haps == 0: 
                continue

            # print("Precision")
            # print("region : ", region)
            # print("hap_id : ", hap_id)
            # print("n_exact_haps : ", n_exact_haps)
            # print("n_assigned_haps : ", n_assigned_haps)
          
            precision_hap_region = min(n_exact_haps / n_assigned_haps, 1)

            # print("Precision_hap_region", precision_hap_region)

            # Add precision for this haplotype to the dictionary
            if hap_id not in precision_dict.keys():
                precision_dict[hap_id] = []
            precision_dict[hap_id].append(precision_hap_region)

            # print("precision_dict", precision_dict)

    # Compute mean precision for each haplotype
    for hap_id in precision_dict.keys():
        precision_dict[hap_id] = sum(precision_dict[hap_id]) / len(precision_dict[hap_id])

    return precision_dict

def calculate_recall(number_of_assigned_haplotypes_per_region: dict, true_num_of_haps_per_region: dict):
    """
    Calculates the recall for multiple haplotypes in each region. Recall is defined as the ratio of the 
    number of true reconstructed haplotypes to the number of true haplotypes. It is calculated for each
    haplotype in each region.

    Args:
    number_of_assigned_haplotypes_per_region: dict, number of haplotypes per region dictionary
    true_num_of_haps_per_region: dict, number of haplotypes for each sample per region

    Returns:
    recall_dict: dict, where keys are haplotype names and values are recall scores for that haplotype
    """
    recall_dict = {}
    
    # Loop over all regions
    for region in number_of_assigned_haplotypes_per_region:

        for hap_id in number_of_assigned_haplotypes_per_region[region]:

            if hap_id not in true_num_of_haps_per_region[region].keys():
                print(hap_id, " not in the true_num_of_haps_per_region")
                continue
           
            true_hap_count = true_num_of_haps_per_region[region][hap_id]
            n_found_haps = number_of_assigned_haplotypes_per_region[region][hap_id]
            
            if true_hap_count == 0: 
                print("for ", hap_id, " true_hap_count is 0")
                continue

            print("Recall")
            print("n_found_haps", n_found_haps)
            print("true_hap_count", true_hap_count)

            recall_hap_region = min(n_found_haps / true_hap_count, 1)

            print('recall_hap_region', recall_hap_region)

            # Add recall for this haplotype to the dictionary
            if hap_id not in recall_dict.keys():
                recall_dict[hap_id] = []
            recall_dict[hap_id].append(recall_hap_region)
    
    print("recall dict", recall_dict)
    # Compute mean recall for each haplotype
    for hap_id in recall_dict.keys():
        recall_dict[hap_id] = sum(recall_dict[hap_id]) / len(recall_dict[hap_id])
    
    print("recall dict")
    print(recall_dict)
    return recall_dict

def calculate_relative_absolute_abundance_error(abundance_per_region: dict, true_abs_per_hap_per_sample_region: dict):
    """
    Calculates the average relative absolute abundance error across all haplotypes. The error is computed 
    for each haplotype and averaged over all regions and haplotypes.

    Args:
    abundance_per_region: dict, abundance per region dictionary
    true_abs_per_hap_per_sample_region: dict, true abundances for each sample per region
    sample_name: str, the name of the sample to consider

    Returns:
    rel_ab_error: float, average relative absolute abundance error across all haplotypes
    """
    haplotype_errors = {}

    # Loop over all regions
    for region in abundance_per_region:
        for hap_id in abundance_per_region[region]:
            # Get the true abundance for the haplotype in this region and sample
            if hap_id in true_abs_per_hap_per_sample_region[region].keys():
                true_abundance = true_abs_per_hap_per_sample_region[region][hap_id]
            else:
                assert(abundance_per_region[region][hap_id] == 0)
                continue
                

            assigned_abundance = abundance_per_region[region][hap_id]

            # Calculate relative absolute abundance error, ensuring we avoid division by zero
            rel_ab_error = abs(assigned_abundance - true_abundance) / (true_abundance + 1e-6)

            if hap_id not in haplotype_errors:
                haplotype_errors[hap_id] = []
            
            haplotype_errors[hap_id].append(rel_ab_error)
    
    # average across regions per haplotype
    all_hap_errors = []
    for hap_id, errors in haplotype_errors.items():
        haplotype_avg_error = sum(errors) / len(errors)
        all_hap_errors.append(haplotype_avg_error)

    # average across haplotypes
    rel_ab_error = sum(all_hap_errors) / len(all_hap_errors)

    return rel_ab_error

def calculate_duplication_ratio(number_of_haplotypes_per_region: dict, true_num_of_haps_per_region: dict):
    dr = 0
    
    for region in number_of_haplotypes_per_region.keys(): 
        dr += sum(number_of_haplotypes_per_region[region].values()) / sum(true_num_of_haps_per_region[region].values())
    
    dr = dr/len( number_of_haplotypes_per_region.keys())
    
    return dr

def calculate_average_number_of_haplotypes(number_of_haplotypes_per_region: dict):
    avg = 0
    for region in number_of_haplotypes_per_region.keys(): 
        avg += sum(number_of_haplotypes_per_region[region].values())
    avg = avg / len(number_of_haplotypes_per_region.keys())
    return avg

def calculate_average_metric(metric_dictionary: dict) -> float:
    """
    Calculates the average metric across haplotypes and regions.
    
    Parameters:
    - metric_dictionary: A dictionary structured as {region -> hap -> []}
    
    Returns:
    - The overall average metric across all haplotypes and regions.
    """
    
    # Step 1: Calculate the average metric per haplotype across all regions
    haplotype_avg_ed = {}  # Dictionary to store average metric for each haplotype
    
    for region, hap_dict in metric_dictionary.items():
        for hap_id, edit_distances in hap_dict.items():
            # Calculate the average metric for the current haplotype in this region
            hap_avg_ed = sum(edit_distances) / len(edit_distances)
            
            # Add this region's average metric to the haplotype's running total
            if hap_id not in haplotype_avg_ed:
                haplotype_avg_ed[hap_id] = []
            haplotype_avg_ed[hap_id].append(hap_avg_ed)
    
    # Step 2: Calculate the average across haplotypes
    total_avg = 0
    num_haplotypes = 0
    
    for hap_id, hap_avg_ed_list in haplotype_avg_ed.items():
        # Calculate the average of averages for this haplotype
        hap_total_avg = sum(hap_avg_ed_list) / len(hap_avg_ed_list)
        total_avg += hap_total_avg
        num_haplotypes += 1
    
    # Step 3: Calculate the overall average across all haplotypes
    overall_avg = total_avg / num_haplotypes if num_haplotypes > 0 else 0
    
    return overall_avg

def closest_haplotype(seq: str, templates: dict, rel_ab: float, true_abundances_per_region: dict, abundance_per_region: dict, region: str):
    """
    Assigns seq to the closest haplotype based on edit distance and relative abundance error.
    
    Parameters:
    - seq: The sequence to be assigned to a haplotype.
    - templates: Dictionary of haplotypes with keys as haplotype IDs and values as haplotype sequences.
    - rel_ab: The relative abundance of the haplotype being considered.
    - true_abundances_per_region: Nested dictionary of true abundances per region and haplotype.
    - abundance_per_region: Nested dictionary of current relative abundances per region and haplotype.
    - region: The specific region to which the sequence belongs.
    
    Returns:
    - closest_hap: The closest haplotype sequence based on edit distance and relative abundance.
    - ed_from_hap: The edit distance of the assigned haplotype.
    """

    # Step 1: Calculate edit distance to each haplotype in 'templates'
    edit_distances = {hap_id: editdistance.eval(seq, hap[region]) for hap_id, hap in templates.items() if hap_id in true_abundances_per_region[region].keys()}
    
    # Step 2: Find the minimum edit distance(s)
    min_ed = min(edit_distances.values())
    closest_haps = [hap_id for hap_id, ed in edit_distances.items() if ed == min_ed and hap_id in true_abundances_per_region[region]]
    
    # Step 3: If there's more than one closest haplotype, minimize the relative abundance error
    if len(closest_haps) > 1:
        for hap_id in closest_haps:
            abundance_per_region[region][hap_id] += rel_ab  # Adjust abundance for the region
        
            # Calculate relative abundance errors for each of the closest haplotypes
            rel_ab_errors = {
                hap_id: abs(abundance_per_region[region][hap_id] - true_abundances_per_region[region][hap_id])/(true_abundances_per_region[region][hap_id] + 1e-6)
            }
        
        # Step 4: Select the haplotype with the minimal relative abundance error
        min_ab_error_idx = min(rel_ab_errors, key=rel_ab_errors.get)
        closest_hap = min_ab_error_idx
        ed_from_hap = min_ed
    else:
        # If there's only one closest haplotype, no tie-breaking is needed
        closest_hap = closest_haps[0]
        ed_from_hap = min_ed
    
    return closest_hap, ed_from_hap

def load_json(file_path: str):
    with open(file_path, 'r') as file:
        return json.load(file)

def get_templates(mixture_json: dict, data_dir: str):
    """
    Returns a dictionary the templates (true haplotypes) for each sequence. 

    Args:
    mixture_file: JSON containing the mixture seq identifier and simulated abundance
    data_dir: directory containing the simulated data per sample
    """ 
    templates_names = get_seq_ids(mixture_json)

    templates_files = [(os.path.join(data_dir, f"{name}.template"), name) for name in templates_names]
    result = {}

    for templates_file, seq_id in templates_files: 
    
        with open(templates_file, 'r') as file:
            genomic_region = None
            sequence = None
            
            for line in file:
                line = line.strip()

                if line.startswith('>'):
                    # Parse the header to extract seq_id, genomic_region, and num
                    header = line[1:]  # Skip '>'
                    parts = header.split(':')
                    seq_id = parts[0][2:]

                    genomic_region = ':'.join(parts[1:2])
                    num = int(parts[-1])  # we use num to identify which sequence
                    
                    # Initialize the dictionary structure if not yet
                    if seq_id not in result:
                        result[seq_id] = {}
                    
                    # Prepare for the next sequence line
                    sequence = ""
                    
                    # If num is 0, initialize genomic_region in the dictionary
                    if num == 0:
                        if genomic_region not in result[seq_id]:
                            result[seq_id][genomic_region] = sequence
                            
                else:
                    # Add the sequence to the current sequence
                    sequence = line
                    # Store sequence when num == 0
                    if num == 0:
                        result[seq_id][genomic_region] += sequence

    return result

def get_seq_ids(mixture_json: dict):
    """
    Returns the sequence identifiers

    Args:
    mixture_json: JSON containing the mixture seq identifier and simulated abundance
    """ 

    return list(mixture_json.keys())

def get_n_haps(mixture_json: dict):
    """
    Returns the number of haplotypes for the sample

    Args:
    mixture_file: JSON containing the mixture seq identifier and simulated abundance
    """

    return len(mixture_json.keys())

def true_abundances_and_haplotypes_per_region(genomic_regions: List[str] , templates_seq: dict, templates_ab: dict):
    """
    Returns the true abundances and number of haplotypes per region for the sample
    
    Args:
    genomic_regions: list of strings, genomic regions
    templates_seq: dictionary with seq: ACTG-template sequence
    templates_ab: dictionary with seq: true abundance
    """

    true_abundances_per_region = {}
    true_num_of_haps_per_region = {}


    for region in genomic_regions:
        region = str(region[0]) + "_" +  str(region[1])
        true_n_haps = 0

        unique_seqs = {}
        region_abs = {}
        true_num_of_haps_per_region[region] = {}

        for seq_identifier, region_seq_dict in templates_seq.items():
            
            seq = region_seq_dict[region]

            if seq not in unique_seqs.keys(): 
                unique_seqs[seq] = seq_identifier
                region_abs[seq_identifier] = templates_ab[seq_identifier]
                true_num_of_haps_per_region[region][seq_identifier] = 1
            else: 
                existing_seq_id = unique_seqs[seq]
                region_abs[existing_seq_id] += templates_ab[seq_identifier]
                # true_num_of_haps_per_region[region] = {seq_identifier: 0}

        # true_num_of_haps_per_region[region] = len(unique_seqs.keys())
        true_abundances_per_region[region] = region_abs

    return true_abundances_per_region, true_num_of_haps_per_region    

def update_edit_distance_from_closest_hap_per_region(edit_distance_from_closest_hap_per_region:dict, haplotype:str, region:str, distance:int):
    """
    Updates the edit distance from the closest consensus per region dictionary

    Args:
    edit_distance_from_closest_hap_per_region: dict, edit distance from the closest consensus per region dictionary
    haplotype: str, closest haplotype
    region: str, genomic region
    distance: int, edit distance between the input sequence and the closest haplotype

    Returns:
    edit_distance_from_closest_hap_per_region: dict, updated edit distance from the closest consensus per region dictionary
    """
    edit_distance_from_closest_hap_per_region[region][haplotype].append(distance)

    return edit_distance_from_closest_hap_per_region

def init_dict_per_region(genomic_regions:List[str], seq_ids: List[str], init_value):
    """
    Initializes the dictionaries tracking the number of haplotypes, abundance, discard haplotypes per region dictionary

    Args:
    genomic_regions: list of strings, genomic regions
    seq_ids: sequence identifiers of the haplotypes
    init_value: value with which to initalize the dictionary

    Returns:
    initialized_dict: dict, an initialised dict
    """

    return_dict = {}

    for region in genomic_regions:
        
        region = str(region[0]) + "_" + str(region[1])

        if seq_ids != None: 
            return_dict[region] = {seq_id: init_value for seq_id in seq_ids}
        else: 
            return_dict[region] = init_value

    return return_dict

def update_number_of_haplotypes_per_region(number_of_haplotypes_per_region:dict, haplotype:str, region:str):
    """ 
    Updates the number of haplotypes per region dictionary

    Args:
    number_of_haplotypes_per_region: dict, number of haplotypes per region dictionary
    haplotype: str, closest haplotype
    region: str, genomic region

    Returns:
    number_of_haplotypes_per_region: dict, updated number of haplotypes per region dictionary
    """
    
    number_of_haplotypes_per_region[region][haplotype] += 1
    
    return number_of_haplotypes_per_region

def update_abundance_per_region(abundance_per_region:dict, haplotype:str, region:str, rel_ab:float):
    """
    Updates the abundance per region dictionary

    Args:
    abundance_per_region: dict, abundance per region dictionary
    haplotype: str, closest haplotype
    region: str, genomic region
    rel_ab: float, relative abundance of the input sequence

    Returns:
    abundance_per_region: dict, updated abundance per region dictionary
    """

    abundance_per_region[region][haplotype] += rel_ab
    return abundance_per_region

def normalize_abundance_per_region(abundance_per_region:dict) -> dict:
    """
    Normalizes the abundance per region such that the sum of the relative abundances of the Wuhan and Omicron haplotypes is 1

    Args:
    abundance_per_region: dict, abundance per region dictionary

    Returns:
    abundance_per_region: dict, normalized abundance per region dictionary
    """

    for region in abundance_per_region.keys():
        total_ab = sum(abundance_per_region[region].values())
        if total_ab == 0: 
            print("No abundance for region: ", region)
            print("Abundance per region: ", abundance_per_region)
            continue
        
        for seq_id in abundance_per_region[region].keys():
            abundance_per_region[region][seq_id] = abundance_per_region[region][seq_id] / total_ab 
        
        total_ab_check = sum(abundance_per_region[region].values())
        print("Total abundance check: ", total_ab_check)
        assert total_ab_check > 0.99

    return abundance_per_region

def main():
    parser = argparse.ArgumentParser(description='Evaluate HRT output, standard format')
    parser.add_argument('--input_path', type=str, help='')
    parser.add_argument('--data_dir', type=str, help='')
    parser.add_argument('--templates', type=str, help='')
    parser.add_argument('--output', type=str, help='Output file')
    parser.add_argument('--primers', type=str, help='')
    parser.add_argument('--mixture_file', type=str, help='')
    parser.add_argument('--sample_name', type=str, help='Sample name')
    args = parser.parse_args()

    print("Processing sample: ", args.sample_name)
    mixture_json = load_json(args.mixture_file)

    input_df = pd.read_csv(args.input_path, sep='\t', header=0)
    templates = get_templates(mixture_json, args.data_dir)
    print("templates")
    print(templates)
    seq_ids = get_seq_ids(mixture_json)

    # defines the threshold for which a haplotype ends up in the discarded haplotypes
    dissimilarity_cutoff = 0.2

    # ensure that for each genomic region relative abundance sums up to 1
    for region in input_df['region'].unique():
        region_df = input_df[input_df['region'] == region]
        total_rel_ab = region_df['rel_abundance'].sum()
        if total_rel_ab != 1:
            # iterate over the rows of the region and update the relative abundance
            for index, row in input_df.iterrows():
                if row['region'] == region:
                    input_df.at[index, 'rel_abundance'] = row['rel_abundance'] / total_rel_ab
    
    genomic_regions = get_genomic_regions(args.primers, mode = 2)

    true_abundances_per_region, true_num_of_haps_per_region = true_abundances_and_haplotypes_per_region(genomic_regions, templates, mixture_json)
    
    # print("true_num_of_haps_per_region", true_num_of_haps_per_region)
    # print("true_abundances_per_region", true_abundances_per_region)
    ##### metrics tracking #####
    discarded_haplotypes_per_region = init_dict_per_region(genomic_regions, None, 0)
    edit_distance_from_closest_hap_per_region = init_dict_per_region(genomic_regions, seq_ids, [])
    percent_identity_from_closest_hap_per_region = init_dict_per_region(genomic_regions, seq_ids, [])
    number_of_all_haplotypes_per_region = init_dict_per_region(genomic_regions, seq_ids, 0)
    number_of_haplotypes_per_region = init_dict_per_region(genomic_regions, seq_ids, 0)
    number_of_exact_haplotypes_per_region = init_dict_per_region(genomic_regions, seq_ids, 0)
    abundance_per_region = init_dict_per_region(genomic_regions, seq_ids, 0)
    ###########################

    # iterate over the input file
    for index, row in input_df.iterrows():

        haplotype_id = row.iloc[0]
        region = row.iloc[1]
        rel_ab = float(row.iloc[2])
        seq = row.iloc[3]

        closest_hap, ed_from_hap = closest_haplotype(seq, templates, rel_ab, true_abundances_per_region, abundance_per_region, region)
        print(closest_hap, " is the closest hap with edit distance ", ed_from_hap," with rel ab ", rel_ab,  " region ", region)
        
        percent_identity = percent_identity_on_overlap(seq, templates[closest_hap][region])
  
        cutoff_percent_identity = (1-percent_identity) * 100

        print(f"Cutoff percent identity {cutoff_percent_identity}")
        number_of_all_haplotypes_per_region = update_number_of_haplotypes_per_region(number_of_all_haplotypes_per_region, closest_hap, region)
        # update metrics
        if cutoff_percent_identity > dissimilarity_cutoff:
            print("Discarding haplotype ", haplotype_id, " with edit distance ", ed_from_hap, " region ", region, "and length ", len(seq))
            # if the normalised edit distance is higher than the cutoff, the haplotype is discarded
            print(f"Cutoff percent identity {cutoff_percent_identity}")
            discarded_haplotypes_per_region[region] += 1
            continue

        edit_distance_from_closest_hap_per_region = update_edit_distance_from_closest_hap_per_region(edit_distance_from_closest_hap_per_region, closest_hap, region, ed_from_hap)
        percent_identity_from_closest_hap_per_region = update_edit_distance_from_closest_hap_per_region(percent_identity_from_closest_hap_per_region, closest_hap, region, percent_identity)
        number_of_haplotypes_per_region = update_number_of_haplotypes_per_region(number_of_haplotypes_per_region, closest_hap, region)
        abundance_per_region = update_abundance_per_region(abundance_per_region, closest_hap, region, rel_ab)

        if ed_from_hap == 0: 
            number_of_exact_haplotypes_per_region = update_number_of_haplotypes_per_region(number_of_exact_haplotypes_per_region, closest_hap, region)

    # normalize the abundance per region to account for the discarded haplotypes
    abundance_per_region = normalize_abundance_per_region(abundance_per_region)
    print("Abundance per region: ", abundance_per_region)
    # calculate summary statistics for the whole sample
    average_edit_distance = calculate_average_metric(edit_distance_from_closest_hap_per_region)
    average_percent_identity = calculate_average_metric(percent_identity_from_closest_hap_per_region)
    average_number_of_haplotypes = calculate_average_number_of_haplotypes(number_of_haplotypes_per_region)
    recalls = calculate_recall(number_of_haplotypes_per_region, true_num_of_haps_per_region)
    recalls = [recall for seq, recall in recalls.items()]
    recall =  sum(recalls) / len(recalls)
    precisions = calculate_precision(number_of_exact_haplotypes_per_region, number_of_haplotypes_per_region)
    precisions = [precision for seq, precision in precisions.items()]
    precision =  sum(precisions) / len(precisions)
    f1_score = calculate_f1_score(precision, recall)
    duplication_ratio = calculate_duplication_ratio(number_of_haplotypes_per_region, true_num_of_haps_per_region)
    average_number_of_haplotypes_discard = sum(discarded_haplotypes_per_region.values()) / len(discarded_haplotypes_per_region.keys())
    avg_rel_abs_ab_error = calculate_relative_absolute_abundance_error(abundance_per_region, true_abundances_per_region)
    
    # if the output file does not exist, create it
    if not os.path.exists(args.output):
        with open(args.output, 'w') as f:
            f.write("sample_name\tf1_score\trecall\tprecision\taverage_edit_distance\tavg_num_haps_in_discard\tavg_rel_abs_ab_error\taverage_percent_identity\tduplication_ratio\n")
            f.close()
    
    # write the summary statistics to the output file
    with open(args.output, 'a') as f:
        f.write(f"{args.sample_name}\t{round(f1_score,3)}\t{round(recall,3)}\t{round(precision,3)}\t{round(average_edit_distance,3)}\t{round(average_number_of_haplotypes_discard, 3)}\t{round(avg_rel_abs_ab_error,3)}\t{round(average_percent_identity,3)}\t{round(duplication_ratio,3)}\n")
        f.close()

if __name__ == '__main__':
    main()
   