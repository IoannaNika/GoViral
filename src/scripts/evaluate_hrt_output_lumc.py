import argparse
import pandas as pd
from ..utils.utils import cut_amplicon
from Bio import SeqIO
import editdistance
import os
from typing import List, Tuple

def is_the_coverage_sufficient(reads: pd.DataFrame, gr: str, low_limit: int = 100) -> bool:
    """
    Checks if the coverage is sufficient for a given genomic region

    Args:
    reads: pd.DataFrame, reads file
    gr: str, genomic region
    low_limit: int, low limit for the coverage

    Returns:
    bool: True if the coverage is sufficient, False otherwise
    """

    start_gr = int(gr.split('_')[0])
    end_gr = int(gr.split('_')[1])
    reads_region = reads[(reads[2] == start_gr) & (reads[3] == end_gr)]
    reads_coverage = len(reads_region)
    if reads_coverage >= low_limit:
        return True
    else:
        return False

def get_true_ab_per_sample() -> dict:
    """
    Returns the true abundances for each sample
    """ 

    return {
        '01_100': {"Wuhan": 1, "Omicron": 0},
        '02_100': {"Wuhan": 1, "Omicron": 0},
        '03_50': {"Wuhan": 0.5, "Omicron": 0.5},
        '04_75': {"Wuhan": 0.75, "Omicron": 0.25},
        '05_90': {"Wuhan": 0.9, "Omicron": 0.1},
        '06_95': {"Wuhan": 0.95, "Omicron": 0.05},
        '07_98': {"Wuhan": 0.98, "Omicron": 0.02},
        '08_0': {"Wuhan": 0, "Omicron": 1},
        '09_0': {"Wuhan": 0, "Omicron": 1}
    }

def get_n_haps_per_sample() -> dict:
    """
    Returns the number of haplotypes for each sample
    """

    return {
        '01_100': {"Wuhan": 1, "Omicron": 0},
        '02_100': {"Wuhan": 1, "Omicron": 0},
        '03_50': {"Wuhan": 1, "Omicron": 1},
        '04_75': {"Wuhan": 1, "Omicron": 1},
        '05_90': {"Wuhan": 1, "Omicron": 1},
        '06_95': {"Wuhan": 1, "Omicron": 1},
        '07_98': {"Wuhan": 1, "Omicron": 1},
        '08_0': {"Wuhan": 0, "Omicron": 1},
        '09_0': {"Wuhan": 0, "Omicron": 1}
    }

def true_abundances_and_haplotypes_per_sample_per_region(genomic_regions: List[Tuple[int, int]] , wuhan_consensus:str, omicron_consensus:str, ref_seq:str) -> Tuple[dict, dict]:
    """
    Returns the true abundances and number of haplotypes for each sample per region

    Args:
    genomic_regions: list of tuples, genomic regions
    wuhan_consensus: str, wuhan consensus sequence
    omicron_consensus: str, omicron consensus sequence
    ref_seq: str, reference sequence

    Returns:
    true_abundances_per_sample_per_region: dict, true abundances for each sample per region
    true_number_haps_per_sample_region: dict, number of haplotypes for each sample per region
    """

    true_abundances_per_sample_per_region = {}
    samples = ['01_100', '02_100' '03_50', '04_75', '05_90', "06_95", '07_98', '08_0', '09_0']
    true_number_haps_per_sample_region = {}

    for region in genomic_regions:

        wuhan_amplicon = cut_amplicon(wuhan_consensus, ref_seq, region[0], region[1])
        omicron_amplicon = cut_amplicon(omicron_consensus, ref_seq, region[0], region[1])
        
        if editdistance.eval(wuhan_amplicon, omicron_amplicon) > 0:

            for sample in samples:
                region_key = f"{region[0]}_{region[1]}"
                true_abundances_per_sample_per_region[sample] = {f"{region_key}": get_true_ab_per_sample()[sample]}
                true_number_haps_per_sample_region[sample] = {f"{region_key}": get_n_haps_per_sample()[sample]}

        else:

            for sample in samples:
                region_key = f"{region[0]}_{region[1]}"

                if sample in ['08_0', '09_0']:

                    true_abundances_per_sample_per_region[sample] = {f"{region_key}": {"Wuhan": 0, "Omicron": 1}}
                    true_number_haps_per_sample_region[sample] = {f"{region_key}": {"Wuhan": 0, "Omicron": 1}}
                else: 
                    true_abundances_per_sample_per_region[sample] = {f"{region_key}": {"Wuhan": 1, "Omicron": 0}}
                    true_number_haps_per_sample_region[sample] = {f"{region_key}": {"Wuhan": 1, "Omicron": 0}}

    return true_abundances_per_sample_per_region, true_number_haps_per_sample_region    

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
    wuhan_distance = editdistance.eval(wuhan_amplicon, seq)
    omicron_distance = editdistance.eval(omicron_amplicon, seq)

    
    if wuhan_distance < omicron_distance:
        return 'Wuhan', wuhan_distance
    elif wuhan_distance > omicron_distance:
        return 'Omicron', omicron_distance
    else:
        updated_ab_wuhan = abundance_per_region[region]['Wuhan'] + rel_ab
        updated_ab_omicron = abundance_per_region[region]['Omicron'] + rel_ab

        rel_ab_wuhan = abs(updated_ab_wuhan - true_abundances_per_sample_per_region[sample_name][region]['Wuhan']) / (true_abundances_per_sample_per_region[sample_name][region]['Wuhan'] + 1e-6)
        rel_ab_omicron = abs(updated_ab_omicron - true_abundances_per_sample_per_region[sample_name][region]['Omicron']) / (true_abundances_per_sample_per_region[sample_name][region]['Omicron'] + 1e-6)

        if rel_ab_wuhan < rel_ab_omicron:
            return 'Wuhan', wuhan_distance
        else:
            return 'Omicron', omicron_distance

def init_edit_distance_from_closest_consensus_per_region(genomic_regions:List[Tuple[int, int]]) -> dict:
    """
    Initializes the edit distance from the closest consensus per region dictionary

    Args:
    genomic_regions: list of tuples, genomic regions

    Returns:
    edit_distance_from_closest_consensus_per_region: dict, edit distance from the closest consensus per region dictionary
    """

    edit_distance_from_closest_consensus_per_region = {}

    for region in genomic_regions:
        edit_distance_from_closest_consensus_per_region[region] = {'Wuhan': [], 'Omicron': []}

    return edit_distance_from_closest_consensus_per_region

def update_edit_distance_from_closest_consensus_per_region(edit_distance_from_closest_consensus_per_region:dict, haplotype:str, region:str, distance:int) -> dict:
    """
    Updates the edit distance from the closest consensus per region dictionary

    Args:
    edit_distance_from_closest_consensus_per_region: dict, edit distance from the closest consensus per region dictionary
    haplotype: str, closest haplotype
    region: str, genomic region
    distance: int, edit distance between the input sequence and the closest haplotype

    Returns:
    edit_distance_from_closest_consensus_per_region: dict, updated edit distance from the closest consensus per region dictionary
    """

    edit_distance_from_closest_consensus_per_region[region][haplotype].append(distance)
    return edit_distance_from_closest_consensus_per_region

def init_number_of_haplotypes_per_region(genomic_regions:List[Tuple[int, int]]) -> dict:
    """
    Initializes the number of haplotypes per region dictionary

    Args:
    genomic_regions: list of tuples, genomic regions

    Returns:
    number_of_haplotypes_per_region: dict, number of haplotypes per region dictionary
    """

    number_of_haplotypes_per_region = {}

    for region in genomic_regions:
        number_of_haplotypes_per_region[region] = {'Wuhan': 0, 'Omicron': 0}

    return number_of_haplotypes_per_region


def update_number_of_haplotypes_per_region(number_of_haplotypes_per_region:dict, haplotype:str, region:str) -> dict:
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

def init_abundance_per_region(genomic_regions:List[Tuple[int, int]]) -> dict:
    """
    Initializes the abundance per region dictionary

    Args:
    genomic_regions: list of tuples, genomic regions

    Returns:
    abundance_per_region: dict, abundance per region dictionary
    """

    abundance_per_region = {}

    for region in genomic_regions:
        abundance_per_region[region]['Wuhan'] = 0
        abundance_per_region[region]['Omicron'] = 0

    return abundance_per_region

def update_abundance_per_region(abundance_per_region:dict, haplotype:str, region:str, rel_ab:float) -> dict:
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
        average_edit_distance_wuhan += edit_distance_from_closest_consensus_per_region[region]['Wuhan']
        average_edit_distance_omicron += edit_distance_from_closest_consensus_per_region[region]['Omicron']
    
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
            recall_wuhan_region = max (number_of_haplotypes_per_region[region]['Wuhan'] / true_number_haps_per_sample_region[sample_name][region]['Wuhan'], 1)
        if true_number_haps_per_sample_region[sample_name][region]['Omicron'] == 0:
            recall_omicron_region = 1
        else:
            recall_omicron_region = max (number_of_haplotypes_per_region[region]['Omicron'] / true_number_haps_per_sample_region[sample_name][region]['Omicron'], 1)

        recall_wuhan.append(recall_wuhan_region)
        recall_omicron.append(recall_omicron_region)

    recall_wuhan = sum(recall_wuhan) / len(recall_wuhan)
    recall_omicron = sum(recall_omicron) / len(recall_omicron)

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

    duplication_ratio = sum(duplication_ratios) / len(duplication_ratios)

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


def main():
    parser = argparse.ArgumentParser(description='Evaluate HRT output, standard format')
    parser.add_argument('--input', type=str, help='Input file')
    parser.add_argument('--ref_seq', type=str, help='Reference sequence')
    parser.add_argument('--output', type=str, help='Output file')
    parser.argument('--sample_name', type=str, help='Sample name')
    parser.add_argument('--reads', type=str, help='Reads file')
    args = parser.parse_args()

    reads_tsv = pd.read_csv(args.reads, sep='\t', header=None)

    wuhan_consensus_path = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/LUMC/consensus/wuhan.fasta"
    omicron_consensus_path = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/LUMC/consensus/omicron.fasta"

    wuhan_consensus = str(SeqIO.read(wuhan_consensus_path, "fasta").seq.upper())
    omicron_consensus = str(SeqIO.read(omicron_consensus_path, "fasta").seq.upper())

    true_abs_per_hap_per_sample_region, true_n_haps_per_sample_region = true_abundances_and_haplotypes_per_sample_per_region(genomic_regions, wuhan_consensus, omicron_consensus, ref_seq)

    ref_seq = args.ref_seq
    ref_seq = str(SeqIO.read(ref_seq, "fasta").seq.upper())
    input_path = args.input
    input_df = pd.read_csv(input_path, sep='\t', header=None)

    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                    (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                    (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                    (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]

    # filter genomic regions with insufficient coverage
    genomic_regions = [region for region in genomic_regions if is_the_coverage_sufficient(reads_tsv, f"{region[0]}_{region[1]}")]

    ##### metrics tracking #####
    edit_distance_from_closest_consensus_per_region = init_edit_distance_from_closest_consensus_per_region(genomic_regions)
    number_of_haplotypes_per_region = init_number_of_haplotypes_per_region(genomic_regions)
    abundance_per_region = init_abundance_per_region(genomic_regions)
    ###########################

    # iterate over the input file
    for index, row in input_df.iterrows():
        haplotype_id = row['haplotype_id']
        region = row['region']
        rel_ab = float(row['rel_abundance'])
        seq = row['sequence']

        # is it a region with sufficient coverage?
        if not is_the_coverage_sufficient(reads_tsv, region):
            # dont evaluate regions with insufficient coverage
            continue

        start = region.split(':')[0]
        end = region.split(':')[1]
        # extract wuhan and omicron consensus sequences for the specific genomic region
        wuhan_amplicon = cut_amplicon(wuhan_consensus, ref_seq, start, end)
        omicron_amplicon = cut_amplicon(omicron_consensus, ref_seq, start, end)

        closest_hap, ed_from_hap = closest_haplotype(seq, rel_ab, omicron_amplicon, wuhan_amplicon, true_abs_per_hap_per_sample_region, region, abundance_per_region, args.sample_name.split("-")[0])

        # update metrics
        edit_distance_from_closest_consensus_per_region = update_edit_distance_from_closest_consensus_per_region(edit_distance_from_closest_consensus_per_region, closest_hap, region, ed_from_hap)
        number_of_haplotypes_per_region = update_number_of_haplotypes_per_region(number_of_haplotypes_per_region, closest_hap, region)
        abundance_per_region = update_abundance_per_region(abundance_per_region, closest_hap, region, rel_ab)

    
    # calculate summary statistics for the whole sample
    average_edit_distance = calculate_average_edit_distance(edit_distance_from_closest_consensus_per_region)
    average_number_of_haplotypes = calculate_average_number_of_haplotypes(number_of_haplotypes_per_region)
    recall_omicron, recall_wuhan = calculate_recall(number_of_haplotypes_per_region, true_n_haps_per_sample_region,args.sample_name.split("-")[0])
    recall = (recall_omicron + recall_wuhan) / 2
    duplication_ratio = calculate_duplication_ratio(number_of_haplotypes_per_region, true_n_haps_per_sample_region, args.sample_name.split("-")[0])
    avg_rel_abs_ab_error = calculate_relative_absolute_abundance_error(abundance_per_region, true_abs_per_hap_per_sample_region, args.sample_name.split("-")[0])

    # if the output file does not exist, create it
    if not os.path.exists(args.output):
        with open(args.output, 'w') as f:
            f.write("sample_name\taverage_edit_distance\taverage_number_of_haplotypes\trecall\trecall_wuhan\trecall_omicron\tduplication_ratio\tavg_rel_abs_ab_error\n")
    
    # write the summary statistics to the output file
    with open(args.output, 'a') as f:
        f.write(f"{args.sample_name}\t{average_edit_distance}\t{average_number_of_haplotypes}\t{recall}\t{recall_wuhan}\t{recall_omicron}\t{duplication_ratio}\t{avg_rel_abs_ab_error}\n")

if __name__ == '__main__':
    main()
   