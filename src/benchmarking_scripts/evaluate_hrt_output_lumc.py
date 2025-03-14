import argparse
import pandas as pd
from utils.utils import cut_amplicon
from Bio import SeqIO
import editdistance
import os
from typing import List, Tuple
from utils.evaluation import calculate_average_edit_distance, calculate_average_number_of_haplotypes, calculate_recall, calculate_duplication_ratio, closest_haplotype, calculate_relative_absolute_abundance_error, normalised_edit_distance_on_overlap, calculate_precision, calculate_f1_score, percent_identity_on_overlap

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
        print(start_gr, end_gr, "Enough coverage")
        return True
    else:
        print(start_gr, end_gr, "Not enough coverage")
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
    samples = ['01_100', '02_100', '03_50', '04_75', '05_90', "06_95", '07_98', '08_0', '09_0']
    true_number_haps_per_sample_region = {}


    for region in genomic_regions:

        wuhan_amplicon = cut_amplicon(wuhan_consensus, ref_seq, region[0], region[1])
        omicron_amplicon = cut_amplicon(omicron_consensus, ref_seq, region[0], region[1])
        
        if editdistance.eval(wuhan_amplicon, omicron_amplicon) > 0:

            for sample in samples:

                region_key = f"{region[0]}_{region[1]}"

                if sample not in true_abundances_per_sample_per_region.keys():
                    true_abundances_per_sample_per_region[sample] = {}
                if sample not in true_number_haps_per_sample_region.keys():
                    true_number_haps_per_sample_region[sample] = {}
                    
                true_abundances_per_sample_per_region[sample].update({f"{region_key}": get_true_ab_per_sample()[sample]})
                true_number_haps_per_sample_region[sample].update({f"{region_key}": get_n_haps_per_sample()[sample]})

        else:

            for sample in samples:

                region_key = f"{region[0]}_{region[1]}"

                if sample not in true_abundances_per_sample_per_region.keys():
                    true_abundances_per_sample_per_region[sample] = {}
                if sample not in true_number_haps_per_sample_region.keys():
                    true_number_haps_per_sample_region[sample] = {}

                if sample in ['08_0', '09_0']:
                    true_abundances_per_sample_per_region[sample].update({f"{region_key}": {"Wuhan": 0, "Omicron": 1}})
                    true_number_haps_per_sample_region[sample].update({f"{region_key}": {"Wuhan": 0, "Omicron": 1}})
                else: 
                    true_abundances_per_sample_per_region[sample].update({f"{region_key}": {"Wuhan": 1, "Omicron": 0}})
                    true_number_haps_per_sample_region[sample].update({f"{region_key}": {"Wuhan": 1, "Omicron": 0}})
    
    print("True abundances per sample per region")
    print(true_abundances_per_sample_per_region)
    
    return true_abundances_per_sample_per_region, true_number_haps_per_sample_region    


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
        gr_key = str(region[0]) + "_" + str(region[1])
        edit_distance_from_closest_consensus_per_region[gr_key] = {'Wuhan': [], 'Omicron': []}

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
        gr_key = str(region[0]) + "_" + str(region[1])
        number_of_haplotypes_per_region[gr_key] = {'Wuhan': 0, 'Omicron': 0}

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
        gr_key = str(region[0]) + "_" + str(region[1])
        abundance_per_region[gr_key] = {"Wuhan": 0, "Omicron":0}

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

def init_discarded_haplotypes_per_region(genomic_regions:List[Tuple[int, int]]) -> dict:
    """
    Initializes the discarded haplotypes per region dictionary. The dictionary contains the number of discarded haplotypes per region.
    Discarded haplotypes are the ones that are not sufficiently similar to the Wuhan or Omicron consensus sequences.

    Args:
    genomic_regions: list of tuples, genomic regions

    Returns:
    discarded_haplotypes_per_region: dict, discarded haplotypes per region dictionary
    """

    discarded_haplotypes_per_region = {}

    for region in genomic_regions:
        gr_key = str(region[0]) + "_" + str(region[1])
        discarded_haplotypes_per_region[gr_key] = 0

    return discarded_haplotypes_per_region

def normalize_abundance_per_region(abundance_per_region:dict) -> dict:
    """
    Normalizes the abundance per region such that the sum of the relative abundances of the Wuhan and Omicron haplotypes is 1

    Args:
    abundance_per_region: dict, abundance per region dictionary

    Returns:
    abundance_per_region: dict, normalized abundance per region dictionary
    """

    for region in abundance_per_region.keys():
        total_ab = abundance_per_region[region]['Wuhan'] + abundance_per_region[region]['Omicron']
        if total_ab == 0: 
            print("No abundance for region: ", region)
            print("Abundance per region: ", abundance_per_region)
            continue
        abundance_per_region[region]['Wuhan'] = abundance_per_region[region]['Wuhan'] / total_ab 
        abundance_per_region[region]['Omicron'] = abundance_per_region[region]['Omicron'] / total_ab
        
        print("Total abundance check: ", (abundance_per_region[region]['Wuhan'] + abundance_per_region[region]['Omicron']))
        assert (abundance_per_region[region]['Wuhan'] + abundance_per_region[region]['Omicron']) > 0.99

    return abundance_per_region

def main():
    parser = argparse.ArgumentParser(description='Evaluate HRT output, standard format')
    parser.add_argument('--input', type=str, help='Input file')
    parser.add_argument('--ref_seq', type=str, help='Reference sequence')
    parser.add_argument('--output', type=str, help='Output file')
    parser.add_argument('--sample_name', type=str, help='Sample name')
    parser.add_argument('--reads', type=str, help='Reads file')
    args = parser.parse_args()

    print("Processing sample: ", args.sample_name)

    reads_tsv = pd.read_csv(args.reads, sep='\t', header=None)

    wuhan_consensus_path = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/LUMC/consensus/wuhan.fasta"
    omicron_consensus_path = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/LUMC/consensus/omicron.fasta"

    wuhan_consensus = str(SeqIO.read(wuhan_consensus_path, "fasta").seq.upper())
    omicron_consensus = str(SeqIO.read(omicron_consensus_path, "fasta").seq.upper())

    ref_seq = args.ref_seq
    ref_seq = str(SeqIO.read(ref_seq, "fasta").seq.upper())
    input_path = args.input
    input_df = pd.read_csv(input_path, sep='\t', header=0)

    # defines the threshold for which a haplotype ends up in the discarded haplotypes
    dissimilarity_cutoff = 0.2

    # remove sequences that are shorter than 800 bp
    # input_df = input_df[input_df['sequence'].apply(lambda x: len(x) > 800)]

    # ensure that for each genomic region relative abundance sums up to 1
    for region in input_df['region'].unique():
        region_df = input_df[input_df['region'] == region]
        total_rel_ab = region_df['rel_abundance'].sum()
        if total_rel_ab != 1:
            # iterate over the rows of the region and update the relative abundance
            for index, row in input_df.iterrows():
                if row['region'] == region:
                    input_df.at[index, 'rel_abundance'] = row['rel_abundance'] / total_rel_ab

    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                    (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                    (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                    (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]

    # filter genomic regions with insufficient coverage
    genomic_regions = [region for region in genomic_regions if is_the_coverage_sufficient(reads_tsv, f"{region[0]}_{region[1]}")]

    true_abs_per_hap_per_sample_region, true_n_haps_per_sample_region = true_abundances_and_haplotypes_per_sample_per_region(genomic_regions, wuhan_consensus, omicron_consensus, ref_seq)
    
    ##### metrics tracking #####
    discarded_haplotypes_per_region = init_discarded_haplotypes_per_region(genomic_regions)
    edit_distance_from_closest_consensus_per_region = init_edit_distance_from_closest_consensus_per_region(genomic_regions)
    percent_identity_from_closest_consensus_per_region = init_edit_distance_from_closest_consensus_per_region(genomic_regions)
    number_of_assigned_haplotypes_per_region = init_number_of_haplotypes_per_region(genomic_regions)
    number_of_all_haplotypes_per_region = init_number_of_haplotypes_per_region(genomic_regions)
    number_of_exact_haplotypes_per_region = init_number_of_haplotypes_per_region(genomic_regions)
    abundance_per_region = init_abundance_per_region(genomic_regions)
    ###########################

    # iterate over the input file
    for index, row in input_df.iterrows():
        haplotype_id = row.iloc[0]
        region = row.iloc[1]
        rel_ab = float(row.iloc[2])
        seq = row.iloc[3]

        # is it a region with sufficient coverage?
        if not is_the_coverage_sufficient(reads_tsv, region):
            # dont evaluate regions with insufficient coverage
            continue

        start = int(region.split('_')[0])
        end = int(region.split('_')[1])
        
        # extract wuhan and omicron consensus sequences for the specific genomic region
        wuhan_amplicon = cut_amplicon(wuhan_consensus, ref_seq, start, end)
        omicron_amplicon = cut_amplicon(omicron_consensus, ref_seq, start, end)

        closest_hap, ed_from_hap = closest_haplotype(seq, rel_ab, omicron_amplicon, wuhan_amplicon, true_abs_per_hap_per_sample_region, region, abundance_per_region, args.sample_name.split("-")[0])
        print(closest_hap, " is the closest hap with edit distance ", ed_from_hap," with rel ab ", rel_ab,  " region ", region)
        
        if closest_hap == 'Wuhan':
            percent_identity = percent_identity_on_overlap(seq, wuhan_amplicon)
        else:
            percent_identity = percent_identity_on_overlap(seq, omicron_amplicon)
        
        cutoff_percent_identity = (1-percent_identity) * 100
        print(f"Cutoff percent identity {cutoff_percent_identity}")
        number_of_all_haplotypes_per_region = update_number_of_haplotypes_per_region(number_of_all_haplotypes_per_region, closest_hap, region)
        # update metrics
        if cutoff_percent_identity > dissimilarity_cutoff:
            try:
                print("Discarding haplotype ", haplotype_id, " with edit distance ", ed_from_hap, " region ", region, "and length ", len(seq))
            except:
                print("Discarding haplotype ", haplotype_id, " with edit distance ", ed_from_hap, " region ", region, "and seq ", seq)
            # if the normalised edit distance is higher than the cutoff, the haplotype is discarded
            print(f"Cutoff percent identity {cutoff_percent_identity}")
            discarded_haplotypes_per_region[region] += 1
            continue

        edit_distance_from_closest_consensus_per_region = update_edit_distance_from_closest_consensus_per_region(edit_distance_from_closest_consensus_per_region, closest_hap, region, ed_from_hap)
        percent_identity_from_closest_consensus_per_region = update_edit_distance_from_closest_consensus_per_region(percent_identity_from_closest_consensus_per_region, closest_hap, region, percent_identity)
        number_of_assigned_haplotypes_per_region = update_number_of_haplotypes_per_region(number_of_assigned_haplotypes_per_region, closest_hap, region)
        abundance_per_region = update_abundance_per_region(abundance_per_region, closest_hap, region, rel_ab)
    
        if ed_from_hap == 0: 
            number_of_exact_haplotypes_per_region = update_number_of_haplotypes_per_region(number_of_exact_haplotypes_per_region, closest_hap, region)

    # normalize the abundance per region to account for the discarded haplotypes
    abundance_per_region = normalize_abundance_per_region(abundance_per_region)
    print("Abundance per region: ", abundance_per_region)
    # calculate summary statistics for the whole sample
    average_edit_distance = calculate_average_edit_distance(edit_distance_from_closest_consensus_per_region)
    average_percent_identity = calculate_average_edit_distance(percent_identity_from_closest_consensus_per_region)
    average_number_of_haplotypes = calculate_average_number_of_haplotypes(number_of_assigned_haplotypes_per_region)
    recall_wuhan, recall_omicron = calculate_recall(number_of_assigned_haplotypes_per_region, true_n_haps_per_sample_region, args.sample_name.split("-")[0])
    recalls = [recall_wuhan, recall_omicron]
    recalls = [x for x in recalls if str(x) != 'nan']
    recall = round(sum(recalls)/len(recalls),3)
    precision_omicron, precision_wuhan = calculate_precision(number_of_assigned_haplotypes_per_region = number_of_exact_haplotypes_per_region, number_of_all_haplotypes_per_region = number_of_assigned_haplotypes_per_region)
    precisions = [precision_omicron, precision_wuhan]
    precisions = [x for x in precisions if str(x) != 'nan']
    precision = round(sum(precisions)/len(precisions),3)
    f1_score = calculate_f1_score(precision, recall)
    duplication_ratio = calculate_duplication_ratio(number_of_assigned_haplotypes_per_region, true_n_haps_per_sample_region, args.sample_name.split("-")[0])
    average_number_of_haplotypes_discard = sum(discarded_haplotypes_per_region.values()) / len(discarded_haplotypes_per_region.keys())
    avg_rel_abs_ab_error = calculate_relative_absolute_abundance_error(abundance_per_region, true_abs_per_hap_per_sample_region, args.sample_name.split("-")[0])
    
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
   