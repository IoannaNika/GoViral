import argparse
import os
import pandas as pd



def main():
    parser = argparse.ArgumentParser("Summarise evaluation across samples")
    parser.add_argument("--results_file", type=str, help="Path to the results file")
    parser.add_argument("--output_dir", type=str, help="Path to the output directory")
    args = parser.parse_args()

    results_file = args.results_file

    results = pd.read_csv(results_file, sep="\t", header=0)
    # sample_name	average_edit_distance	average_number_of_haplotypes	recall	recall_wuhan	recall_omicron	duplication_ratio	avg_rel_abs_ab_error
    # sample name example: 01_100-original-haplodmf-whole_genome

    # split sample name into parts based on -
    # if everything after the first dash is the same, we group by them 
    # the first part can vary becase is the sample which we want to summarise by

    results["sample_name"] = results["sample_name"].apply(lambda x: "-".join(x.split("-")[1:]))

    # group by the sample name
    grouped = results.groupby("sample_name")

    # for each group, calculate the mean of the metrics
    summary = grouped.mean()

    # save the summary
    summary.to_csv(os.path.join(args.output_dir, "result_summary.tsv"), sep="\t", index=True)
  
if __name__ == '__main__':
    main()