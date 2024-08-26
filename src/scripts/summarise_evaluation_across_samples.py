import argparse
import os
import pandas as pd



def main():
    parser = argparse.ArgumentParser("Summarise evaluation across samples")
    parser.add_argument("--results_file", type=str, help="Path to the results file")
    parser.add_argument("--output_file", type=str, help="Path to the output directory")
    parser.add_argument("--split_per_region_and_wg",action="store_true", help="If true, the results are split per region and whole genome")
    args = parser.parse_args()

    results_file = args.results_file

    results = pd.read_csv(results_file, sep="\t", header=0)
    # sample_name	average_edit_distance	average_number_of_haplotypes	recall	recall_wuhan	recall_omicron	duplication_ratio	avg_rel_abs_ab_error
    # sample name example: 01_100-original-haplodmf-whole_genome

    # split sample name into parts based on -
    # if everything after the first dash is the same, we group by them 
    # the first part can vary becase is the sample which we want to summarise by

    # replace in the names the first part with "pures" if the first part is "01_100", "02_100", "08_0", "09_0", otherwise replace it with "mixed"
    results["sample_name"] = results["sample_name"].apply(lambda x: "pures" + "-" + "-".join(x.split("-")[1:]) if x.startswith("01_100") or x.startswith("02_100") or x.startswith("08_0") or x.startswith("09_0") else "mixed" + "-" + "-".join(x.split("-")[1:]))

    # group by the sample name and keep the column "sample_name"
    grouped = results.groupby("sample_name")

    # for each group, calculate the mean of the metrics
    summary = grouped.mean()

    # make output directory if it doesn't exist
    output_path = "/".join(args.output_file.split("/")[:-1])
    print(output_path)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if args.split_per_region_and_wg:
        per_region_summary = summary[summary.index.str.contains("region")]
        # append _per_region to the output file name
        pr_file_name = args.output_file.split(".")[0] + "_per_region." + args.output_file.split(".")[1]
        per_region_summary.to_csv(pr_file_name, sep="\t", index=True)

        whole_genome_summary = summary[summary.index.str.contains("whole_genome")]

        # append _whole_genome to the output file name
        wg_file_name = args.output_file.split(".")[0] + "_whole_genome." + args.output_file.split(".")[1]
        whole_genome_summary.to_csv(wg_file_name, sep="\t", index=True)

    else: 
        # save the summary
        summary.to_csv(args.output_file, sep="\t", index=True)
  
if __name__ == '__main__':
    main()