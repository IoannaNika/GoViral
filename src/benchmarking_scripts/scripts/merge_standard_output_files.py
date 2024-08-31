import argparse
import os
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='Merge standard output files')
    parser.add_argument('--results_dir', dest='results', type=str, required=True, help='Results directory')
    args = parser.parse_args()

    # create the output file 
    outfile_path = os.path.join(args.results, "merged_standard_output.tsv")
    outfile = open(outfile_path, "w")
    outfile.close()


    # get the list of all files that are named standard_output.tsv
    # within the results directory merge all the files into one file
    # save the file as merged_standard_output.tsv in the results director

    results_path = Path(args.results)
    all_files = list(results_path.rglob("standard_output.tsv")) 

    for file_path in all_files:
       # read tsv file
        df = pd.read_csv(file_path, sep="\t", header=0)
        # if the outfile is empty, include the header in the output file
        # else do not include the header
        if os.path.getsize(outfile_path) == 0:
            df.to_csv(outfile_path, sep="\t", index=False)
        else:
            df.to_csv(outfile_path, sep="\t", index=False, header=False, mode='a')
             
if __name__ == '__main__':
    main()


    
