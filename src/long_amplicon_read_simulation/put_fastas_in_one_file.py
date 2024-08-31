import os
import pandas as pd
import editdistance
import argparse



def main():

    parser = argparse.ArgumentParser(description="Make a fasta file with all reads")
    parser.add_argument('--input_dir', dest = 'input_dir', required=True, type=str, help="path to input file, fasta format")
    args = parser.parse_args()


    # iterate though all fasta files in the input directory
    input_dir = args.input_dir
    output_file = input_dir + "/reads.fasta"
    # remove the output file if it exists
    if os.path.exists(output_file):
        os.remove(output_file)

    # list the number of files in the input directory
    reads_dir = os.path.join(input_dir, "reads")
    num_files = len(os.listdir(reads_dir))
    print("Number of files in the input directory: ", num_files)

    cnt = 0

    lines = []
   
    # iterate through all files in the input directory
    for file in os.listdir(reads_dir):
        cnt += 1
        print("File: ", cnt, " out of ", num_files)

        with open(reads_dir + "/" + file, "r") as f:
            lines.extend(f.readlines())


    with open(output_file, "w") as f:
        f.writelines(lines)

    print("Done")

        
            

if __name__ == "__main__":
    main()
