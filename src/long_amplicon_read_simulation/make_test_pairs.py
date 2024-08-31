import argparse
import os
import sys
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--reads', dest = 'reads', required=True, type=str, help="TSV file with single reads")
    parser.add_argument('--output', dest = 'output', required=True, type=str, help="Output file for read pairs")
    args = parser.parse_args()

    reads = args.reads
    output = args.output

    # if output file exists, delete it
    if os.path.exists(output):
        os.remove(output)
        # create new file
    with open(output, "w") as f:
        f.write("read_1\tread_2\tlabel\tgenomic_region\n")

    # read tsv file with list of reads
    df = pd.read_csv(reads, delimiter="\t")

    # get all unique combinations of reads from the same genomic region
    genomic_regions = df["genomic_region"].unique()

    # for each genomic region, get all unique combinations of reads
    for region in genomic_regions:
        reads = df[df["genomic_region"] == region]
        reads = reads["read_id"].unique()

        # get all unique combinations of reads
        for i in range(len(reads)):
            for j in range(i+1, len(reads)):
                with open(args.output, "a") as f:
                    if df[df["read_id"] == reads[i]]["label"].values[0] == df[df["read_id"] == reads[j]]["label"].values[0]:
                        f.write(reads[i] + "\t" + reads[j] + "\t" + "1" + "\t" + region + "\n")
                        f.close()
                    else:
                        f.write(reads[i] + "\t" + reads[j] + "\t" + "0" + "\t" + region + "\n")
                        f.close()
       
    

if __name__ == "__main__":
    sys.exit(main())