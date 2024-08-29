import argparse
import sys
import os
import pandas as pd

def write_sequence_to_file(full_read_id, sequence, outdir_read_dir):
    # check if file exists
    if os.path.exists(outdir_read_dir + "/" + full_read_id + ".fasta"):
        os.remove(outdir_read_dir + "/" + full_read_id + ".fasta")

    if not os.path.exists(outdir_read_dir + "/" + full_read_id + ".fasta"):
        with open(outdir_read_dir + "/" + full_read_id + ".fasta", "w") as f:
            f.write(">" + full_read_id + "\n")
            f.write(sequence + "\n")
    return

def main(): 
    parser = argparse.ArgumentParser(description="Subsample reads and make all pairs")
    parser.add_argument('--input_dir', dest = 'input_dir', required=True, type=str, help="directory containing the reads.tsv file corresponding to the original dataset")
    parser.add_argument('--seed', dest = 'seed', required=True, type=int, help="seed random sampling of reads, with replacement")
    parser.add_argument('--gr_start', dest = 'gr_start', required=True, type=int, help="datasets_for_genomic_region")
    parser.add_argument('--coverage', dest = 'coverage', default=100, required=False, type=int, help="coverage: how many reads to keep")
    args = parser.parse_args()

    reads_tsv_file = args.input_dir + "/reads.tsv"
    reads = pd.read_csv(reads_tsv_file, sep='\t')

    output_dir = args.input_dir + "/seed_" + str(args.seed) + "/subsampled_reads_" + str(args.coverage) + "_" + str(args.gr_start)

    if os.path.exists(output_dir):
        os.system("rm -r " + output_dir)

    os.makedirs(output_dir)

    output_file = output_dir + "/reads.tsv"
    # make output file with 'w' mode
    open(output_file, "w").close()

    reads.columns = ["read_id", "sequence", "start", "end", "strand"]

    reads = reads[reads["start"] == args.gr_start]

    if len(reads) > args.coverage:
        # subsample the reads and leave only 100 reads per genomic region
        reads = reads.groupby(["start", "end"]).sample(n=args.coverage, random_state=args.seed).reset_index(drop=True)
    
    reads.to_csv(output_file, sep='\t', index=False)

    # for all reads in the reads.tsv file, write the read to a file
    for index, row in reads.iterrows():
        read_id = row["read_id"].replace("/", "_")
        sequence = row["sequence"]
        start = row["start"]
        end = row["end"]
        strand = row["strand"]

        if strand == "-":
            sequence = sequence[::-1].translate(str.maketrans("ATGC", "TACG"))
       
        seq_output_dir = output_dir + "/reads"

        if not os.path.exists(seq_output_dir):
            os.makedirs(seq_output_dir)
        
        write_sequence_to_file(read_id, sequence, seq_output_dir)

    all_test_pairs = output_dir + "/all_test_pairs.tsv"
    # make output file with 'w' mode
    open(all_test_pairs, "w").close()

    with open(all_test_pairs, "w") as f:
        f.write("id1\tid2\tstart\tend\n")
    
    # for all reads in the reads.tsv file, make all pairs
    for index, row in reads.iterrows():
        read_id = row["read_id"].replace("/", "_")
        sequence = row["sequence"]
        start = row["start"]
        end = row["end"]
        strand = row["strand"]

        for index2, row2 in reads.iterrows():
            read_id2 = row2["read_id"].replace("/", "_")
            sequence2 = row2["sequence"]
            start2 = row2["start"]
            end2 = row2["end"]
            strand2 = row2["strand"]

            if index2 <= index:
                continue

            if start != start2:
                continue

            if read_id == read_id2:
                continue
            
            if strand2 == "-":
                sequence2 = sequence2[::-1].translate(str.maketrans("ATGC", "TACG"))
            if strand == "-":
                sequence = sequence[::-1].translate(str.maketrans("ATGC", "TACG"))
           
            with open(all_test_pairs, "a") as f:
                f.write(read_id + "\t" + read_id2 + "\t" + str(start) + "\t" + str(end) +"\n")
                f.close()


if __name__ == "__main__":
    main()

