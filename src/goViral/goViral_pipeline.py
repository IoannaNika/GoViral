import os
import pandas as pd
import argparse
import sys

def main(): 
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--directory', dest = 'directory', required=True, type=str, help="")
    parser.add_argument('--input_fastq', dest = 'input_fastq', required=True, type=str, help="Input file containing the fastq reads")
    parser.add_argument('--primers', dest = 'primers', required=True, type=str, help="File containing primer positions. To be used for mapping the amplicon reads to the correct genomic regions")
    parser.add_argument('--ref_seq', dest = 'ref_seq', required=True, type=str, help="Reference sequence. To be used for mapping the amplicon reads to the correct genomic regions")
    parser.add_argument('--gr_start', dest = 'gr_start', required=True, type=int, help="")
    parser.add_argument('--coverage_limit', dest = 'coverage_limit', default=100, required=False, type=int, help="")
    parser.add_argument('--seed_limit', dest = 'seed_limit', required=False, default=100, type=int, help="")
    parser.add_argument('--ab_threshold', dest = 'ab_threshold', required=False, default=0.1, type=float, help="")
    
    args = parser.parse_args()

    # os.system(f"python -m goViral.make_dataset --fastq {args.input_fastq} --primers {args.primers} --ref_seq {args.ref_seq} --out {args.directory}/reads.tsv")
    
    merged_file_path = f"{args.directory}/merged_standard_output.tsv"
    
    if os.path.exists(merged_file_path):
        os.remove(merged_file_path)

    with open(merged_file_path, 'w') as file:
        file.write('haplotype_id\tregion\trel_abundance\tsequence\n')  
        file.close()
    
    output_df = pd.read_csv(merged_file_path, sep = "\t", header = 0)
    
    current_output_length = len(output_df)
    seed = 1

    while True: 

        seed_and_region_sepcific_output = f"{args.directory}/seed_{str(seed)}/subsampled_reads_{str(args.coverage_limit)}_{str(args.gr_start)}"
        
        os.system(f"python -m goViral.subsample_and_make_pairs --input_dir {args.directory} --seed {seed} --gr_start {args.gr_start} --coverage {args.coverage_limit}")
        os.system(f"python -m goViral.run_goViral_model --primers {args.primers} --outdir {seed_and_region_sepcific_output} --path_to_dataset {seed_and_region_sepcific_output}/all_test_pairs.tsv")
        os.system(f"python -m goViral.make_graph --results {seed_and_region_sepcific_output}/predictions.tsv --output {seed_and_region_sepcific_output}/communities.tsv")
        os.system(f"python -m goViral.make_consensus --communities {seed_and_region_sepcific_output}/communities.tsv --output {seed_and_region_sepcific_output}/standard_output.tsv")
        os.system(f"python -m goViral.process_seeds --directory {args.directory} --gr_start {args.gr_start} --seed_limit {seed} --ab_threshold {args.ab_threshold} --coverage {args.coverage_limit}")

        output_df = pd.read_csv(merged_file_path, sep = "\t", header = 0)
        new_output_length = len(output_df)

        # if new_output_length == current_output_length:
        #     print("No new sequences identified")
        #     break
        
        if seed >= args.seed_limit: 
            print("Seed limit reached")
            break
        
        seed +=1
        current_output_length = len(output_df)

if __name__ == "__main__":
    sys.exit(main())  