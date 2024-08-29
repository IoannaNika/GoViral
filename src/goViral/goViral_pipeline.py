import os
import pandas as pd
import argparse
import sys
from utils.utils import get_genomic_regions

def main(): 
    parser = argparse.ArgumentParser(description="Pipeline for running goViral")
    parser.add_argument('--directory', dest = 'directory', required=True, type=str, help="Directory for input and output files")
    parser.add_argument('--input_fastq', dest = 'input_fastq', required=True, type=str, help="Input file containing the fastq reads")
    parser.add_argument('--primers', dest = 'primers', required=True, type=str, help="File containing primer positions. To be used for mapping the amplicon reads to the correct genomic regions")
    parser.add_argument('--ref_seq', dest = 'ref_seq', required=True, type=str, help="Reference sequence. To be used for mapping the amplicon reads to the correct genomic regions")
    parser.add_argument('--coverage_limit', dest = 'coverage_limit', default=100, required=False, type=int, help="Coverage limit for subsampling. How many reads to consider in each subsample")
    parser.add_argument('--seed_limit', dest = 'seed_limit', required=False, default=100, type=int, help="Seed limit for subsampling. How many subsamples to consider")
    parser.add_argument('--follow_reccomendations', action='store_true', help="Follow reccomendations for seed limit. Overrides seed limit if set")
    parser.add_argument('--ab_threshold', dest = 'ab_threshold', required=False, default=0.1, type=float, help="Abundance threshold for filtering out low abundance sequences")
    
    args = parser.parse_args()

    os.system(f"python -m goViral.make_dataset --fastq {args.input_fastq} --primers {args.primers} --ref_seq {args.ref_seq} --out {args.directory}/reads.tsv")
    
    merged_file_path = f"{args.directory}/merged_standard_output.tsv"
    
    if os.path.exists(merged_file_path):
        os.remove(merged_file_path)

    with open(merged_file_path, 'w') as file:
        file.write('haplotype_id\tregion\trel_abundance\tsequence\n')  
        file.close()

    reads_tsv_file = os.path.join(args.directory, "reads.tsv")
    reads = pd.read_csv(reads_tsv_file, sep='\t')
    
    genomic_regions = get_genomic_regions(args.primers)

    for gr_start, gr_end in genomic_regions:

        seed = 1

        gr_specific_reads = reads[reads["start"] == gr_start]

        upper_limit = int(len(gr_specific_reads)/args.coverage_limit)

        if args.seed_limit > upper_limit:
            print("Not enough reads to reach seed limit. Seed limit set to: ", upper_limit)
            seed_limit = upper_limit
        else:
            print("Enough reads to reach seed limit. Seed limit set to: ", args.seed_limit)
            seed_limit = args.seed_limit

        # recommended limit for seed
        low_limit = int(len(gr_specific_reads)/args.coverage_limit)

        if args.seed_limit < low_limit and not args.follow_reccomendations:
            print("Seed limit might be too low. Recommended seed limit is: ", low_limit, " for the given coverage limit: ", args.coverage_limit, "and read count: ", len(reads), " at genomic region start: ", gr_start)

        if args.follow_reccomendations:
            print("Following reccomendations for seed limit, setting seed limit to: ", low_limit)
            seed_limit = low_limit

        while True: 

            seed_and_region_sepcific_output = f"{args.directory}/seed_{str(seed)}/subsampled_reads_{str(args.coverage_limit)}_{str(gr_start)}"
            
            os.system(f"python -m goViral.subsample_and_make_pairs --input_dir {args.directory} --seed {seed} --gr_start {gr_start} --coverage {args.coverage_limit}")
            os.system(f"python -m goViral.run_goViral_model --gr_start {gr_start} --outdir {seed_and_region_sepcific_output} --path_to_dataset {seed_and_region_sepcific_output}/all_test_pairs.tsv")
            os.system(f"python -m goViral.make_graph --results {seed_and_region_sepcific_output}/predictions.tsv --output {seed_and_region_sepcific_output}/communities.tsv")
            os.system(f"python -m goViral.make_consensus --communities {seed_and_region_sepcific_output}/communities.tsv --output {seed_and_region_sepcific_output}/standard_output.tsv")
            
            if seed >= seed_limit: 
                print("Seed limit reached")
                os.system(f"python -m goViral.process_seeds --directory {args.directory} --gr_start {gr_start} --seed_limit {seed} --ab_threshold {args.ab_threshold} --coverage {args.coverage_limit}")
                break
            
            seed +=1

if __name__ == "__main__":
    sys.exit(main())  