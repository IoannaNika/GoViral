import argparse
import sys
import os

def main(): 
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--input_fasta', dest = 'input_fasta', required=True, type=str, help="")
    parser.add_argument('--mixture', dest = 'mixture', required=True, type=str, help="")
    parser.add_argument('--coverage', dest = 'coverage', required=True, type=str, help="")
    parser.add_argument('--primers_file', dest = 'primers_file', required=True, type=str, help="")
    parser.add_argument('--strategy', dest = 'strategy', required=True, type=str, help="pacbio-hifi or ONT")
    parser.add_argument('--outdir', dest = 'outdir', required=True, type=str, help="The output directory")
    args = parser.parse_args()

    if args.strategy != "pacbio-hifi" or args.strategy != "ONT": 
        print(f"{args.strategy} is an invalid option for strategy, must be pacbio-hifi or ONT")
        exit()

    os.system(f"python -m long_amplicon_read_simulation.organise_input_fasta --input_fasta {args.input_fasta}")
    os.system(f"python -m long_amplicon_read_simulation.simulate_amplicon_reads_ab_cov --data_dir {args.outdir} --strategy {args.strategy} --ids_and_ab {args.mixture} --coverage {args.coverage} --outdir {args.outdir} --primers {args.primers_file} --input_fasta {args.input_fasta}")
    all_reads_output = os.path.join(args.outdir, "output")
    
    if strategy == "pacbio-hifi":
        os.system(f"python -m long_amplicon_read_simulation.make_all_reads --data_dir {args.outdir} --out {all_reads_output}")
    else: 
        os.system(f"python -m long_amplicon_read_simulation.make_all_reads --data_dir {args.outdir} --out {all_reads_output} --ont")
    
    test_pairs_output = os.path.join(all_reads_output, "all_test_pairs.tsv")
    all_reads_file = os.path.join(all_reads_output, "all_reads.tsv")
    os.system(f"python -m long_amplicon_read_simulation.make_test_pairs --reads {all_reads_file} --output {test_pairs_output}")
    os.system(f"python -m long_amplicon_read_simulation.put_fastas_in_one_file --input_dir {all_reads_output}")

if __name__ == "__main__":
    sys.exit(main())