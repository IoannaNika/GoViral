import argparse
import sys
import os

def main(): 
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--input_fasta', dest = 'input_fasta', required=True, type=str, help="")
    parser.add_argument('--primers_file', dest = 'primers_file', required=True, type=str, help="")
    parser.add_argument('--strategy', dest = 'strategy', required=True, type=str, help="pacbio-hifi or ONT")
    parser.add_argument('--outdir', dest = 'outdir', required=True, type=str, help="The output directory")
    parser.add_argument('--n', dest = 'n', required=True, type=int, help="Number of tuples to create")
    args = parser.parse_args()

    os.system(f"python -m long_amplicon_read_simulation.organise_input_fasta --input_fasta {args.input_fasta}")

    data_dir = "/".join(args.input_fasta.split("/")[:-1])
    data_dir = os.path.join(data_dir, "processed")

    os.system(f"python -m long_amplicon_read_simulation.make_amplicon_templates --dir {data_dir} --primers {args.primers_file} --input_fasta {args.input_fasta}")

    os.system(f"python -m long_amplicon_read_simulation.simulate_amplicon_reads --dir {data_dir} --strategy {args.strategy}") # at this point the .fastq files contain the simulated long amplicon reads
    
    outdir = os.path.join(args.outdir, "tuples_dataset")
    os.system(f"mkdir -p {outdir}")

    os.system(f"python -m long_amplicon_read_simulation.make_tuple_training_set --n {args.n} --out_dir {outdir} --data_dir {data_dir} --primers {args.primers_file} --strategy {args.strategy} --primer_file_mode 2")  

if __name__ == "__main__":
    sys.exit(main())