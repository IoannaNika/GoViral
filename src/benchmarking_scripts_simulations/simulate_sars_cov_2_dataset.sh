data_dir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations"
dataset_mixtures=("n_1" "n_5" "n_10" "n_15" "n_20")
primer_file="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/merged.bed"
coverage=100

for dataset in "${dataset_mixtures[@]}"; do
    input_fasta_data_dir="${data_dir}/sarscov2/${dataset}/sequences.fasta"
    mixture_file="${data_dir}/sarscov2/mixture_files/${dataset}.json"
    outdir="${data_dir}/sarscov2/simulated_data/original/${dataset}"

    mkdir -p $outdir

    python -m long_amplicon_read_simulation.run_benchmark_creation --input_fasta $input_fasta_data_dir --strategy "pacbio-hifi" --mixture $mixture_file --coverage $coverage --primers_file $primer_file --outdir $outdir
done
