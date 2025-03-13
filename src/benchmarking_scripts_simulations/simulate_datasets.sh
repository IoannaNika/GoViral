data_dir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations"
xsim_100=("ab_3_97" "ab_30_70" "ab_50_50" "ab_70_30" "ab_97_3")
xsim_1000=("ab_3_997" "ab_30_970" "ab_300_700" "ab_500_500" "ab_700_300" "ab_970_30" "ab_997_3")


for coverage in 100 1000; do
    if [ $coverage -eq 100 ]; then

        for sim_ab in "${xsim_100[@]}"; do
            
            for virus in "hcv1b95" "hiv1"; do

                input_fasta_data_dir="${data_dir}/${virus}/sequences.fasta"
                primer_file="${data_dir}/${virus}/primers.bed"
                mixture_file="${data_dir}/${virus}/mixture_files/${coverage}/${sim_ab}.json"
                outdir="${data_dir}/${virus}/simulated_data/${coverage}/original/${sim_ab}"

                mkdir -p $outdir

                python -m long_amplicon_read_simulation.run_benchmark_creation --input_fasta $input_fasta_data_dir --strategy "pacbio-hifi" --mixture $mixture_file --coverage $coverage --primers_file $primer_file --outdir $outdir
            done

        done

    elif [ $coverage -eq 1000 ]; then

        for sim_ab in "${xsim_1000[@]}"; do

            for virus in "hcv1b95" "hiv1"; do

                input_fasta_data_dir="${data_dir}/${virus}/sequences.fasta"
                primer_file="${data_dir}/${virus}/primers.bed"
                mixture_file="${data_dir}/${virus}/mixture_files/${coverage}/${sim_ab}.json"
                outdir="${data_dir}/${virus}/simulated_data/${coverage}/original/${sim_ab}"
                
                mkdir -p $outdir

                python -m long_amplicon_read_simulation.run_benchmark_creation --input_fasta $input_fasta_data_dir --strategy "pacbio-hifi" --mixture $mixture_file --coverage $coverage --primers_file $primer_file --outdir $outdir


            done

        done
    fi
done