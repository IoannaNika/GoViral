# sars cov 2

results_dir=/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/results_simulations_2
ref_seq="../data/LUMC/ref/nCoV-2019.reference.fasta"

for coverage in 100; do
      
    for n_haps in "n_1" "n_5" "n_10" "n_15" "n_20"; do

        for ec_tool in "hicanu" "original"; do

            for hrt in "cliquesnv" "haplodmf" "rvhaplo"; do

                python -m benchmarking_scripts_simulations.standardise_output_tools_wg --results_dir $results_dir --hrt $hrt --ec $ec_tool --sample $n_haps --ref_seq $ref_seq --coverage $coverage --virus "sarscov2"
            
            done

        done

    done

done