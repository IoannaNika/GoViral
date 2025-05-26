# sars cov 2
results_dir=/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/results_simulations_3
ref_seq_hcv1b="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hcv1b95/reference.fasta"
ref_seq_hiv1="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hiv1/reference.fasta"
primers_hcv1b="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hcv1b95/primers.bed"
primers_hiv1="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hiv1/primers.bed"

for coverage in 100; do
    
    for virus in "hcv1b95" "hiv1"; do

        if [[ "$virus" == "hcv1b95" ]]; then
            ref_seq="${ref_seq_hcv1b}"
            primers="${primers_hcv1b}"
        else
            ref_seq="${ref_seq_hiv1}"
            primers="${primers_hiv1}"
        fi
      
        for sample in "ab_3_97" "ab_30_70" "ab_50_50" "ab_70_30" "ab_97_3"; do

            for ec_tool in "hicanu" "original"; do

                for hrt in "cliquesnv" "haplodmf" "rvhaplo"; do

                    python -m benchmarking_scripts_simulations.standardise_output_tools_wg --results_dir $results_dir --hrt $hrt --ec $ec_tool --sample $sample --ref_seq $ref_seq --coverage $coverage --virus $virus --primers $primers
                
                done

            done

        done

    done

done