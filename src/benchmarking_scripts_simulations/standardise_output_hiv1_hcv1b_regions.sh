results_dir=/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/results_simulations_3
ref_seq_hcv1b="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hcv1b95/reference.fasta"
ref_seq_hiv1="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hiv1/reference.fasta"
primers_hcv1b="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hcv1b95/primers.bed"
primers_hiv1="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hiv1/primers.bed"
genomic_regions_hiv1=$(awk '/LEFT/ { left_end = $3 } /RIGHT/ { right_start = $2; print "(" left_end "," right_start ")" }' /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hiv1/primers.bed | paste -sd' ' -)
genomic_regions_hcv1b=$(awk '/LEFT/ { left_end = $3 } /RIGHT/ { right_start = $2; print "(" left_end "," right_start ")" }' /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hcv1b95/primers.bed | paste -sd' ' -)



for coverage in 100; do
    
    for virus in "hcv1b95"; do

        if [[ "$virus" == "hcv1b95" ]]; then
            ref_seq="${ref_seq_hcv1b}"
            primers="${primers_hcv1b}"
        else
            ref_seq="${ref_seq_hiv1}"
            primers="${primers_hiv1}"
        fi

        if [ "$virus" == "hiv1" ]; then
            genomic_regions=("${genomic_regions_hiv1[@]}")
        elif [ "$virus" == "hcv1b95" ]; then
                genomic_regions=("${genomic_regions_hcv1b[@]}")
        fi

        formatted_regions=$(echo "$genomic_regions" | sed -E 's/[()]//g' | sed -E 's/,/_/g')

        for region in $formatted_regions; do
    
            for sample in "ab_3_97" "ab_30_70" "ab_50_50" "ab_70_30" "ab_97_3"; do

                for ec_tool in "hicanu" "original"; do

                    for hrt in "cliquesnv" "haplodmf" "rvhaplo"; do

                        python -m benchmarking_scripts_simulations.standardise_output_tools_regions --results_dir $results_dir --hrt $hrt --ec $ec_tool --sample $sample --ref_seq $ref_seq --coverage $coverage --virus $virus --region $region
                    
                    done
                done
            done

        done

    done

done