results_dir=/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/results_simulations_3
ref_seq="../data/LUMC/ref/nCoV-2019.reference.fasta"
primers="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/merged.bed"
genomic_regions=$(awk '/LEFT/ { left_end = $3 } /RIGHT/ { right_start = $2; print "(" left_end "," right_start ")" }' /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/merged.bed | paste -sd' ' -)
formatted_regions=$(echo "$genomic_regions" | sed -E 's/[()]//g' | sed -E 's/,/_/g')
virus="sarscov2"


for coverage in 100; do
    
    for region in $formatted_regions; do
    
        for sample in "n_1" "n_5" "n_10" "n_15" "n_20"; do

            for ec_tool in "hicanu" "original"; do

                for hrt in "cliquesnv" "haplodmf" "rvhaplo"; do

                    python -m benchmarking_scripts_simulations.standardise_output_tools_regions --results_dir $results_dir --hrt $hrt --ec $ec_tool --sample $sample --ref_seq $ref_seq --coverage $coverage --virus $virus --region $region
                    
                done
            done
        done

    done

done