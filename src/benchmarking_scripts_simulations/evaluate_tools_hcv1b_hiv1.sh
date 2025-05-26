results_dir="results_simulations_3"
for coverage in 100; do
      
    for sample in "ab_3_97" "ab_30_70" "ab_50_50" "ab_70_30" "ab_97_3"; do
        
        for virus in "hiv1" "hcv1b95"; do

            for ec_tool in "hicanu" "original"; do

                for hrt in "cliquesnv" "haplodmf" "rvhaplo"; do


                    mixture_file="../data/Simulations/${virus}/mixture_files/${coverage}/${sample}.json"
                    input_path="${results_dir}/${hrt}/${virus}/${coverage}/whole_genome/${ec_tool}/${sample}/standard_output.tsv"
                    templates="../data/Simulations/${virus}/mixture_files/${sample}.json"
                    output="${results_dir}/results_hcv1b_hiv1.tsv"
                    primers="../data/Simulations/${virus}/primers.bed"
                    data_dir="../data/Simulations/${virus}/simulated_data/${coverage}/whole_genome/original/${sample}"
                    
                    if [[ ! -f "${input_path}" ]]; then
                        echo "${input_path} does not exist, skipping..."
                        continue
                    fi

                    echo "Processing sample ${input_path}"

                    python -m benchmarking_scripts_simulations.evaluate_hrt_output_simulations --input_path $input_path --data_dir $data_dir --templates $templates --output $output --primers $primers --mixture_file $mixture_file --sample_name "${sample}-${ec_tool}-${hrt}-whole_genome-${virus}"
                
                done

            done

        done

    done

done
