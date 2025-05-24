# sars cov 2
for coverage in 100; do
      
    for n_haps in "n_1" "n_5" "n_10" "n_15" "n_20"; do

        for ec_tool in "hicanu" "original"; do

            for hrt in "cliquesnv" "haplodmf" "rvhaplo"; do


                mixture_file="../data/Simulations/sarscov2/mixture_files/${n_haps}.json"
                input_path="results_simulations_2/${hrt}/sarscov2/${coverage}/whole_genome/${ec_tool}/${n_haps}/standard_output.tsv"
                templates="../data/Simulations/sarscov2/mixture_files/${n_haps}.json"
                output="results_simulations_2/results_sars_cov_2_tools.tsv"
                primers="../data/Simulations/sarscov2/merged.bed"
                data_dir="../data/Simulations/sarscov2/simulated_data/${coverage}/whole_genome/original/${n_haps}"
                
                if [[ ! -f "${input_path}" ]]; then
                    echo "${input_path} does not exist, skipping..."
                    continue
                fi

                echo "Processing sample ${input_path}"

                python -m benchmarking_scripts_simulations.evaluate_hrt_output_simulations --input_path $input_path --data_dir $data_dir --templates $templates --output $output --primers $primers --mixture_file $mixture_file --sample_name "${n_haps}-${ec_tool}-${hrt}-whole_genome"
            
            done

        done

    done

done
