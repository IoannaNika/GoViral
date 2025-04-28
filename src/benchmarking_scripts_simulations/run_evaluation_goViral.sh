primers="../data/Simulations/sarscov2/merged.bed"
output="results_simulations/goViral/sarscov2/100/whole_genome/original/results.tsv"

for n_haps in "n_1" "n_5" "n_10" "n_15" "n_20"; do

    mixture_file="../data/Simulations/sarscov2/mixture_files/${n_haps}.json"
    input_path="results_simulations/goViral/sarscov2/100/whole_genome/original/${n_haps}/merged_standard_output.tsv"
    data_dir="../data/Simulations/sarscov2/simulated_data/100/whole_genome/original/${n_haps}"
    templates="../data/Simulations/sarscov2/mixture_files/${n_haps}.json"

    python -m benchmarking_scripts_simulations.evaluate_hrt_output_simulations --input_path $input_path --data_dir $data_dir --templates $templates --output $output --primers $primers --mixture_file $mixture_file --sample_name "${n_haps}-original-goViral-per-region"

done 


xsim_100=("ab_3_97" "ab_30_70" "ab_50_50" "ab_70_30" "ab_97_3")
coverage=100

for virus in "hiv1", "hcv1b95"; do 

    primers="../data/Simulations/${virus}/primers.bed"
    output="results_simulations/goViral/${virus}/results.tsv"

    for sim_ab in "${xsim_100[@]}"; do
        mixture_file="../data/Simulations/${virus}/mixture_files/${coverage}/${sim_ab}.json"
        input_path="results_simulations/goViral/${virus}/original/${coverage}/${sim_ab}/merged_standard_output.tsv"
        data_dir="../data/Simulations/${virus}/simulated_data/${coverage}/original/${sim_ab}"
        templates="../data/Simulations/${virus}/mixture_files/${coverage}/${sim_ab}.json"
        python -m benchmarking_scripts_simulations.evaluate_hrt_output_simulations --input_path $input_path --data_dir $data_dir --templates $templates --output $output --primers $primers --mixture_file $mixture_file --sample_name "${sim_ab}-original-goViral-per-region"
    done

done