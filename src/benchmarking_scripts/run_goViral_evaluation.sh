results_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/results"
reference="../data/LUMC/ref/nCoV-2019.reference.fasta"
output="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/results/results.tsv"

for sample in "01_100" "02_100" "03_50" "04_75" "05_90" "06_95" "07_98" "08_0" "09_0"; do   
    
    input="${results_path}/goViral/sample_${sample}/merged_standard_output.tsv"
    reads="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data/natural_mixtures/${sample}/reads.tsv"
    sample_name="${sample}-original-goViral-per_region"
    
    python -m benchmarking_scripts.evaluate_hrt_output_lumc --input $input --ref_seq $reference --output $output --sample_name $sample_name --reads $reads
done