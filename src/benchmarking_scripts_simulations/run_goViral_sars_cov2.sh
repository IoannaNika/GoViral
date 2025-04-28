accelerator="gpu"
virus="sarscov2"
out_dir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src"
data_dir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations"
genomic_regions_sars_cov_2=$(awk '$4 ~ /_LEFT/ {print $3}' /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/merged.bed)

for coverage in 100; do
      
    for n_haps in "n_1" "n_5" "n_10" "n_15" "n_20"; do
        
        merged_outdir="${out_dir}/results_simulations/goViral/${virus}/${coverage}/whole_genome/original/${n_haps}"
        old_output_file="${merged_outdir}/merged_standard_output.tsv"
        
        if [ -f "$old_output_file" ] ; then
            rm "$old_output_file"
        fi  

        for start_gr in $genomic_regions_sars_cov_2; do 

            outdir="${out_dir}/results_simulations/goViral/${virus}/${coverage}/whole_genome/original/${n_haps}/gr_${start_gr}"
            datadir="${data_dir}/${virus}/simulated_data/${coverage}/whole_genome/original/${n_haps}/output/all_test_pairs.tsv"
            mkdir -p $outdir
     
            python  -m goViral.run_goViral_model --outdir $outdir --path_to_dataset $datadir --gr_start $start_gr --accelerator $accelerator
            python -m goViral.make_graph --results ${outdir}/predictions.tsv --output ${outdir}/communities.tsv
            python -m goViral.make_consensus --communities  ${outdir}/communities.tsv --output  ${outdir}/standard_output.tsv
            python -m benchmarking_scripts_simulations.process_regions_simulations --directory $merged_outdir --gr_start $start_gr --ab_threshold 0
            
        done
    
    done
   
done

