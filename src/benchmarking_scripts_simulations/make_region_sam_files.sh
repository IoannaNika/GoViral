for coverage in 100; do
      
    for n_haps in "n_1" "n_5" "n_10" "n_15" "n_20"; do

        for ec_tool in "hicanu" "original"; do

            sam_file="../data/Simulations/sarscov2/simulated_data/${coverage}/whole_genome/${ec_tool}/${n_haps}/sam_bam/reads_${ec_tool}.sam"
            primers="../data/Simulations/sarscov2/merged.bed"
            output="../data/Simulations/sarscov2/simulated_data/${coverage}/regions/${ec_tool}/${n_haps}/sam_bam"
            
            mkdir -p $output
            python -m benchmarking_scripts_simulations.make_region_sam_files --sam $sam_file --primers $primers --output $output

            done
        
        done
    
    done

done


for coverage in 100; do

    for virus in "hiv1" "hcv1b95"; do
      
        for sample in "ab_3_97" "ab_30_70" "ab_50_50" "ab_70_30" "ab_97_3"; do

            for ec_tool in "hicanu" "original"; do
            
                sam_file="../data/Simulations/${virus}/simulated_data/${coverage}/whole_genome/${ec_tool}/${sample}/sam_bam/reads_${ec_tool}.sam"
                primers="../data/Simulations/${virus}/primers.bed"
                output="../data/Simulations/${virus}/simulated_data/${coverage}/regions/${ec_tool}/${sample}/sam_bam"
                
                mkdir -p $output
                python -m benchmarking_scripts_simulations.make_region_sam_files --sam $sam_file --primers $primers --output $output

                done
            
            done

        done
    
    done

done