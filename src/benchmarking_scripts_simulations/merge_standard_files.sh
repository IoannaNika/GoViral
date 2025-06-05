
for virus in "hiv1" "hcv1b95"; do 
    
    for sample in "ab_3_97" "ab_30_70" "ab_50_50" "ab_70_30" "ab_97_3"; do

        for ec_tool in "hicanu" "original"; do
        
            for hrt_tool in "cliquesnv" "haplodmf" "rvhaplo"; do 

                # merge region standardised tsv files called "standard_output.tsv" into one file
                python -m benchmarking_scripts.merge_standard_output_files --results_dir results_simulations_3/${hrt_tool}/${virus}/100/regions/${ec_tool}/${sample}/
            
            done
        done
    done
done


for virus in "sarscov2"; do 
    
    for sample in "n_1" "n_5" "n_10" "n_15" "n_20"; do

        for ec_tool in "hicanu" "original"; do
        
            for hrt_tool in "cliquesnv" "haplodmf" "rvhaplo"; do 

                # merge region standardised tsv files called "standard_output.tsv" into one file
                python -m benchmarking_scripts.merge_standard_output_files --results_dir results_simulations_3/${hrt_tool}/${virus}/100/regions/${ec_tool}/${sample}/
            
            done
        done
    done
done