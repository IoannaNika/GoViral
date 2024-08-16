genomic_regions="(54,1183),(1128,2244),(2179,3235),(3166,4240),(4189,5337),
                    (5286,6358),(6307,7379),(7328,8363),(8282,9378),(9327,10429),
                    (10370,11447),(11394,12538),(12473,13599),(13532,14619),
                        (14568,15713),(15634,16698),(16647,17732),(17649,18684),
                    (18618,19655),(19604,20676),(20581,21620),(21562,22590),
                        (22537,23609),(23544,24714),(24658,25768),(25712,26835),
                        (26766,27872),(27808,28985),(28699,29768),(29768,29790)"



for hrt_tool in "cliquesnv" "rvhaplo" "haplodmf"; do
        
    for ec_tool in "lorma" "canu" "hifiasm" "original"; do
        
        for sample in "03_50"  "01_100" "02_100" "04_75" "05_90" "06_95" "07_98" "08_0" "09_0"; do 

            for region in $genomic_regions; do
                    
                # remove the parentheses and split the region into start and end
                region=$(echo $region | tr -d '()')
                start=$(echo $region | cut -d ',' -f 1)
                end_region=$(echo $region | cut -d ',' -f 2)
                region_final="${start}_${end_region}"

                python scripts/standardise_output_tools_regions.py --results_dir results/ --hrt $hrt_tool --ec $ec_tool --region $region_final
            
            done
        
        # merge region standardised tsv files called "standard_output.tsv" into one file
        python scripts/merge_standard_output.py --results_dir results/${hrt_tool}/${ec_tool}/per_region/
  
        done
    
    done

done

