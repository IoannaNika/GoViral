reference="../data/LUMC/ref/nCoV-2019.reference.fasta"
results_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/results"
output="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/results/results.tsv"

for sample in "01_100" "02_100" "03_50" "04_75" "05_90" "06_95" "07_98" "08_0" "09_0"; do   

    reads="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data/natural_mixtures/${sample}/reads.tsv"

    for ec_tool in "original" "canu" "hifiasm" "lorma"; do

        for hr_tool in "haplodmf" "cliquesnv" "rvhaplo"; do

            for mode in "whole_genome" "per_region"; do

                sample_name="${sample}-${ec_tool}-${hr_tool}-${mode}"
                
                echo ${sample_name}

                if [[ "${mode}" == "whole_genome" ]]; then
                    input="${results_path}/${hr_tool}/${ec_tool}/${mode}/${sample}/standard_output.tsv"
                else
                    input="${results_path}/${hr_tool}/${ec_tool}/${mode}/${sample}/merged_standard_output.tsv"
                fi

                echo $input

                # if the input file does not exist, skip
                if [[ ! -f $input ]]; then
                    continue
                fi
                
                echo "file exists"

                python -m benchmarking_scripts.evaluate_hrt_output_lumc --input $input --ref_seq $reference --output $output --sample_name $sample_name --reads $reads
            
            done
        done
    done
done

python scripts/summarise_evaluation_across_samples.py --results_file results/results.tsv --output_file results/summary.tsv --split_per_region_and_wg
