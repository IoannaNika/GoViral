#!/bin/bash

#### inputs, each tool requires a different environment to be loaded. 
tool_to_run=$1

#### functions
function_run_cliquesnv(){
    sample=$1
    ec_method=$2
    coverage=$3
    results_dir=$4
    start=$5
    end_region=$6
    virus=$7

    in_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/${virus}/simulated_data/${coverage}/regions/${ec_method}/${sample}/sam_bam/region_${start}_${end_region}.sam"
    out_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/${results_dir}/cliquesnv/${virus}/${coverage}/regions/${ec_method}/${sample}/${start}_${end_region}"
    
    if [[ ! -e "$in_path" ]]; then
        echo "File does not exist: $in_path"
        return 1  
    fi

    java -jar ../CliqueSNV/clique-snv.jar -m snv-pacbio -tf 0.0 -t 1 -log -in $in_path -outDir $out_path -sp ${start} -ep ${end_region}
}


function_run_rvhaplo(){
    sample=$1
    ec_method=$2
    reference=$3
    coverage=$4
    results_dir=$5
    start=$6
    end_region=$7
    virus=$8

    in_path="/mnt/data/Simulations/${virus}/simulated_data/${coverage}/regions/${ec_method}/${sample}/sam_bam/region_${start}_${end_region}.sam"
    out_path="/mnt/src/${results_dir}/rvhaplo/${virus}/${coverage}/regions/${ec_method}/${sample}/${start}_${end_region}"
    
    in_path_actual="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/${virus}/simulated_data/${coverage}/regions/${ec_method}/${sample}/sam_bam/region_${start}_${end_region}.sam"
    
    if [[ ! -e "$in_path_actual" ]]; then
        echo "File does not exist: $in_path_actual"
        return 1  
    fi

    cd ../RVHaplo
    image_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/images/rvhaplo_image.sif"
    mountdir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking"
    apptainer exec --bind $mountdir:/mnt $image_path ./rvhaplo.sh -i $in_path -r $reference -o $out_path --error_rate 0.01 --abundance 0 -sp ${start} -ep ${end_region}
    cd ../src
}

function_run_haplodmf(){
    sample=$1
    ec_method=$2
    reference=$3
    coverage=$4
    results_dir=$5
    start=$6
    end_region=$7
    virus=$8
    
    in_path="/mnt/data/Simulations/${virus}/simulated_data/${coverage}/regions/${ec_method}/${sample}/sam_bam/region_${start}_${end_region}.sam"
    out_path="/mnt/src/${results_dir}/haplodmf/${virus}/${coverage}/regions/${ec_method}/${sample}/${start}_${end_region}"
    
    in_path_actual="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/${virus}/simulated_data/${coverage}/regions/${ec_method}/${sample}/sam_bam/region_${start}_${end_region}.sam"

    if [[ ! -e "$in_path_actual" ]]; then
        echo "File does not exist: $in_path_actual"
        return 1  
    fi

    mkdir -p /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/${results_dir}/haplodmf/${virus}/${coverage}/regions/${ec_method}/${sample}/${start}_${end_region}

    cd /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/HaploDMF
    image_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/images/haplodmf.sif"
    mountdir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking"
    apptainer exec --bind $mountdir:/mnt $image_path ./haplodmf.sh -i $in_path -r $reference -o $out_path -sp ${start} -ep ${end_region}
    cd /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src
}

#### logic

results_dir="results_simulations_3"

ref_seq_hcv1b="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hcv1b95/reference.fasta"
ref_seq_hiv1="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hiv1/reference.fasta"

genomic_regions_hiv1=$(awk '/LEFT/ { left_end = $3 } /RIGHT/ { right_start = $2; print "(" left_end "," right_start ")" }' /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hiv1/primers.bed | paste -sd' ' -)
genomic_regions_hcv1b=$(awk '/LEFT/ { left_end = $3 } /RIGHT/ { right_start = $2; print "(" left_end "," right_start ")" }' /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hcv1b95/primers.bed | paste -sd' ' -)

for coverage in 100; do
      
    for sample in "ab_3_97" "ab_30_70" "ab_50_50" "ab_70_30" "ab_97_3"; do

        for ec_tool in "hicanu" "original"; do
            
            for virus in "hcv1b95" "hiv1"; do 

                if [[ "$virus" == "hcv1b95" ]]; then
                    ref_seq="${ref_seq_hcv1b}"
                else
                    ref_seq="${ref_seq_hiv1}"
                fi

                if [ "$virus" == "hiv1" ]; then
                    genomic_regions=("${genomic_regions_hiv1[@]}")
                elif [ "$virus" == "hcv1b95" ]; then
                    genomic_regions=("${genomic_regions_hcv1b[@]}")
                fi

                for region in $genomic_regions; do
                    # Remove the parentheses and split the region into start and end
                    region=$(echo $region | tr -d '()')
                
                    # Extract the first and last number using a comma as the delimiter
                    start=$(echo $region | cut -d',' -f1)
                    end_region=$(echo $region | cut -d',' -f2)

                
                    if [[ "$tool_to_run" == "cliquesnv" ]]; then
                        function_run_cliquesnv $sample $ec_tool $coverage $results_dir $start $end_region $virus
                    fi

                    if [[ "$tool_to_run" == "rvhaplo" ]]; then
                        function_run_rvhaplo $sample $ec_tool $ref_seq $coverage $results_dir $start $end_region $virus
                    fi

                    if [[ "$tool_to_run" == "haplodmf" ]]; then
                        function_run_haplodmf $sample $ec_tool $ref_seq $coverage $results_dir $start $end_region $virus
                    fi
                    
                done

            done

        done
    
    done

done
