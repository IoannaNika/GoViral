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

    in_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/simulated_data/${coverage}/regions/${ec_method}/${sample}/sam_bam/region_${start}_${end_region}.sam"
    out_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/${results_dir}/cliquesnv/sarscov2/${coverage}/regions/${ec_method}/${sample}/${start}_${end_region}"
    
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

    in_path="/mnt/data/Simulations/sarscov2/simulated_data/${coverage}/regions/${ec_method}/${sample}/sam_bam/region_${start}_${end_region}.sam"
    out_path="/mnt/src/${results_dir}/rvhaplo/sarscov2/${coverage}/regions/${ec_method}/${sample}/${start}_${end_region}"
    
    in_path_actual="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/simulated_data/${coverage}/regions/${ec_method}/${sample}/sam_bam/region_${start}_${end_region}.sam"
    
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
    
    in_path="/mnt/data/Simulations/sarscov2/simulated_data/${coverage}/regions/${ec_method}/${sample}/sam_bam/region_${start}_${end_region}.sam"
    out_path="/mnt/src/${results_dir}/haplodmf/sarscov2/${coverage}/regions/${ec_method}/${sample}/${start}_${end_region}"
    
    in_path_actual="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/simulated_data/${coverage}/regions/${ec_method}/${sample}/sam_bam/region_${start}_${end_region}.sam"

    if [[ ! -e "$in_path_actual" ]]; then
        echo "File does not exist: $in_path_actual"
        return 1  
    fi

    mkdir -p /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/${results_dir}/haplodmf/sarscov2/${coverage}/regions/${ec_method}/${sample}/${start}_${end_region}

    cd /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/HaploDMF
    image_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/images/haplodmf.sif"
    mountdir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking"
    apptainer exec --bind $mountdir:/mnt $image_path ./haplodmf.sh -i $in_path -r $reference -o $out_path -sp ${start} -ep ${end_region}
    cd /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src
}

#### logic

results_dir="results_simulations_3"

ref_seq="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/LUMC/ref/nCoV-2019.reference.fasta"

genomic_regions="(54,1183) (1128,2244) (2179,3235) (3166,4240) (4189,5337) (5286,6358) (6307,7379) (7328,8363) (8282,9378) (9327,10429) (10370,11447) (11394,12538) (12473,13599) (13532,14619) (14568,15713) (15634,16698) (16647,17732) (17649,18684) (18618,19655) (19604,20676) (20581,21620) (21562,22590) (22537,23609) (23544,24714) (24658,25768) (25712,26835) (26766,27872) (27808,28985) (28699,29768) (29768,29790)"

for coverage in 100; do
      
    for n_haps in "n_1" "n_5" "n_10" "n_15" "n_20"; do

        for ec_tool in "hicanu" "original"; do

            for region in $genomic_regions; do
                # Remove the parentheses and split the region into start and end
                region=$(echo $region | tr -d '()')
            
                # Extract the first and last number using a comma as the delimiter
                start=$(echo $region | cut -d',' -f1)
                end_region=$(echo $region | cut -d',' -f2)

               
                if [[ "$tool_to_run" == "cliquesnv" ]]; then
                    function_run_cliquesnv $n_haps $ec_tool $coverage $results_dir $start $end_region
                fi

                if [[ "$tool_to_run" == "rvhaplo" ]]; then
                    function_run_rvhaplo $n_haps $ec_tool $ref_seq $coverage $results_dir $start $end_region
                fi

                if [[ "$tool_to_run" == "haplodmf" ]]; then
                    function_run_haplodmf $n_haps $ec_tool $ref_seq $coverage $results_dir $start $end_region
                fi
                
            done

        done
    
    done

done
