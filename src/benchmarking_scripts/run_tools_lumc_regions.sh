#!/bin/bash

#### inputs, each tool requires a different environment to be loaded. 
tool_to_run=$1

#### functions
function_run_cliquesnv(){
    sample=$1
    start=$2
    end_region=$3
    ec_method=$4

    in_path="../data/LUMC/per_region/${ec_method}/${sample}/${sample}_${start}_${end_region}.sam"
    out_path="results_debug/cliquesnv/${ec_method}/per_region/${sample}/${start}_${end_region}/"
    
    mkdir -p $out_path

    java -jar ../CliqueSNV/clique-snv.jar -m snv-pacbio -log -in $in_path -outDir $out_path -sp ${start} -ep ${end_region} -t 1 -tf 0
}

function_run_haplodmf(){
    sample=$1
    start=$2
    end_region=$3
    ec_method=$4
    reference=$5
    
    in_path="/mnt/data/LUMC/per_region/${ec_method}/${sample}/${sample}_${start}_${end_region}.sam"
    out_path="/mnt/src/results_debug/haplodmf/${ec_method}/per_region/${sample}/${start}_${end_region}/"
    
    mkdir -p "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/results_debug/haplodmf/${ec_method}/per_region/${sample}/${start}_${end_region}/"

    cd ../HaploDMF/
    image_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/images/haplodmf.sif"
    mountdir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking"
    apptainer exec --bind $mountdir:/mnt $image_path ./haplodmf.sh -i $in_path -r $reference -o $out_path  -sp ${start} -ep ${end_region}
    cd ../src
}

function_run_rvhaplo(){
    sample=$1
    start=$2
    end_region=$3
    ec_method=$4
    reference=$5

    in_path="/mnt/data/LUMC/per_region/${ec_method}/${sample}/${sample}_${start}_${end_region}.sam"
    out_path="/mnt/src/results_debug/rvhaplo/${ec_method}/per_region/${sample}/${start}_${end_region}/"
    
    cd ../RVHaplo/
    image_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/images/rvhaplo_image.sif"
    mountdir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking"
    apptainer exec --bind $mountdir:/mnt $image_path ./rvhaplo.sh -i $in_path -r $reference -o $out_path --error_rate 0.01 -sp ${start} -ep ${end_region} -a 0
    cd ../src
}


#### logic

reference="../data/LUMC/ref/nCoV-2019.reference.fasta"

genomic_regions="(54,1183) (1128,2244) (2179,3235) (3166,4240) (4189,5337) (5286,6358) (6307,7379) (7328,8363) (8282,9378) (9327,10429) (10370,11447) (11394,12538) (12473,13599) (13532,14619) (14568,15713) (15634,16698) (16647,17732) (17649,18684) (18618,19655) (19604,20676) (20581,21620) (21562,22590) (22537,23609) (23544,24714) (24658,25768) (25712,26835) (26766,27872) (27808,28985) (28699,29768) (29768,29790)"


for sample in "01_100" "02_100" "03_50" "04_75" "05_90" "06_95" "07_98" "08_0" "09_0"
do 
    # Iterate over each genomic region
    for region in $genomic_regions; do
        # Remove the parentheses and split the region into start and end
        region=$(echo $region | tr -d '()')
    
        # Extract the first and last number using a comma as the delimiter
        start=$(echo $region | cut -d',' -f1)
        end_region=$(echo $region | cut -d',' -f2)

        for ec_tool in "canu" "hifiasm" "lorma" "original"; do

            echo "ec_tool ${ec_tool}"
            echo "region: ${start} ${end_region}"
        
            if [[ "$tool_to_run" == "cliquesnv" ]]; then
                function_run_cliquesnv $sample $start $end_region $ec_tool
            fi
            if [[ "$tool_to_run" == "haplodmf" ]]; then
                function_run_haplodmf $sample $start $end_region $ec_tool $reference
            fi
            if [[ "$tool_to_run" == "rvhaplo" ]]; then
                function_run_rvhaplo $sample $start $end_region $ec_tool $reference
            fi
            
        done

    done
done