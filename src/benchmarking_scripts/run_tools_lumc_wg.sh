#!/bin/bash

#### inputs, each tool requires a different environment to be loaded. 
tool_to_run=$1

#### functions
function_run_cliquesnv(){
    sample=$1
    ec_method=$2

    in_path="../data/LUMC/whole_genome/${ec_method}/${sample}/${sample}.sam"
    out_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/results_debug/cliquesnv/${ec_method}/whole_genome/${sample}/"

    java -jar ../CliqueSNV/clique-snv.jar -m snv-pacbio -log -in $in_path -outDir $out_path -t 1 -tf 0
}

function_run_haplodmf(){
    sample=$1
    ec_method=$2
    reference=$3
    
    in_path="/mnt/data/LUMC/whole_genome/${ec_method}/${sample}/${sample}.sam"
    out_path="/mnt/src/results_debug/haplodmf/${ec_method}/whole_genome/${sample}/"

    mkdir -p "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/results_debug/haplodmf/${ec_method}/whole_genome/${sample}/"
    
    cd /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/HaploDMF
    image_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/images/haplodmf.sif"
    mountdir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking"
    apptainer exec --bind $mountdir:/mnt $image_path ./haplodmf.sh -i $in_path -r $reference -o $out_path
    cd /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src
}

function_run_rvhaplo(){
    sample=$1
    ec_method=$2
    reference=$3

    in_path="/mnt/data/LUMC/whole_genome/${ec_method}/${sample}/${sample}.sam"
    out_path="/mnt/src/results_debug/rvhaplo/${ec_method}/whole_genome/${sample}/"
    
    cd ../RVHaplo
    image_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/images/rvhaplo_image.sif"
    mountdir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking"
    apptainer exec --bind $mountdir:/mnt $image_path ./rvhaplo.sh -i $in_path -r $reference -o $out_path --error_rate 0.01 -a 0
    cd ../src
}

function_run_viloca(){
    sample=$1
    ec_method=$2
    reference="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/LUMC/ref/nCoV-2019.reference.fasta"

    in_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/LUMC/whole_genome/${ec_method}/${sample}/${sample}.sorted.bam"
    out_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/results_debug/viloca/${ec_method}/whole_genome/${sample}/"
    file="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/LUMC/primers/nCoV-2019.bed"
    
    mkdir -p $out_path

    # samtools sort $in_path -o $in_path_sorted

    cd $out_path
    viloca run -f $reference -b $in_path -z $zfile 
    cd /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src

}

#### logic

reference="../data/LUMC/ref/nCoV-2019.reference.fasta"

for sample in "01_100" "02_100" "03_50" "04_75" "05_90" "06_95" "07_98" "08_0" "09_0"; do 

    for ec_tool in "hifiasm" "lorma" "original" "canu"; do
        echo $sample
        echo $ec_tool
   
        if [[ "$tool_to_run" == "cliquesnv" ]]; then
            function_run_cliquesnv $sample $ec_tool
        fi
        if [[ "$tool_to_run" == "haplodmf" ]]; then
            function_run_haplodmf $sample $ec_tool $reference
        fi
        if [[ "$tool_to_run" == "rvhaplo" ]]; then
            function_run_rvhaplo $sample $ec_tool $reference
        fi
        if [[ "$tool_to_run" == "viloca" ]]; then
            function_run_viloca $sample $ec_tool
        fi
    
    done

done