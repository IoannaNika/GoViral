#!/bin/bash

#### inputs, each tool requires a different environment to be loaded. 
tool_to_run=$1

#### functions
function_run_cliquesnv(){
    sample=$1
    ec_method=$2

    in_path="../data/LUMC/whole_genome/${ec_method}/${sample}/${sample}.sam"
    out_path="results/cliquesnv/${ec_method}/whole_genome/${sample}/"

    java -jar CliqueSNV/clique-snv.jar -m snv-pacbio -log -in $in_path -outDir $out_path
}

function_run_haplodmf(){
    sample=$1
    ec_method=$2
    reference=$3
    
    in_path="../data/LUMC/whole_genome/${ec_method}/${sample}/.sam"
    out_path="results/whole_genome/haplodmf/${ec_method}/per_region/${sample}/"

    /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/HaploDMF/haplodmf.sh -i $in_path -r $reference -o $out_path
}

function_run_rvhaplo(){
    sample=$1
    ec_method=$2
    reference=$3

    in_path="../data/LUMC/whole_genome/${ec_method}/${sample}/${sample}.sam"
    out_path="results/whole_genome/rvhaplo/${ec_method}/per_region/${sample}/"

    /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/RVHaplo/rvhaplo.sh -i $in_path -r $reference -o $out_path --error_rate 0.01
}

#### logic

reference="../data/LUMC/ref/nCoV-2019.reference.fasta"

for sample in "03_50" "01_100" "02_100" "04_75" "05_90" "06_95" "07_98" "08_0" "09_0"; do 
  
    for tool in "original" "canu" "hifiasm" "lorma"; do
        if tool_to_run="cliquesnv"; then
            function_run_cliquesnv $sample $tool
        fi
        if tool_to_run="haplodmf"; then
            function_run_haplodmf $sample $tool $reference
        fi
        if tool_to_run="rvhaplo"; then
            function_run_rvhaplo $sample $tool $reference
        fi
        
    done   

done