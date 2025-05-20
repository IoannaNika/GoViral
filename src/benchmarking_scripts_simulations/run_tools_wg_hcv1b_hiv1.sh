#!/bin/bash

#### inputs, each tool requires a different environment to be loaded. 
tool_to_run=$1

function_produce_sam_bam_bed(){
    ref_seq=$1
    input=$2
    ec_tool=$3
    sample=$4
    cov=$5
    virus=$6

    echo $ref_seq
    echo $input
    echo $ec_tool
    echo $cov
    echo $virus

    output="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/${virus}/simulated_data/${coverage}/whole_genome/${ec_tool}/${sample}/sam_bam"
    # echo $output
    mkdir -p $output

    output_file="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/${virus}/simulated_data/${coverage}/whole_genome/${ec_tool}/${sample}/sam_bam/reads_${ec_tool}"
    # echo $output_file

    # align reads to reference using minimap2
    minimap2 -ax map-pb $ref_seq $input > $output_file.sam 
    samtools view -bS $output_file.sam > $output_file.bam
    samtools sort $output_file.bam -o $output_file.sorted.bam
    samtools index $output_file.sorted.bam
    bedtools bamtobed -i $output_file.sorted.bam > $output_file.bed
}

#### functions
function_run_cliquesnv(){
    sample=$1
    ec_method=$2
    coverage=$3
    virus=$4
    results_dir=$5

    in_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/${virus}/simulated_data/${coverage}/whole_genome/${ec_method}/${sample}/sam_bam/reads_${ec_tool}.sam"
    out_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/${results_dir}/cliquesnv/${virus}/${coverage}/whole_genome/${ec_method}/${sample}/"

    java -jar ../CliqueSNV/clique-snv.jar -m snv-pacbio -tf 0.0 -t 1 -log -in $in_path -outDir $out_path
}


function_run_rvhaplo(){
    sample=$1
    ec_method=$2
    reference=$3
    coverage=$4
    virus=$5
    results_dir=$6

    in_path="/mnt/data/Simulations/${virus}/simulated_data/${coverage}/whole_genome/${ec_method}/${sample}/sam_bam/reads_${ec_tool}.sam"
    out_path="/mnt/src/${results_dir}/rvhaplo/${virus}/${coverage}/whole_genome/${ec_method}/${sample}/"
    
    cd ../RVHaplo
    image_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/images/rvhaplo_image.sif"
    mountdir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking"
    apptainer exec --bind $mountdir:/mnt $image_path ./rvhaplo.sh -i $in_path -r $reference -o $out_path --error_rate 0.01 --abundance 0
    cd ../src
}

function_run_haplodmf(){
    sample=$1
    ec_method=$2
    reference=$3
    coverage=$4
    virus=$5
    results_dir=$6
    
    in_path="/mnt/data/Simulations/${virus}/simulated_data/${coverage}/whole_genome/${ec_method}/${sample}/sam_bam/reads_${ec_tool}.sam"
    out_path="/mnt/src/${results_dir}/haplodmf/${virus}/${coverage}/whole_genome/${ec_method}/${sample}/"
    
    cd /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/HaploDMF
    image_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/images/haplodmf.sif"
    mountdir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking"
    apptainer exec --bind $mountdir:/mnt $image_path ./haplodmf.sh -i $in_path -r $reference -o $out_path
    cd /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src
}

#### logic

results_dir="results_simulations_2"

ref_seq_hcv1b="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hcv1b95/reference.fasta"
ref_seq_hiv1="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hiv1/reference.fasta"

for coverage in 100; do

    for virus in "hcv1b95" "hiv1"; do 
        
        if [[ "$virus" == "hcv1b95" ]]; then
            ref_seq="${ref_seq_hcv1b}"
        else
            ref_seq="${ref_seq_hiv1}"
        fi
        
        for sample in "ab_3_97" "ab_30_70" "ab_50_50" "ab_70_30" "ab_97_3"; do
                        
            for ec_tool in "hicanu" "original"; do

                echo "${virus} ${sample} ${ec_tool}"
                
                in_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/${virus}/simulated_data/${coverage}/whole_genome/${ec_tool}/${sample}/output.canu.correctedReads.fasta"
                
                function_produce_sam_bam_bed $ref_seq $in_path $ec_tool $sample $coverage $virus 

                if [[ "$tool_to_run" == "cliquesnv" ]]; then
                    function_run_cliquesnv $sample $ec_tool $coverage $virus $results_dir
                fi

                if [[ "$tool_to_run" == "rvhaplo" ]]; then
                    function_run_rvhaplo $sample $ec_tool $ref_seq $coverage $virus $results_dir
                fi

                if [[ "$tool_to_run" == "haplodmf" ]]; then
                    function_run_haplodmf $sample $ec_tool $ref_seq $coverage $virus $results_dir
                fi
            done
    
        done
    done
done
