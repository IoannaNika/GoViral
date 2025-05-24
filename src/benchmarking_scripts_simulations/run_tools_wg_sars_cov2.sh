#!/bin/bash

#### inputs, each tool requires a different environment to be loaded. 
tool_to_run=$1

function_produce_sam_bam_bed(){
    ref_seq=$1
    input=$2
    ec_tool=$3
    sample=$4
    cov=$5

    echo $ref_seq
    echo $input
    echo $ec_tool
    echo $cov

    output="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/simulated_data/${cov}/whole_genome/${ec_tool}/${sample}/sam_bam"
    # echo $output
    mkdir -p $output

    output_file="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/simulated_data/${cov}/whole_genome/${ec_tool}/${sample}/sam_bam/reads_${ec_tool}"
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
    results_dir=$4

    in_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/simulated_data/${coverage}/whole_genome/${ec_method}/${sample}/sam_bam/reads_${ec_tool}.sam"
    out_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/${results_dir}/cliquesnv/sarscov2/${coverage}/whole_genome/${ec_method}/${sample}/"

    java -jar ../CliqueSNV/clique-snv.jar -m snv-pacbio -tf 0.0 -t 1 -log -in $in_path -outDir $out_path
}


function_run_rvhaplo(){
    sample=$1
    ec_method=$2
    reference=$3
    $coverage=$4
    results_dir=$5

    in_path="/mnt/data/Simulations/sarscov2/simulated_data/${coverage}/whole_genome/${ec_method}/${sample}/sam_bam/reads_${ec_method}.sam"
    out_path="/mnt/src/${results_dir}/rvhaplo/sarscov2/${coverage}/whole_genome/${ec_method}/${sample}/"
    
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
    results_dir=$5
    
    in_path="/mnt/data/Simulations/sarscov2/simulated_data/${coverage}/whole_genome/${ec_method}/${sample}/sam_bam/reads_${ec_method}.sam"
    out_path="/mnt/src/${results_dir}/haplodmf/sarscov2/${coverage}/whole_genome/${ec_method}/${sample}/"
    
    mkdir -p /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src/${results_dir}/haplodmf/sarscov2/${coverage}/whole_genome/${ec_method}/${sample}

    cd /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/HaploDMF
    image_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/images/haplodmf.sif"
    mountdir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking"
    apptainer exec --bind $mountdir:/mnt $image_path ./haplodmf.sh -i $in_path -r $reference -o $out_path
    cd /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src
}

#### logic

results_dir="results_simulations_3"

ref_seq="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/LUMC/ref/nCoV-2019.reference.fasta"

for coverage in 100; do
      
    for n_haps in "n_1" "n_5" "n_10" "n_15" "n_20"; do

        for ec_tool in "hicanu" "original"; do

            if [[ "$ec_tool" == "hicanu" ]]; then   
                in_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/simulated_data/${coverage}/whole_genome/${ec_tool}/${n_haps}/output.canu.correctedReads.fasta"
            else 
                in_path="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/simulated_data/${coverage}/whole_genome/${ec_tool}/${n_haps}/output/reads.fasta"
            fi 
            
            function_produce_sam_bam_bed $ref_seq $in_path $ec_tool $n_haps  $coverage
   
            if [[ "$tool_to_run" == "cliquesnv" ]]; then
                function_run_cliquesnv $n_haps $ec_tool $coverage $results_dir
            fi

            if [[ "$tool_to_run" == "rvhaplo" ]]; then
                function_run_rvhaplo $sample $ec_tool $ref_seq $coverage $results_dir
            fi

            if [[ "$tool_to_run" == "haplodmf" ]]; then
                function_run_haplodmf $sample $ec_tool $ref_seq $coverage $results_dir
            fi
        done
    
    done

done
