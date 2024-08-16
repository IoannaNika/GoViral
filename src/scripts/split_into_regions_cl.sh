function_produce_sam_bam_bed(){
    ref_seq=$1
    input=$2
    ec_tool=$3
    sample=$4

    # echo $ref_seq
    # echo $input
    # echo $ec_tool

    output="../data/LUMC/whole_genome/${ec_tool}/${sample}/"
    # echo $output
    mkdir $output

    output_file="../data/LUMC/whole_genome/${ec_tool}/${sample}/${sample}"
    # echo $output_file

    # align reads to reference using minimap2
    minimap2 -ax map-pb $ref_seq $input.fastq > $output_file.sam 
    samtools view -bS $output_file.sam > $output_file.bam
    samtools sort $output_file.bam -o $output_file.sorted.bam
    samtools index $output_file.sorted.bam
    bedtools bamtobed -i $output_file.sorted.bam > $output_file.bed
}

ref_seq="../data/LUMC/ref/nCoV-2019.reference.fasta"

for sample in "03_50" # "01_100" "02_100" "04_75" "05_90" "06_95" "07_98" "08_0" "09_0"
do
    ### Lorma
    tool="lorma"
    echo $tool
    in_path_lorma="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/error_corrected_data/lorma/lumc/$sample/corrected"
    bed_file_lorma="../data/LUMC/whole_genome/${tool}/${sample}/${sample}"
    out_path_lorma="../data/LUMC/per_region/${tool}/${sample}/${sample}"

    mkdir "../data/LUMC/per_region/${tool}/${sample}/"    

    function_produce_sam_bam_bed $ref_seq $in_path_lorma $tool $sample

    python scripts/split_into_regions_cl.py --input_path $in_path_lorma --bed $bed_file_lorma --out_path $out_path_lorma --ref_seq $ref_seq

    ### Canu
    tool="canu"
    echo $tool
    in_path_canu="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/error_corrected_data/hicanu/lumc/$sample/output.canu.correctedReads"
    bed_file_canu="../data/LUMC/whole_genome/${tool}/${sample}/${sample}"
    out_path_canu="../data/LUMC/per_region/${tool}/${sample}/${sample}"

    mkdir "../data/LUMC/per_region/${tool}/${sample}/"

    function_produce_sam_bam_bed $ref_seq $in_path_canu $tool $sample

    python scripts/split_into_regions_cl.py --input_path $in_path_canu --bed $bed_file_canu --out_path $out_path_canu --ref_seq $ref_seq


    ### Hifiasm
    tool="hifiasm"
    echo $tool
    in_path_hifiasm="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/error_corrected_data/hifiasm/lumc/$sample/assembly.asm.ec"
    bed_file_hifiasm="../data/LUMC/whole_genome/${tool}/${sample}/${sample}"
    out_path_hifiasm="../data/LUMC/per_region/${tool}/${sample}/${sample}"

    mkdir "../data/LUMC/per_region/${tool}/${sample}/"

    function_produce_sam_bam_bed $ref_seq $in_path_hifiasm $tool $sample

    python scripts/split_into_regions_cl.py --input_path $in_path_hifiasm --bed $bed_file_hifiasm --out_path $out_path_hifiasm --ref_seq $ref_seq

    
    ### Original
    tool="original"
    echo $tool
    ## Non-error corrected
    in_path="../data/LUMC/whole_genome/original/trimmed_$sample"
    bed_file="../data/LUMC/whole_genome/${tool}/${sample}/${sample}"
    out_path="../data/LUMC/per_region/${tool}/${sample}/${sample}"

    mkdir "../data/LUMC/per_region/${tool}/${sample}/"

    function_produce_sam_bam_bed $ref_seq $in_path $tool $sample
    
    python scripts/split_into_regions_cl.py --input_path $in_path --bed $bed_file --out_path $out_path --ref_seq $ref_seq

done