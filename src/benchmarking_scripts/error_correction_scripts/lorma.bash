
for sample in "01_100" "02_100" "04_75" "05_90" "06_95" "07_98" "08_0" "09_0" # 03_50
do 
    LoRMA -reads data/lumc_data/natural_mixtures/$sample/trimmed_$sample.fastq -output error_corrected_data/lorma/lumc/$sample/corrected.fasta -discarded error_corrected_data/lorma/lumc/$sample/discarded.fasta
    
done