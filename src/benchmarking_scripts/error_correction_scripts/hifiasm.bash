
for sample in "01_100" "02_100" "03_50" "04_75" "05_90" "06_95" "07_98" "08_0" "09_0"
do 
    hifiasm/hifiasm --write-ec -o error_corrected_data/hifiasm/lumc/$sample/assembly.asm data/lumc_data/natural_mixtures/$sample/trimmed_$sample.fastq

done