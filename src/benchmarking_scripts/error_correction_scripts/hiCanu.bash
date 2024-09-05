
for sample in "03_50" "01_100" "02_100" "04_75" "05_90" "06_95" "07_98" "08_0" "09_0"
do 
    canu-2.2/bin/canu -correct -d error_corrected_data/hicanu/lumc/$sample/ -p output.canu useGrid=false maxMemory=90g genomeSize=29.8k saveReads=True stopOnLowCoverage=0 -pacbio data/lumc_data/natural_mixtures/$sample/trimmed_$sample.fastq

    gzip -d -c error_corrected_data/hicanu/lumc/$sample/output.canu.correctedReads.fasta.gz > error_corrected_data/hicanu/lumc/$sample/output.canu.correctedReads.fasta

done
