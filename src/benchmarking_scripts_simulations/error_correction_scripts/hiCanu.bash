
for sample in "n_1" "n_5" "n_10" "n_15" "n_20"
do 
    canu -correct -d data/Simulations/sarscov2/simulated_data/whole_genome/original/${sample}/output/reads.fasta
    gzip -d -c data/Simulations/sarscov2/simulated_data/whole_genome/hicanu/${sample}/output.canu.correctedReads.fasta.gz > /data/Simulations/sarscov2/simulated_data/whole_genome/hicanu/${sample}/output.canu.correctedReads.fasta
done
