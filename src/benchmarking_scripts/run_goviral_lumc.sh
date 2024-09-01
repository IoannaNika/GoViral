
lumc_sample=$1

outdir="results/goViral/sample_${lumc_sample}"
input_fastq="../data/LUMC/whole_genome/original/trimmed_${lumc_sample}.fastq"
primers="../data/LUMC/primers/nCoV-2019.bed"
ref_seq="../data/LUMC/ref/nCoV-2019.reference.fasta"

python -m goViral.goViral_pipeline --directory ${outdir} --input_fastq ${input_fastq} --primers ${primers} --ref_seq ${ref_seq} --seed_limit 10
