accelerator="gpu"
out_dir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/src"
data_dir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations"
xsim_100=("ab_3_97" "ab_30_70" "ab_50_50" "ab_70_30" "ab_97_3")
xsim_1000=("ab_3_997" "ab_30_970" "ab_300_700" "ab_500_500" "ab_700_300" "ab_970_30" "ab_3_997")

genomic_regions_hiv1=$(awk '$4 ~ /_LEFT/ {print $3}' /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hiv1/primers.bed)
genomic_regions_hcv1b=$(awk '$4 ~ /_LEFT/ {print $3}' /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/hcv1b95/primers.bed)

for coverage in 100 1000; do

    for virus in "hiv1" "hcv1b95"; do

        if [ "$virus" == "hiv1" ]; then
            genomic_regions=("${genomic_regions_hiv1[@]}")
        elif [ "$virus" == "hcv1b95" ]; then
            genomic_regions=("${genomic_regions_hcv1b[@]}")
        fi
            
        if [ $coverage -eq 100 ]; then

            for sim_ab in "${xsim_100[@]}"; do
        
                outdir="${out_dir}/results_simulations/goViral/${virus}/original/${coverage}/${sim_ab}"
                datadir="${data_dir}/${virus}/simulated_data/${coverage}/original/${sim_ab}/output/reads.fasta"
                ref_seq="${data_dir}/${virus}/reference.fasta"
                primers="${data_dir}/${virus}/primers.bed"
                
                mkdir -p $outdir
                
                python  -m goViral.goViral_pipeline --directory $outdir --input_fastq $datadir --ab_threshold 0 --follow_reccomendations --ref_seq $ref_seq --primers $primers

            done

        elif [ $coverage -eq 1000 ]; then

            for sim_ab in "${xsim_1000[@]}"; do
                    
                outdir="${out_dir}/results_simulations/goViral/${virus}/original/${coverage}/${sim_ab}/"
                datadir="${data_dir}/${virus}/simulated_data/${coverage}/original/${sim_ab}/output/reads.fasta"
                ref_seq="${data_dir}/${virus}/reference.fasta"
                primers="${data_dir}/${virus}/primers.bed"
                mkdir -p $outdir
                
                python  -m goViral.goViral_pipeline --directory $outdir --input_fastq $datadir --ab_threshold 0 --follow_reccomendations --ref_seq $ref_seq --primers $primers

            done

        fi

    done
    
done

