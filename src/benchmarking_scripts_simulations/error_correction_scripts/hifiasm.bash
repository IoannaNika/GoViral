
for sample in "n_1" "n_5" "n_10" "n_15" "n_20"
do 
    inputfasta="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/simulated_data/100/whole_genome/original/${sample}/output/reads.fasta"
    outputdir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/sarscov2/simulated_data/100/whole_genome/hifiasm/${sample}"

    echo "inputfasta=$inputfasta"
    echo "outputdir=$outputdir"

    /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/hifiasm/hifiasm --write-ec -o ${outputdir}/assembly.asm ${inputfasta}

done

for virus in "hiv1" "hcv1b95"; do 
    
    for ab_sample in "ab_3_97" "ab_30_70" "ab_50_50" "ab_70_30" "ab_97_3"; do 
        
        inputfasta="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/${virus}/simulated_data/100/whole_genome/original/${ab_sample}/output/reads.fasta"
        outputdir="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Benchmarking/data/Simulations/${virus}/simulated_data/100/whole_genome/hifiasm/${ab_sample}"

        /tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/hifiasm/hifiasm --write-ec -o ${outputdir}/assembly.asm ${inputfasta}

    done

done