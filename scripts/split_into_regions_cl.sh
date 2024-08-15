ref_seq="data/LUMC/ref/nCoV-2019.reference.fasta"

### Lorma
echo "LORMA"
in_path_lorma="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/error_corrected_data/lorma/lumc/03_50/corrected"
out_path_lorma="data/LUMC/error_corrected_per_region/lorma/03_50"

python scripts/split_into_regions_cl.py --in_path $in_path_lorma --out_path $out_path_lorma --ref_seq $ref_seq

### Canu
echo "CANU"
in_path_canu="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/error_corrected_data/hicanu/lumc/03_50/output.canu.correctedReads"
out_path_canu="data/LUMC/error_corrected_per_region/canu/03_50"

python scripts/split_into_regions_cl.py --in_path $in_path_canu --out_path $out_path_canu --ref_seq $ref_seq


### Hifiasm
echo "HIFIASM"
in_path_hifiasm="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/error_corrected_data/hifiasm/lumc/03_50/assembly.asm.ec"
out_path_hifiasm="data/LUMC/error_corrected_per_region/hifiasm/03_50"

python scripts/split_into_regions_cl.py --in_path $in_path_hifiasm --out_path $out_path_hifiasm --ref_seq $ref_seq
