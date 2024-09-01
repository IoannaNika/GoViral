# GoViral pipeline

GoViral is a viral haplotype reconstruction pipeline for long amplicon PacBio-Hifi sequencing reads. In this repository we provide a simulator for long amplicon sequencing PacBio-Hifi and ONT reads. The simulator is build using the pbsim3 and the ccs tools.

## How to use

#### GoViral pipeline
```
python -m goViral.goViral_pipeline --directory <OUTPUT_DIRECTORY> --input_fastq <INPUT_FASTQ_FILE> --primers <.BED FILE WITH PRIMERS AND THEIR POSITIONS> --ref_seq <REFERENCE_SEQUENCE>
```

Other options include: 

`--coverage_limit` How many reads to consider in each subsample, default is 100.

`--seed_limit` Seed limit for subsampling. How many subsamples to consider.

`--follow_reccomendation` sets the seet limit based on the coverage of the dataset provided.

Increasing the coverage limit or the seed limit will increase the run-time. The run-time can increase also when the follow_reccomendation option is set.

`--ab_threshold` Abundance threshold for filtering out low abundandant reconstructed haplotypes.

#### Long amplicon read simulation

1. Benchmark dataset creation

```
python -m long_amplicon_read_simulation.run_benchmark_creation --input_fasta <INPUT_FASTA_SEQUENCES> --mixture <JSON FILE_WITH_SEQUENCE_IDS_AND_THEIR_DESIRED_FRACTION_> --coverage <NUMBER_SPECIFYING_COVERAGE> --primers_file <.BED FILE WITH PRIMERS AND THEIR POSITIONS> --strategy <pacbio-hifi or ONT> --outdir <OUTPUT_DIRECTORY>
```

2. Training dataset creation

```
python -m long_amplicon_read_simulation.run_training_set_creation --input_fasta <INPUT_FASTA_SEQUENCES> --primers_file <.BED FILE WITH PRIMERS AND THEIR POSITIONS> --strategy <pacbio-hifi or ONT> --outdir <OUTPUT_DIRECTORY> --n <NUMBER OF EXAMPLES TO CREATE>
```

## Example

For the file formats consult the files provided in this example.




