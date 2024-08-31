
import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO
import subprocess

def simulate_ONT_HQ_reads(directory, identifier, cores): 
    os.system("pbsim --strategy templ --method errhmm --errhmm data/error_models/ERRHMM-ONT-HQ.model --template {}/{}.template --prefix {}/{}".format(directory, identifier, directory, identifier))
    return

def simulate_hifi_reads(directory, identifier, cores):
    os.system("pbsim --strategy templ --method errhmm --errhmm data/error_models/ERRHMM-SEQUEL.model  --template {}/{}.template --pass-num 10 --prefix {}/{}".format(directory, identifier, directory, identifier))
    os.system("samtools view -bS {}/{}.sam > {}/{}.bam".format(directory, identifier, directory, identifier))
    os.system("ccs {}/{}.bam -j {} {}/{}.fastq".format(directory, identifier, cores, directory, identifier))
    return

def parse_R_output(output):
    # parse output from R script
    output = output.split("\n")
    start = output[1].split(" ")[1]
    end = output[4].split(" ")[1]
    return start, end

def get_amplicon_positions(Fprob, Rprob, seq_path):
    # call r script to get amplicon positions
    results = subprocess.run(["Rscript match_primers.R " + Fprob + " " + Rprob + " " + seq_path], shell=True, capture_output=True, text=True)
    results = parse_R_output(results.stdout)
    start = int(results[0]) + len(Fprob)
    end = int(results[1]) - len(Rprob)
    return start, end


def simulate_ONT_HQ_reads(directory, identifier, cores): 
    os.system("pbsim --strategy templ --method errhmm --errhmm data/error_models/ERRHMM-ONT-HQ.model --template {}/{}.template --prefix {}/{}".format(directory, identifier, directory, identifier))
    return

def simulate_hifi_reads(directory, identifier, cores):
    os.system("pbsim --strategy templ --method errhmm --errhmm data/error_models/ERRHMM-SEQUEL.model  --template {}/{}.template --pass-num 10 --prefix {}/{}".format(directory, identifier, directory, identifier))
    os.system("samtools view -bS {}/{}.sam > {}/{}.bam".format(directory, identifier, directory, identifier))
    os.system("ccs {}/{}.bam -j {} {}/{}.fastq".format(directory, identifier, cores, directory, identifier))
    return


def create_template(input_fasta, seq_path, s_id, primers, n_templates):
    primers_df = pd.read_csv(template, sep='\t', header=None)
    primer_df.columns = ["chr", "start", "end", "name_1", "score", "strand", "primer"]
    final_template = "" 

    # read template file two lines at a time (positive & negative strand info)
    for i in range(0, len(primer_df), 2):
        # get positive strand info for reference genome
        pos_strand = primer_df.iloc[i]
        # get negative strand info for reference genome
        neg_strand = primer_df.iloc[i+1]

        seq_start = pos_strand["end"]
        seq_end = neg_strand["start"]

        Fprob = pos_strand["primer"].strip()
        # print(Fprob)
        Rprob = neg_strand["primer"].strip()
        # print(Rprob)
        
        try:
            # get amplicon positions
            start, end = get_amplicon_positions(Fprob, Rprob, seq_path)

        except:
            print("Something went wrong with the amplicon positions.")
            print("Will get amplicon with MSAs instead")
            # get amplicon with MSAs
            try: 
                parent_dir = "/".join(seq_path.split("/")[: -1])
                print(input_fasta, parent_dir)
                amplicon = subprocess.run(['python -m long_amplicon_read_simulation.find_amplicons_with_mafft --input_fasta {} --directory {} --start {} --end {} --record_id {}'.format(input_fasta, parent_dir, seq_start, seq_end, s_id)], shell=True, capture_output=True, text=True)
                print("Amplicon: ", amplicon)
                amplicon = str(amplicon.stdout).strip()
                cnt = 0
                while cnt < int(n_templates/2):
                    # write hald with + strand and half with - strand
                    final_template += ">" + "+_" + s_id + ":" + str(seq_start) + "_" + str(seq_end) + ":" + str(cnt) + "\n"
                    final_template += amplicon + "\n"
                    cnt += 1
                leftover_templates = n_templates - int(n_templates/2)
                cnt2 = 0
                while cnt2 < leftover_templates:
                    final_template += ">" + "-_" + s_id + ":" + str(seq_start) + "_" + str(seq_end) + ":" + str(cnt+cnt2) + "\n"
                    final_template += str(amplicon[::-1].translate(str.maketrans("ATGC", "TACG"))) + "\n"
                    cnt2 += 1
                
                continue
            
            except Exception as e:
                print(e)
                print("Something went wrong with the MSAs.")
                print("Skipping amplicon {}, {} for sequence {}".format(seq_start, seq_end, s_id))
                continue

        # parse the fasta file
        record = SeqIO.parse(seq_path, "fasta")
        seq = next(record).seq
        
        amplicon = seq[start:end]
        cnt = 0
        
        while cnt < int(n_templates/2):
            # write hald with + strand and half with - strand
            final_template += ">" + "+_" + s_id + ":" + str(seq_start) + "_" + str(seq_end) + ":" + str(cnt) + "\n"
            final_template += str(amplicon) + "\n"
            cnt += 1
        
        leftover_templates = n_templates - int(n_templates/2)
        
        cnt2 = 0
        
        while cnt2 < leftover_templates:
            final_template += ">" + "-_" + s_id + ":" + str(seq_start) + "_" + str(seq_end) + ":" + str(cnt+cnt2) + "\n"
            final_template += str(str(amplicon)[::-1].translate(str.maketrans("ATGC", "TACG"))) + "\n"
            cnt2 += 1

        
        print("Final template: ", final_template)

    return final_template