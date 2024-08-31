import pandas as pd
import os
import sys
import argparse
import random
import editdistance
import csv
from utils.utils import get_genomic_regions, remove_directory_contents
from utils.file_parsing import parse_fastq_ONT

def parse_template_fasta(template_file):
    with open(template_file, "r") as f:
        lines = f.readlines()

    lines = [x.strip() for x in lines]
    l = 0
    reads = {}
    while l < len(lines):
        if lines[l].startswith(">"):
            read_name = lines[l][1:-2] # e.g. >+_EPI_ISL_2683873:54_1183:1 -> +_EPI_ISL_2683873:54_1183
            l+=1
            read = lines[l]
            l+=1
            
            l+=18 # in the template file there are 10 same reads, so skip them

            if read == "": 
                continue

            reads[read_name] = read

    return reads

def check_template_sequence(read_anchor_gisaid_id, read_negative_gisaid_id, genomic_region, directory_acnhor, directory_negative, strand_anch, strand_neg):
    # load the template sequence from which the read was simulated
    template_file_anchor = os.path.join(directory_acnhor, read_anchor_gisaid_id + ".template")
    
    template_reads = parse_template_fasta(template_file_anchor)
    
    template_file_negative = os.path.join(directory_negative, read_negative_gisaid_id + ".template")
    template_reads_negative = parse_template_fasta(template_file_negative)
    
    id_anch = strand_anch + "_" + read_anchor_gisaid_id + ":" + str(genomic_region[0]) + "_" + str(genomic_region[1])
    id_neg = strand_neg + "_" + read_negative_gisaid_id + ":" + str(genomic_region[0]) + "_" + str(genomic_region[1])
    
    if (id_neg not in template_reads_negative) or (id_anch not in template_reads):
        print(f"Either or both {id_neg} and {id_anch} are missing from templates.")
        return False

    if strand_anch == "-":
        template_reads[id_anch] = template_reads[id_anch][::-1].translate(str.maketrans("ATGC", "TACG"))
    if strand_neg == "-":
        template_reads_negative[id_neg] = template_reads_negative[id_neg][::-1].translate(str.maketrans("ATGC", "TACG"))

    ed = editdistance.eval(template_reads[id_anch],template_reads_negative[id_neg])
    if ed <= 0:
        return False

    return True

def parse_maf_file(maf_file):
    maf_dict = {}
    with open(maf_file, "r") as f:
        lines = f.readlines()

    cnt = 0
    while cnt < len(lines):
        if lines[cnt].startswith("a"):
            next_line = lines[cnt+1].strip().split(" ")
            genomic_region = next_line[1].split(":")[1]
            strand = next_line[1][0]
            next_line = lines[cnt+2].strip().split(" ")
            sim_read_id = "/".join(next_line[1].split("/")[:-1])
            if genomic_region not in maf_dict.keys():
                maf_dict[genomic_region] = set()
                maf_dict[genomic_region].add(strand + "_" + sim_read_id)
            else:
                maf_dict[genomic_region].add(strand + "_" + sim_read_id)

            cnt += 4

    return maf_dict

def parse_maf_file_ONT(maf_file):
    maf_dict = {}
    with open(maf_file, "r") as f:
        lines = f.readlines()

    cnt = 0
    while cnt < len(lines):
        if lines[cnt].startswith("a"):
            next_line = lines[cnt+1].strip().split(" ")
            genomic_region = next_line[1].split(":")[1]
            strand = next_line[1][0]
            next_line = lines[cnt+2].strip().split(" ")
            sim_read_id = next_line[1].strip()
            if genomic_region not in maf_dict.keys():
                maf_dict[genomic_region] = set()
                maf_dict[genomic_region].add(strand + "_" + sim_read_id)
            else:
                maf_dict[genomic_region].add(strand + "_" + sim_read_id)

            cnt += 4
    
    return maf_dict

def parse_fastq(fastq_file):

    reads = {}

    with open(fastq_file, "r") as f:
        lines = f.readlines()

    lines = [x.strip() for x in lines]

    l = 0
    while(l!= len(lines)):
        # @S/1/ccs -> lines[l].split(" ")[0]
        # S/1/ccs -> lines[l].split(" ")[0][1:]
        # [S, 1, css] lines[l].split(" ")[0][1:].split("/")
        # [S, 1] -> lines[l].split(" ")[0][1:].split("/")[:-1]
        # S/1 -> "/".join(lines[l].split(" ")[0][1:].split("/")[:-1])
        read_name = "/".join(lines[l].split(" ")[0][1:].split("/")[:-1])
        l+=1
        read = lines[l]
        l+=2
        # quality_score = lines[l]
        l+=1
        reads[read_name] = read
    
    return reads


def parse_fastq_ONT(fastq_file):

    reads = {}

    with open(fastq_file, "r") as f:
        lines = f.readlines()

    lines = [x.strip() for x in lines]

    l = 0
    while(l!= len(lines)):
        # @S_1
        read_name = lines[l].strip()[1:]
        l+=1
        read = lines[l]
        l+=2
        # quality_score = lines[l]
        l+=1
        reads[read_name] = read

    return reads

def sample_negative_pair(read_anchor, maf_dict_2, reads2, gr_key):
    # get negative pair of reads from the same genomic region
    # find interection between maf_dict_2[gr_key] and reads2.keys()
    # sample one of the reads from the intersection
    intersection = set([x[2:] for x in  maf_dict_2[gr_key]]).intersection(set(reads2.keys()))
    if len(intersection) == 0:
        raise RuntimeError("No reads in the intersection")

    read_negative_sid = random.sample(list(intersection),1)

    reads_negative = reads2[read_negative_sid[0]]

    # find strand information
    strand_neg = find_strand(maf_dict_2[gr_key], read_negative_sid, single=True)

    if strand_neg == "-":
        reads_negative = reads_negative[::-1].translate(str.maketrans("ATGC", "TACG"))
    
    ed_neg =  editdistance.eval(read_anchor, reads_negative)
    
    return (read_negative_sid[0], reads_negative), ed_neg, strand_neg

def check_if_in_tsv(gisaid_id_1, gisaid_id_2, read_1_sid, read_2_sid, outdir):
    # load tsv file 
    tsv_file = os.path.join(outdir, "samples.tsv")

    out_df = pd.read_csv(tsv_file, sep="\t", header=0)
    complete_read_id_1 = gisaid_id_1 + "_" + read_1_sid.replace("/", "")
    complete_read_id_2 = gisaid_id_2 + "_" + read_2_sid.replace("/", "")

    if complete_read_id_1 in out_df["read_1"].values and complete_read_id_2 in out_df["read_2"].values:
        print("positive pair is already in tsv file")
        return True
    if complete_read_id_1 in out_df["read_2"].values and complete_read_id_2 in out_df["read_1"].values:
        print("positive pair is already in tsv file")
        return True
    
    return False

def find_strand(maf_dict_for_gr, ids, single=False):   
    
    if single == True:
        
        strand = ""
        for item in maf_dict_for_gr:
            if item[2:] == ids[0]:
                strand+=item[0]
                break
        if (strand != "+" and  strand != "-"): 
            raise RuntimeError(f"Strand information not found during negative sampling. Instead found {strand}")
        
        return strand
       
    strand_pos = ""
    strand_anch = ""
    for item in maf_dict_for_gr: 
        if item[2:] == ids[0]: 
            strand_anch+=item[0]
        if item[2:] == ids[1]:
            strand_pos+=item[0]
    
    if (strand_pos != "+" and  strand_pos != "-") or (strand_anch != "+" and  strand_anch != "-"): 
        raise RuntimeError(f"strand information not found during positive sampling. Instead we found: {strand_pos} {strand_anch}")
    
    return strand_anch, strand_pos 



def sample_positive_pair(maf_dict_1, reads1, gr_key):
    # sample positive pair of reads from the same genomic region
    # find intersection between maf_dict_1[gr_key] and reads1.keys()
    # sample two reads from the intersection

    if gr_key not in maf_dict_1.keys(): 
        raise RuntimeError(f"{gr_key} not found in the .maf file. Skipping sampling from this region.")

    intersection = set([x[2:] for x in maf_dict_1[gr_key]]).intersection(set(reads1.keys()))

    if len(intersection) < 2:
        raise RuntimeError("not enough reads in the intersection")
  
    reads_positive = random.sample(list(intersection), 2)
    
    read_anchor = reads1[reads_positive[0]]
    read_anchor_sid = reads_positive[0]
    read_positive = reads1[reads_positive[1]]
    read_pos_sid = reads_positive[1]
    (strand_anch, strands_pos) = find_strand(maf_dict_1[gr_key], reads_positive)

    if strand_anch == "-":
        read_anchor = read_anchor[::-1].translate(str.maketrans("ATGC", "TACG"))
    if strands_pos == "-":
        read_positive = read_positive[::-1].translate(str.maketrans("ATGC", "TACG"))

    ed_pos = editdistance.eval(read_anchor, read_positive)

    return (read_anchor_sid, read_anchor), (read_pos_sid, read_positive), ed_pos, (strand_anch, strands_pos)


def write_to_tsv(gisaid_id_1, gisaid_id_2, read_1_sid, read_2_sid, label, genomic_region, read_length1, read_length2, ed, outdir):
    complete_read_id_1 = gisaid_id_1 + "_" + read_1_sid.replace("/", "")
    complete_read_id_2 = gisaid_id_2 + "_" + read_2_sid.replace("/", "")

    tsv_file = os.path.join(outdir, "samples.tsv")
    # load tsv file
    with open(tsv_file, 'a') as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow([complete_read_id_1, complete_read_id_2, label, genomic_region, read_length1, read_length2, ed])
        tsvfile.close()

    return 


def write_to_fasta(genome_id, read_sid, read, outdir):
    complete_read_id = genome_id + "_" + read_sid.replace("/", "")
    read_file = os.path.join(outdir, "reads", complete_read_id + ".fasta")
    with open(read_file, "w") as f:
        f.write(">" + complete_read_id + "\n")
        f.write(read + "\n")

    return

def sample(genomic_region, sample1, sample2, outdir, strategy):
    
    maf_file_1 = os.path.join(sample1[1], sample1[0] + ".maf")
    fastq_file_1 = os.path.join(sample1[1], sample1[0] + ".fastq")

    maf_file_2 = os.path.join(sample2[1], sample2[0] + ".maf")
    fastq_file_2 = os.path.join(sample2[1], sample2[0] + ".fastq")

    gr_key  = str(genomic_region[0]) + "_" + str(genomic_region[1])

    if strategy == "ONT": 
        maf_dict_1 = parse_maf_file_ONT(maf_file_1)
        maf_dict_2 = parse_maf_file_ONT(maf_file_2)
    else: 
        maf_dict_1 = parse_maf_file(maf_file_1) # genomic region -> simulated read id 
        maf_dict_2 = parse_maf_file(maf_file_2)

    if strategy == "ONT": 
        reads_1 = parse_fastq_ONT(fastq_file_1) # simulated read id -> read
        reads_2 = parse_fastq_ONT(fastq_file_2)
    else: 
        reads_1 = parse_fastq(fastq_file_1) # simulated read id -> read
        reads_2 = parse_fastq(fastq_file_2)

    # get positive pair
    (read_anchor_sid, read_anchor), (read_pos_sid, read_positive), ed_pos, (strand_anch, strands_pos) = sample_positive_pair(maf_dict_1, reads_1, gr_key)
    
    patience = 0
    # check if the positive pair is already in the tsv file
    while check_if_in_tsv(sample1[0], sample1[0], read_anchor_sid, read_pos_sid, outdir):
        # resampling positive pair
        (read_anchor_sid, read_anchor), (read_pos_sid, read_positive), ed_pos, (strand_anch, strands_pos) = sample_positive_pair(maf_dict_1, reads_1, gr_key)
        patience +=1
        if patience >= 2 and check_if_in_tsv(sample1[0], sample1[0], read_anchor_sid, read_pos_sid, outdir): 
            return 0

    # get negative pair for the anchor read
    (read_negative_sid, read_negative), ed_neg, strand_neg = sample_negative_pair(read_anchor, maf_dict_2, reads_2, gr_key)

    patience = 0

    # check if the negative pair is already in the tsv file or the negative pair is actually identical reads
    while check_if_in_tsv(sample1[0], sample2[0], read_anchor_sid, read_negative_sid, outdir) or ed_neg == 0:
        # resampling negative pair
        (read_negative_sid, read_negative), ed_neg, strand_neg = sample_negative_pair(read_anchor, maf_dict_2, reads_2, gr_key)
        patience +=1
        if patience >=2 and check_if_in_tsv(sample1[0], sample2[0], read_anchor_sid, read_negative_sid, outdir) or ed_neg == 0:
            return 0

    if not check_template_sequence(sample1[0], sample2[0], genomic_region, sample1[1], sample2[1], strand_anch, strand_neg):
        print("Either negative pair is identical or the reads do not exist. Will begin sampling a new sample")
        return 0

    if ed_neg == 0:
        print("negative pair is identical, will begin sampling a new sample")
        return 0

    if ed_neg <= ed_pos:
        write_to_tsv(sample1[0], sample1[0], read_anchor_sid, read_pos_sid, "positive", gr_key, len(read_anchor), len(read_positive), ed_pos, outdir)
        write_to_tsv(sample1[0], sample2[0], read_anchor_sid, read_negative_sid, "hard negative", gr_key, len(read_anchor), len(read_negative), ed_neg, outdir)
        write_to_fasta(sample1[0], read_anchor_sid, read_anchor, outdir)
        write_to_fasta(sample1[0], read_pos_sid, read_positive, outdir)
        write_to_fasta(sample2[0], read_negative_sid, read_negative, outdir)
        
        return  2
    
    # write to tsv file both positive and negative pair
    write_to_tsv(sample1[0], sample1[0], read_anchor_sid, read_pos_sid, "positive", gr_key, len(read_anchor), len(read_positive), ed_pos, outdir)
    write_to_tsv(sample1[0], sample2[0], read_anchor_sid, read_negative_sid, "negative", gr_key, len(read_anchor), len(read_negative), ed_neg, outdir)
    
    # write to fasta file
    write_to_fasta(sample1[0], read_anchor_sid, read_anchor, outdir)
    write_to_fasta(sample1[0], read_pos_sid, read_positive, outdir)
    write_to_fasta(sample2[0], read_negative_sid, read_negative, outdir)
    
    return 2
    

def main():
    parser = argparse.ArgumentParser(description="Split multi-fasta files into single fasta files")
    parser.add_argument('--n', dest = 'num', required=True, type=int, help="number of samples to produce")
    parser.add_argument('--data_dir', dest = 'data_dir', required=True, type=str, help="data directory to sample from")
    parser.add_argument('--out_dir', dest = 'out_dir', required=True, type=str, help="output directory to write samples to")
    parser.add_argument('--primers', dest = 'primers', required=True, type=str, help="file containing primers")
    parser.add_argument('--strategy', dest = 'strategy', default="pacbio-hifi", required=True, type=str, help="ONT or pacbio-hifi")    
    parser.add_argument('--primer_file_mode', dest = 'mode', required=True, type=int, help="")    
    args = parser.parse_args()

    data_directory = args.data_dir
    num = args.num
    out_dir = args.out_dir
    strategy = args.strategy

    remove_directory_contents(out_dir)

    # make tsv file to write the samples to
    tsv_file = os.path.join(out_dir, "samples.tsv")

    with open(tsv_file, "w") as f:
        f.write("read_1\tread_2\tlabel\tgenomic_region\tread_length_r1\tread_length_r2\tedit_distance\n")

    read_dir = os.path.join(out_dir, "reads")

    if not os.path.exists(read_dir):
        os.makedirs(read_dir)
   
    genomic_regions = get_genomic_regions(args.primers, mode=args.mode)
    fasta_files = os.listdir(data_directory)
    fasta_files = [os.path.join(data_directory,file) for file in fasta_files if (file != "msa_alignment.fasta" and file.endswith(".fasta"))]

    info = [(f.split("/")[-1][:-6], "/".join(f.split("/")[:-1])) for f in fasta_files]
    
    count = 0

    while count < num:

        genomic_region = genomic_regions[random.randint(0, len(genomic_regions) - 1)]
        
        info_sample = random.sample(info, 2)

        # sample the reads
        try:
            num_samples = sample(genomic_region, info_sample[0], info_sample[1], out_dir, strategy)
        
        except RuntimeError(): 
            print("genomic region not found in one of the genomes selected: ", genomic_region, " continuing to next iteration")
            num_samples = 0

        count += num_samples
        print("number of samples: ", count)

    print("finished sampling")
    print('Total number of samples: ', count)

if __name__ == "__main__":
    sys.exit(main())