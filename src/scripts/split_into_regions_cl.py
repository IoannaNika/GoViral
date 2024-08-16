
import argparse
import os
import sys
import pandas as pd

def map_to_correct_region(start, genomic_regions):
    """
    Function to map a read to the correct genomic region
    
    Args:
    start: int, start position of the read
    genomic_regions: list of tuples, genomic regions

    Returns:
    start: int, start position of the read
    end: int, end position of the read
    """

    distances = {}
    
    for i in range(len(genomic_regions)):
        # calculate distance from start points
        dist_start = abs(genomic_regions[i][0] - start)
        distances[i] = dist_start

    # find the key with minimum distance
    min_key = min(distances, key=distances.get)

    return genomic_regions[min_key][0], genomic_regions[min_key][1]


def main(): 
    parser = argparse.ArgumentParser(description='Split fastq file into genomic regions')
    parser.add_argument('--bed', type=str, help='Path to bed file from whole genome')
    parser.add_argument('--out_path', type=str, help='Path to output fastq file')
    parser.add_argument('--input_path', type=str, help='Path to input fastq file')
    parser.add_argument('--ref_seq', type=str, help='Path to reference sequence')
    args = parser.parse_args()

    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                    (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                    (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                        (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                    (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                        (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                        (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]

    ref_seq = args.ref_seq
    bed = args.bed
    out_path = args.out_path
    input_path= args.input_path

    # map reads to genomic regions
    # read bed file with pandas
    df = pd.read_csv(bed + ".bed", sep="\t", header=None)

    for start,end in genomic_regions:
        if os.path.exists("{}_{}.fastq".format(out_path, str(start)+"_"+str(end))):
            os.system("rm {}_{}.fastq".format(out_path, str(start)+"_"+str(end)))
            os.system("touch {}_{}.fastq".format(out_path, str(start)+"_"+str(end)))

    for index, row in df.iterrows():
        start = row[1]
        end = row[2]
        identifier = row[3]

        start, end = map_to_correct_region(start, genomic_regions)

        # find the idifier in the fastq file
        os.system("grep -A 1 {} {}.fastq >> {}_{}.fastq".format(identifier, input_path, out_path, str(start)+"_"+str(end)))

    for start,end in genomic_regions:
        # create sam/bam files for each region
        os.system("minimap2 -ax map-pb {} {}_{}.fastq > {}_{}.sam".format(ref_seq, out_path, str(start)+"_"+str(end), out_path, str(start)+"_"+str(end)))
        os.system("samtools view -bS {}_{}.sam > {}_{}.bam".format(out_path, str(start)+"_"+str(end), out_path, str(start)+"_"+str(end)))
        os.system("samtools sort {}_{}.bam -o {}_{}.sorted.bam".format(out_path, str(start)+"_"+str(end), out_path, str(start)+"_"+str(end)))
        os.system("samtools index {}_{}.sorted.bam".format(out_path, str(start)+"_"+str(end)))
        
if __name__ == "__main__":   
    main()