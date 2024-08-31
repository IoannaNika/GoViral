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