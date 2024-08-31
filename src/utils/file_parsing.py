from typing import Dict

def parse_fastq_ONT(fastq_file: str) -> Dict[str, str]:
    """
    Function to parse a fastq file and return a dictionary of reads with read name as key and read sequence as value

    Args:
        fastq_file: str, path to the fastq file

    Returns:
        dict, a dictionary with read name as key and read sequence as value
    """

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