#!/usr/bin/env python3

import argparse
import os
import pysam
from utils.utils import get_genomic_regions

def find_closest_region_by_start(pos, regions):
    starts = [r[0] for r in regions]  # extract starts only
    distances = [abs(pos - s) for s in starts]
    return distances.index(min(distances))

def assign_reads(regions, sam_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    out_sams = {}

    infile = pysam.AlignmentFile(sam_file, "r")
    for read in infile:
        if read.is_unmapped:
            continue
        pos = read.reference_start
        region_idx = find_closest_region_by_start(pos, regions)
        region = regions[region_idx]
        out_path = os.path.join(output_dir, f"region_{region[0]}_{region[1]}.sam")
        if region_idx not in out_sams:
            out_sams[region_idx] = pysam.AlignmentFile(out_path, "w", header=infile.header)
        out_sams[region_idx].write(read)

    infile.close()
    for out in out_sams.values():
        out.close()

    print(f"✅ Done. Reads split into {len(out_sams)} SAM files in '{output_dir}'.")

def main():
    parser = argparse.ArgumentParser(description="Split SAM reads by closest region start position.")
    parser.add_argument("--sam", required=True, help="Input SAM file")
    parser.add_argument("--primers", required=True, help="Pickled list of (start, end) tuples")
    parser.add_argument("--output", default="split_sams", help="Output directory")
    args = parser.parse_args()
    regions = get_genomic_regions(primers_file = args.primers, mode= 2)
    assign_reads(regions, args.sam, args.output)

if __name__ == "__main__":
    main()
