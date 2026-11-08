#!/usr/bin/env python3
# where : Process EXTRACT_FLANKED : extract_flanked.py ${clean_sam} ${sample} ${flank_arg}  
# usage : This script extract sequences from a SAM file and can extract only the flanking regions of the alignment, based on the CIGAR string.

import sys
import re

def clip_lengths(cigar):
    left_match = re.match(r'^(\d+)S', cigar)
    left = int(left_match.group(1)) if left_match else 0

    right_match = re.search(r'(\d+)S$', cigar)
    right = int(right_match.group(1)) if right_match else 0

    return left, right


def main(sam_path, sample, flank):
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue

            fields = line.rstrip("\n").split("\t")
            qname, flag, locus, pos, mapq, cigar = fields[0:6]
            seq = fields[9]

            if cigar == "*":
                continue

            if flank is not None:
                left_clip, right_clip = clip_lengths(cigar)
                
                start = max(0, left_clip - flank)
                end = min(len(seq), len(seq) - right_clip + flank)

            else:
                start = 0
                end = len(seq)

            extracted = seq[start:end]

            with open(f"{sample}-{locus}.fasta", "w") as out_f:
                out_f.write(f">{sample}\n{extracted.upper()}\n")


if __name__ == "__main__":
    sam_path = sys.argv[1]
    sample = sys.argv[2]
    flank = int(sys.argv[3]) if len(sys.argv) > 3 else None
    main(sam_path, sample, flank)