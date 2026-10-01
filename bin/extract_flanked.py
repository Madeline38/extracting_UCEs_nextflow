#!/usr/bin/env python3
import sys
import re

def clip_lengths(cigar):
    """Retourne (clip_gauche, clip_droit) : nombre de bases non alignées
    au début et à la fin du contig, d'après le CIGAR."""
    left_match = re.match(r'^(\d+)S', cigar)
    left = int(left_match.group(1)) if left_match else 0

    right_match = re.search(r'(\d+)S$', cigar)
    right = int(right_match.group(1)) if right_match else 0

    return left, right


def main(sam_path, sample, flank=160, one_file_per_locus=True):
    out = None
    if not one_file_per_locus:
        out = open(f"{sample}_extracted.fasta", "w")

    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue

            fields = line.rstrip("\n").split("\t")
            qname, flag, locus, pos, mapq, cigar = fields[0:6]
            seq = fields[9]

            if cigar == "*":
                continue

            left_clip, right_clip = clip_lengths(cigar)

            start = max(0, left_clip - flank)
            end = min(len(seq), len(seq) - right_clip + flank)
            extracted = seq[start:end]

            if one_file_per_locus:
                with open(f"{sample}-{locus}.fasta", "w") as out_f:
                    out_f.write(f">{sample}\n{extracted.upper()}\n")
            else:
                out.write(f">{sample}__{locus}__{qname}\n{extracted}\n")

    if out:
        out.close()


if __name__ == "__main__":
    sam_path = sys.argv[1]
    sample = sys.argv[2]
    flank = int(sys.argv[3]) if len(sys.argv) > 3 else 160
    main(sam_path, sample, flank=flank)