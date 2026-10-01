#!/usr/bin/env python3
import sys
from pathlib import Path
from Bio import AlignIO

fasta_in = Path(sys.argv[1])
phylip_out = fasta_in.stem + ".phy"

AlignIO.convert(str(fasta_in), "fasta", phylip_out, "phylip-relaxed")