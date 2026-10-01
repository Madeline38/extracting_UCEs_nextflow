#!/usr/bin/env python3
import sys
from pathlib import Path
from Bio import AlignIO

fasta_in = Path(sys.argv[1])
nexus_out = fasta_in.stem + ".nexus"

AlignIO.convert(str(fasta_in), "fasta", nexus_out, "nexus", molecule_type="DNA")