#!/usr/bin/env python3
import sys
import csv
from collections import defaultdict
 
 
def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].strip()
                if name.startswith("_R_"):  # préfixe ajouté par mafft --adjustdirection
                    name = name[3:]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            yield name, "".join(seq)
 
 
def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)
 
 
def main(counts_tsv, membership_tsv, out_md, fastas):
    # Totaux sur les locus
    loci = list(csv.DictReader(open(counts_tsv), delimiter="\t"))
    z_total = len(loci)
    x_kept = sum(r["status"] == "kept" for r in loci)
    y_removed = z_total - x_kept
 
    # Par sample : UCEs conservés / retirés
    kept, removed = defaultdict(int), defaultdict(int)
    for r in csv.DictReader(open(membership_tsv), delimiter="\t"):
        (kept if r["status"] == "kept" else removed)[r["sample"]] += 1
    samples = sorted(set(kept) | set(removed))
 
    # %GC par sample (sur toutes les bases de ses UCEs finaux, sans gaps ni N)
    gc, acgt = defaultdict(int), defaultdict(int)
    for fa in fastas:
        for name, seq in read_fasta(fa):
            s = seq.upper()
            gc[name] += s.count("G") + s.count("C")
            acgt[name] += sum(s.count(b) for b in "ACGT")
 
    rows1, rows2 = [], []
    for s in samples:
        total = kept[s] + removed[s]
        pct = 100.0 * removed[s] / total if total else 0.0
        rows1.append([s, kept[s], removed[s], f"{pct:.1f} %"])
        gc_pct = f"{100.0 * gc[s] / acgt[s]:.1f}" if acgt[s] else "NA"
        rows2.append([s, kept[s], gc_pct])
 
    md = [
        "Ce rapport présente les résultats obtenus à la fin du pipeline d'extraction d'UCEs.\n",
        f"On a en tout conservé **{x_kept}** UCEs et on en a retiré **{y_removed}** "
        f"lors du filtrage sur les **{z_total}** initiaux.\n",
        "## Résumé des UCEs par sample\n",
        md_table(["Sample", "Nombre d'UCEs final", "Nombre d'UCEs retirés", "% retiré"], rows1),
        "\n## Stats des UCEs finaux\n",
        md_table(["Sample", "Nombre d'UCEs final", "%GC"], rows2),
        "",
    ]
    with open(out_md, "w") as f:
        f.write("\n".join(md))
 
 
if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])