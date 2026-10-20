"""why does cucumber PRJNA1150737 map poorly to the transcriptome?
200,000 read-1 sequences per library are classified hierarchically with minimap2 (mappy, short-read preset):
rRNA (Ensembl ncRNA rRNA genes) -> chloroplast genome (DQ119058) -> cDNA (ASM407v2) -> genome (intron/intergenic)
-> unmapped. Read quality, adapter content and GC are also summarised. Output: results/extended/lowmap_classification.tsv
"""
import os
import gzip
import sys
from pathlib import Path

import mappy as mp
import numpy as np
import pandas as pd

REF = Path(os.environ.get("PROJECT_DIR", ".") + "/ref")
W = Path(os.environ.get("PROJECT_DIR", ".") + "/lowmap")
OUT = Path(os.environ.get("PROJECT_DIR", ".") + "/results/extended")


def write_rrna():
    f = W / "rrna_cp.fa"
    if f.exists():
        return f
    with open(f, "w") as o:
        keep = False
        for l in gzip.open(REF / "cu_ncrna.fa.gz", "rt"):
            if l.startswith(">"):
                keep = "transcript_biotype:rRNA" in l
            if keep:
                o.write(l)
        for l in open(REF / "cu_cp.fa"):
            o.write(">chloroplast\n" if l.startswith(">") else l)
    return f


def main(runs):
    print("building indexes", flush=True)
    idx_r = mp.Aligner(str(write_rrna()), preset="sr")
    idx_c = mp.Aligner(str(REF / "Cucumis_sativus.ASM407v2.cdna.all.fa.gz"), preset="sr")
    idx_g = mp.Aligner(str(REF / "cu_genome.fa.gz"), preset="sr")
    rows = []
    for run in runs:
        c = dict(rRNA=0, chloroplast=0, mRNA=0, genome_other=0, unmapped=0)
        q, gc, ad, ln = [], [], 0, []
        n = 0
        for name, seq, qual in mp.fastx_read(str(W / f"{run}.fq")):
            n += 1
            ln.append(len(seq))
            q.append(np.mean([ord(x) - 33 for x in qual]))
            gc.append((seq.count("G") + seq.count("C")) / max(1, len(seq)))
            ad += "AGATCGGAAGAGC" in seq
            hits = list(idx_r.map(seq))
            if hits:
                c["chloroplast" if hits[0].ctg == "chloroplast" else "rRNA"] += 1
                continue
            if any(True for _ in idx_c.map(seq)):
                c["mRNA"] += 1
                continue
            if any(True for _ in idx_g.map(seq)):
                c["genome_other"] += 1
                continue
            c["unmapped"] += 1
        r = {"run": run, "reads": n, **{k: v / n * 100 for k, v in c.items()},
             "mean_read_length": np.mean(ln), "mean_base_quality": np.mean(q), "pct_reads_Q<20": np.mean(np.array(q) < 20) * 100,
             "mean_GC": np.mean(gc) * 100, "pct_adapter": ad / n * 100}
        rows.append(r)
        print({k: (round(v, 1) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    pd.DataFrame(rows).to_csv(OUT / "lowmap_classification.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main(sys.argv[1:])
