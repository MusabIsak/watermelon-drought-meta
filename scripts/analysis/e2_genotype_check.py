"""data-driven check of the PRJNA1150737 sample labels.

The BioSample alias, treatment and cultivar fields pair the four libraries one way (DT/DS = drought-treated
WBC-23-2/DGPC-59; CT/CS = controls), whereas the experiment titles pair them another way. Treatment cannot be
checked independently of expression, but genotype can: libraries from the same cultivar share homozygous
sequence variants. Read-1 sequences from each library (first 400,000) are aligned to the cucumber nuclear genome
with minimap2 (mappy, short-read preset); at sites covered >= 8x in all four libraries, the alternative-allele
fraction is computed and libraries are compared by the distance between their allele-fraction profiles at
variable sites.
Output: results/extended/prjna1150737_genotype_check.tsv
"""
import os
import gzip
import itertools
from collections import defaultdict
from pathlib import Path

import mappy as mp
import numpy as np
import pandas as pd

LM = Path(os.environ.get("PROJECT_DIR", ".") + "/lowmap")
O = Path(os.environ.get("PROJECT_DIR", ".") + "/results/extended")
RUNS = ["SRR30555221", "SRR30555222", "SRR30555223", "SRR30555224"]
N = 400000

aligner = mp.Aligner(str(LM / "cu_all.fa"), preset="sr")
ref = {}


def refseq(c):
    if c not in ref:
        ref[c] = aligner.seq(c)
    return ref[c]


alt = {r: defaultdict(lambda: defaultdict(int)) for r in RUNS}
alns = {r: [] for r in RUNS}
for r in RUNS:
    n = 0
    for name, seq, qual in mp.fastx_read(str(LM / f"{r}_1.fq")):
        n += 1
        if n > N:
            break
        h = next(iter(aligner.map(seq, cs=True)), None)
        if h is None or h.mapq < 30 or h.ctg.startswith(("chloroplast", "mito")):
            continue
        alns[r].append((h.ctg, h.r_st, h.r_en))
        pos, cs, i = h.r_st, h.cs, 0
        while i < len(cs):
            op = cs[i]
            j = i + 1
            while j < len(cs) and cs[j] not in ":*+-~":
                j += 1
            arg = cs[i + 1:j]
            if op == ":":
                pos += int(arg)
            elif op == "*":
                alt[r][(h.ctg, pos)][arg[1].upper()] += 1
                pos += 1
            elif op == "-":
                pos += len(arg)
            i = j
    print(r, "nuclear alignments", len(alns[r]), flush=True)

cand = sorted({k for r in RUNS for k, d in alt[r].items() if sum(d.values()) >= 2})
by_chr = defaultdict(list)
for c, p in cand:
    by_chr[c].append(p)
by_chr = {c: np.array(sorted(v)) for c, v in by_chr.items()}
cov = {r: defaultdict(int) for r in RUNS}
for r in RUNS:
    for c, a, b in alns[r]:
        arr = by_chr.get(c)
        if arr is None:
            continue
        lo, hi = np.searchsorted(arr, a), np.searchsorted(arr, b)
        for p in arr[lo:hi]:
            cov[r][(c, int(p))] += 1
sites = [k for k in cand if all(cov[r].get(k, 0) >= 8 for r in RUNS)]
rows = []
for k in sites:
    tot = {r: cov[r][k] for r in RUNS}
    allele = max(set().union(*[alt[r][k].keys() for r in RUNS]), key=lambda b: sum(alt[r][k].get(b, 0) for r in RUNS))
    af = {r: alt[r][k].get(allele, 0) / tot[r] for r in RUNS}
    if max(af.values()) >= 0.8 and min(af.values()) <= 0.2:     # homozygous difference between libraries
        rows.append(dict(chrom=k[0], pos=k[1], alt=allele, **{f"af_{r}": af[r] for r in RUNS},
                         **{f"depth_{r}": tot[r] for r in RUNS}))
df = pd.DataFrame(rows)
df.to_csv(O / "prjna1150737_genotype_sites.tsv", sep="\t", index=False)
print("sites covered >=8x in all four:", len(sites), "; homozygous-differential sites:", len(df))
A = df[[f"af_{r}" for r in RUNS]].values > 0.5
out = []
for a, b in itertools.combinations(range(4), 2):
    out.append(dict(pair=f"{RUNS[a]}-{RUNS[b]}", identical_genotype_fraction=float((A[:, a] == A[:, b]).mean()),
                    n_sites=len(df)))
res = pd.DataFrame(out)
res.to_csv(O / "prjna1150737_genotype_check.tsv", sep="\t", index=False)
print(res.to_string(index=False))
