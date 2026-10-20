"""Supplementary analysis S1: effect of using the first 4 M reads instead of a random 8 M subsample.

All 16 libraries of WM_Leaf_M20Y34_2018 (full files available) were re-quantified from the first 4 M reads
(quant_s1/) and the study-level DE was repeated with the identical design. Outputs: results/robustness/.
"""
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import a1_de

M = Path(os.environ.get("PROJECT_DIR", ".") + "")
R = M / "results"
O = R / "robustness"
STUDY = "WM_Leaf_M20Y34_2018"

design = pd.read_csv(M / "census" / "design.tsv", sep="\t", dtype=str).fillna("")
m = design[design.study == STUDY].set_index("run")
counts, tpm = {}, {}
for run in m.index:
    q = pd.read_csv(M / "quant_s1" / run / "quant.sf", sep="\t", index_col=0)
    g = q.index.map(lambda x: re.sub(r"\.\d+$", "", x))
    counts[run] = q.NumReads.groupby(g).sum()
    tpm[run] = q.TPM.groupby(g).sum()
counts = pd.DataFrame(counts).fillna(0).round().astype(int)
tpm = pd.DataFrame(tpm).fillna(0)

# 1) per-library agreement of expression estimates
tpm0 = pd.read_csv(R / "tpm_watermelon.tsv.gz", sep="\t", index_col=0)[m.index]
rows = []
for run in m.index:
    a, b = np.log2(tpm0[run] + 1), np.log2(tpm[run].reindex(tpm0.index).fillna(0) + 1)
    ok = (tpm0[run] >= 1) | (tpm[run].reindex(tpm0.index).fillna(0) >= 1)
    rows.append(dict(run=run, pearson_log2tpm=stats.pearsonr(a[ok], b[ok])[0],
                     spearman_tpm=stats.spearmanr(a[ok], b[ok])[0], genes=int(ok.sum())))
lib = pd.DataFrame(rows)
lib.to_csv(O / "S1_library_correlation.tsv", sep="\t", index=False)

# 2) differential expression with the identical design
a1_de.R = O  # write DE table into results/robustness/de/
(O / "de").mkdir(exist_ok=True)
info = a1_de.run_deseq(counts, m, a1_de.covariates(m), f"{STUDY}__first4M", dict(study=STUDY))
d0 = pd.read_csv(R / "de" / f"{STUDY}.tsv", sep="\t", index_col=0)
d1 = pd.read_csv(O / "de" / f"{STUDY}__first4M.tsv", sep="\t", index_col=0)
common = d0.index.intersection(d1.index)
x, y = d0.loc[common, "log2FoldChange"], d1.loc[common, "log2FoldChange"]
sig0 = set(d0.index[(d0.padj < 0.05) & (d0.log2FoldChange.abs() >= 1)])
sig1 = set(d1.index[(d1.padj < 0.05) & (d1.log2FoldChange.abs() >= 1)])
cds = pd.read_csv(R / "cds.tsv", sep="\t", index_col=0)
cc = cds.index.intersection(d1.index)
sgn = np.sign(cds.loc[cc, "meta_log2FC"])
summ = dict(
    genes_compared=len(common),
    lfc_pearson=stats.pearsonr(x, y)[0], lfc_spearman=stats.spearmanr(x, y)[0],
    deg_random8M=len(sig0), deg_first4M=len(sig1), deg_overlap=len(sig0 & sig1),
    deg_jaccard=len(sig0 & sig1) / len(sig0 | sig1),
    deg_random8M_recovered=len(sig0 & sig1) / len(sig0),
    cds_tested=len(cc),
    cds_same_direction_first4M=float((np.sign(d1.loc[cc, "log2FoldChange"]) == sgn).mean()),
    cds_same_direction_random8M=float((np.sign(d0.loc[cc, "log2FoldChange"]) == sgn).mean()),
    cds_lfc_spearman=stats.spearmanr(d0.loc[cc, "log2FoldChange"], d1.loc[cc, "log2FoldChange"])[0],
    median_library_pearson=lib.pearson_log2tpm.median(), min_library_pearson=lib.pearson_log2tpm.min())
pd.Series(summ).to_csv(O / "S1_summary.tsv", sep="\t", header=False)
print(pd.Series(summ).round(4).to_string())
