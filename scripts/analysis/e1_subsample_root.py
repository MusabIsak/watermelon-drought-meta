"""PRJNA326331 (paired-end, 6 libraries) quantified from
random 4 M read pairs of the complete, MD5-verified files (quant_s2/) versus the first 4 M pairs used in the main
analysis. Study-level DE is repeated with the identical design; the CDS is rebuilt with the alternative estimates.
Outputs: results/extended/de/WM_Root_PEG_2016__random4M.tsv, results/extended/S1b_root_summary.tsv,
results/extended/S1b_root_libraries.tsv
"""
import os
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import a1_de
import a2_meta as MA

M = Path(os.environ.get("PROJECT_DIR", ".") + "")
R = M / "results"
O = R / "extended"
STUDY = "WM_Root_PEG_2016"

design = pd.read_csv(M / "census" / "design.tsv", sep="\t", dtype=str).fillna("")
m = design[design.study == STUDY].set_index("run")
counts, tpm, lib = {}, {}, []
tpm0 = pd.read_csv(R / "tpm_watermelon.tsv.gz", sep="\t", index_col=0)
for run in m.index:
    q = pd.read_csv(M / "quant_s2" / run / "quant.sf", sep="\t", index_col=0)
    g = q.index.map(lambda x: re.sub(r"\.\d+$", "", x))
    counts[run] = q.NumReads.groupby(g).sum()
    tpm[run] = q.TPM.groupby(g).sum()
    i1 = json.loads((M / "quant_s2" / run / "aux_info" / "meta_info.json").read_text())
    i0 = json.loads((M / "quant" / run / "aux_info" / "meta_info.json").read_text())
    a = np.log2(tpm0[run] + 1)
    b = np.log2(tpm[run].reindex(tpm0.index).fillna(0) + 1)
    ok = (tpm0[run] >= 1) | (tpm[run].reindex(tpm0.index).fillna(0) >= 1)
    lib.append(dict(run=run, mapping_first4M=i0["percent_mapped"], mapping_random4M=i1["percent_mapped"],
                    pearson_log2tpm=stats.pearsonr(a[ok], b[ok])[0]))
counts = pd.DataFrame(counts).fillna(0).round().astype(int)
L = pd.DataFrame(lib)
L.to_csv(O / "S1b_root_libraries.tsv", sep="\t", index=False)
print(L.round(3).to_string(index=False))

a1_de.R = O
(O / "de").mkdir(exist_ok=True)
a1_de.run_deseq(counts, m, a1_de.covariates(m), f"{STUDY}__random4M", dict(study=STUDY))
d0 = pd.read_csv(R / "de" / f"{STUDY}.tsv", sep="\t", index_col=0)
d1 = pd.read_csv(O / "de" / f"{STUDY}__random4M.tsv", sep="\t", index_col=0)
common = d0.index.intersection(d1.index)
x, y = d0.loc[common, "log2FoldChange"], d1.loc[common, "log2FoldChange"]
ok = x.notna() & y.notna()
sig0 = set(d0.index[(d0.padj < 0.05) & (d0.log2FoldChange.abs() >= 1)])
sig1 = set(d1.index[(d1.padj < 0.05) & (d1.log2FoldChange.abs() >= 1)])
cds = pd.read_csv(R / "cds.tsv", sep="\t", index_col=0)
cc = cds.index.intersection(d1.index).intersection(d0.index)
sgn = np.sign(cds.loc[cc, "meta_log2FC"])
# CDS rebuilt with the random-4M root estimates
lfc, se = {}, {}
for s in MA.DISCOVERY:
    d = d1 if s == STUDY else pd.read_csv(R / "de" / f"{s}.tsv", sep="\t", index_col=0)
    lfc[s], se[s] = d.log2FoldChange, d.lfcSE
lfc, se = pd.DataFrame(lfc), pd.DataFrame(se)
keep = lfc.notna().sum(axis=1) >= 3
est, sre, z, p, *_ = MA.dersimonian_laird(lfc[keep].values, se[keep].values)
mm = pd.DataFrame(index=lfc.index[keep])
mm["meta_log2FC"], mm["meta_fdr"] = est, MA.bh(p)
sign = np.sign(lfc[keep])
mm["n_same_dir"] = np.where(est > 0, (sign > 0).sum(axis=1), (sign < 0).sum(axis=1))
sel = mm[(mm.meta_fdr < 0.01) & (mm.meta_log2FC.abs() >= 1) & (mm.n_same_dir >= 3)].index
summ = dict(genes_compared=int(ok.sum()), genes_tested_first4M=len(d0), genes_tested_random4M=len(d1),
            lfc_pearson=stats.pearsonr(x[ok], y[ok])[0], lfc_spearman=stats.spearmanr(x[ok], y[ok])[0],
            deg_first4M=len(sig0), deg_random4M=len(sig1), deg_overlap=len(sig0 & sig1),
            deg_first4M_recovered=len(sig0 & sig1) / max(1, len(sig0)),
            cds_tested=len(cc),
            cds_same_direction_first_vs_random=float((np.sign(d0.loc[cc, "log2FoldChange"]) == np.sign(d1.loc[cc, "log2FoldChange"])).mean()),
            cds_signature_direction_first4M=float((np.sign(d0.loc[cc, "log2FoldChange"]) == sgn).mean()),
            cds_signature_direction_random4M=float((np.sign(d1.loc[cc, "log2FoldChange"]) == sgn).mean()),
            cds_lfc_spearman=stats.spearmanr(d0.loc[cc, "log2FoldChange"], d1.loc[cc, "log2FoldChange"]).correlation,
            rebuilt_cds_size=len(sel), rebuilt_cds_retained=len(set(sel) & set(cds.index)) / len(cds),
            rebuilt_cds_jaccard=len(set(sel) & set(cds.index)) / len(set(sel) | set(cds.index)),
            median_mapping_first4M=L.mapping_first4M.median(), median_mapping_random4M=L.mapping_random4M.median())
pd.Series(summ).to_csv(O / "S1b_root_summary.tsv", sep="\t", header=False)
print(pd.Series(summ).round(4).to_string())
