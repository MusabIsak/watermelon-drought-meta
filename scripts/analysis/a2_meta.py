"""Random-effects meta-analysis across the 4 watermelon drought discovery studies and the
Core Drought Signature (CDS).

CDS rule (pre-specified): meta BH-FDR < 0.01, |meta log2FC| >= 1, same direction in >= 3 of 4 studies.
Outputs: results/meta.tsv.gz, results/cds.tsv
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

R = Path(os.environ.get("PROJECT_DIR", ".") + "/results")
DISCOVERY = ["WM_Leaf_M20Y34_2018", "WM_Leaf_CS_2020", "WM_Root_PEG_2016", "WM_Leaf_PEGtime_2024"]
FDR, LFC, MIN_SAME = 0.01, 1.0, 3


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    n = ok.sum()
    order = np.argsort(p[ok])
    q = p[ok][order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    tmp = np.empty(n)
    tmp[order] = np.minimum(q, 1)
    out[ok] = tmp
    return out


def dersimonian_laird(y, se):
    """Vectorised DL random-effects meta-analysis. y, se: genes x studies (NaN allowed)."""
    w = 1 / se ** 2
    k = np.sum(~np.isnan(y), axis=1)
    sw = np.nansum(w, axis=1)
    fe = np.nansum(w * y, axis=1) / sw
    Qs = np.nansum(w * (y - fe[:, None]) ** 2, axis=1)
    c = sw - np.nansum(w ** 2, axis=1) / sw
    tau2 = np.maximum(0, (Qs - (k - 1)) / c)
    wr = 1 / (se ** 2 + tau2[:, None])
    swr = np.nansum(wr, axis=1)
    est = np.nansum(wr * y, axis=1) / swr
    se_re = np.sqrt(1 / swr)
    z = est / se_re
    p = 2 * stats.norm.sf(np.abs(z))
    I2 = np.where(Qs > 0, np.maximum(0, (Qs - (k - 1)) / Qs), 0) * 100
    Qp = stats.chi2.sf(Qs, np.maximum(k - 1, 1))
    return est, se_re, z, p, tau2, I2, Qs, Qp, k


def meta(studies, min_k=None):
    lfc, se = {}, {}
    for s in studies:
        d = pd.read_csv(R / "de" / f"{s}.tsv", sep="\t", index_col=0)
        lfc[s] = d.log2FoldChange
        se[s] = d.lfcSE
    lfc = pd.DataFrame(lfc)
    se = pd.DataFrame(se)
    min_k = min_k or len(studies)
    keep = lfc.notna().sum(axis=1) >= min_k
    lfc, se = lfc[keep], se[keep]
    est, sre, z, p, tau2, I2, Q, Qp, k = dersimonian_laird(lfc.values, se.values)
    out = lfc.add_prefix("lfc_").copy()
    out["meta_log2FC"], out["meta_se"], out["meta_z"], out["meta_p"] = est, sre, z, p
    out["meta_fdr"] = bh(p)
    out["tau2"], out["I2"], out["Q_p"], out["k"] = tau2, I2, Qp, k
    sign = np.sign(lfc)
    out["n_up"] = (sign > 0).sum(axis=1)
    out["n_down"] = (sign < 0).sum(axis=1)
    out["n_same_dir"] = np.where(est > 0, out.n_up, out.n_down)
    return out


def signature(m, fdr=FDR, lfc=LFC, min_same=MIN_SAME):
    s = m[(m.meta_fdr < fdr) & (m.meta_log2FC.abs() >= lfc) & (m.n_same_dir >= min_same)]
    return list(s.index[s.meta_log2FC > 0]), list(s.index[s.meta_log2FC < 0])


def main():
    m = meta(DISCOVERY, min_k=3)
    m.to_csv(R / "meta.tsv.gz", sep="\t")
    up, dn = signature(m)
    cds = m.loc[up + dn].copy()
    cds["direction"] = np.where(cds.meta_log2FC > 0, "up", "down")
    cds.sort_values("meta_p").to_csv(R / "cds.tsv", sep="\t")
    print(f"genes in meta-analysis: {len(m)}; meta-FDR<0.01: {(m.meta_fdr < 0.01).sum()}")
    print(f"Core Drought Signature: {len(up)} up, {len(dn)} down (total {len(up) + len(dn)})")
    print(f"median I2 in CDS: {cds.I2.median():.1f}%")


if __name__ == "__main__":
    main()
