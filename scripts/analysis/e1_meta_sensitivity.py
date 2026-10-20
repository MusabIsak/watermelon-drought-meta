"""sensitivity of the Core Drought Signature (CDS) to the meta-analysis model,
the heterogeneity estimator, the selection thresholds and the inclusion of individual discovery studies.

For every variant we report the size of the signature, its overlap with the published CDS, and whether three
downstream conclusions hold: (i) paired-ABRE enrichment among induced genes, (ii) direction agreement in the
cucumber PEG study (PRJDB37115, highest-quality cucumber data), (iii) agreement with salt (leaf) but not cold.
Outputs: results/extended/meta_sensitivity.tsv, results/extended/meta_estimators_genes.tsv
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import fisher_exact

import a2_meta as MA

R = Path(os.environ.get("PROJECT_DIR", ".") + "/results")
O = R / "extended"
O.mkdir(exist_ok=True)


def load(studies, de_override=None):
    lfc, se = {}, {}
    for s in studies:
        f = (de_override or {}).get(s, R / "de" / f"{s}.tsv")
        d = pd.read_csv(f, sep="\t", index_col=0)
        lfc[s], se[s] = d.log2FoldChange, d.lfcSE
    return pd.DataFrame(lfc), pd.DataFrame(se)


def tau2_pm(y, v, iters=100):
    """Paule-Mandel tau^2 per gene (vectorised bisection on the generalised Q = k-1 equation)."""
    k = np.sum(~np.isnan(y), axis=1)
    lo = np.zeros(len(y))
    hi = np.full(len(y), 50.0)

    def q(t):
        w = 1 / (v + t[:, None])
        mu = np.nansum(w * y, axis=1) / np.nansum(w, axis=1)
        return np.nansum(w * (y - mu[:, None]) ** 2, axis=1) - (k - 1)
    pos0 = q(lo) > 0
    for _ in range(iters):
        mid = (lo + hi) / 2
        f = q(mid)
        lo = np.where(f > 0, mid, lo)
        hi = np.where(f > 0, hi, mid)
    return np.where(pos0, (lo + hi) / 2, 0.0)


def tau2_reml(y, v, iters=200):
    """REML tau^2 per gene by fixed-point iteration (Viechtbauer 2005), vectorised."""
    t = np.maximum(0, np.nanvar(y, axis=1, ddof=1) - np.nanmean(v, axis=1))
    for _ in range(iters):
        w = 1 / (v + t[:, None])
        sw = np.nansum(w, axis=1)
        mu = np.nansum(w * y, axis=1) / sw
        num = np.nansum(w ** 2 * ((y - mu[:, None]) ** 2 - v), axis=1)
        new = np.maximum(0, num / np.nansum(w ** 2, axis=1) + 1 / sw)
        if np.nanmax(np.abs(new - t)) < 1e-8:
            t = new
            break
        t = new
    return t


def combine(lfc, se, method="DL", min_k=3):
    keep = lfc.notna().sum(axis=1) >= min_k
    y, s = lfc[keep].values, se[keep].values
    v = s ** 2
    k = np.sum(~np.isnan(y), axis=1)
    if method == "FE":
        t2 = np.zeros(len(y))
    elif method in ("DL", "HKSJ"):
        t2 = MA.dersimonian_laird(y, s)[4]
    elif method == "REML":
        t2 = tau2_reml(y, v)
    elif method == "PM":
        t2 = tau2_pm(y, v)
    w = 1 / (v + t2[:, None])
    sw = np.nansum(w, axis=1)
    est = np.nansum(w * y, axis=1) / sw
    se_re = np.sqrt(1 / sw)
    if method == "HKSJ":
        # Hartung-Knapp-Sidik-Jonkman: rescaled variance, t distribution with k-1 df (conservative for small k)
        qhk = np.nansum(w * (y - est[:, None]) ** 2, axis=1) / np.maximum(k - 1, 1)
        se_hk = np.sqrt(np.maximum(qhk, 1) / sw)   # truncated at 1 (Knapp-Hartung ad hoc variance correction)
        p = 2 * stats.t.sf(np.abs(est / se_hk), np.maximum(k - 1, 1))
    else:
        p = 2 * stats.norm.sf(np.abs(est / se_re))
    out = pd.DataFrame(index=lfc.index[keep])
    out["meta_log2FC"], out["meta_p"], out["meta_fdr"], out["tau2"] = est, p, MA.bh(p), t2
    sign = np.sign(lfc[keep])
    out["n_same_dir"] = np.where(est > 0, (sign > 0).sum(axis=1), (sign < 0).sum(axis=1))
    return out


def stouffer(lfc, se, min_k=3):
    """Alternative non-model-based combination: unweighted Stouffer z of per-study Wald statistics."""
    keep = lfc.notna().sum(axis=1) >= min_k
    z = (lfc[keep] / se[keep])
    k = z.notna().sum(axis=1)
    zs = z.sum(axis=1) / np.sqrt(k)
    p = 2 * stats.norm.sf(np.abs(zs))
    out = pd.DataFrame(index=z.index)
    # effect size for the |lfc| filter: inverse-variance (fixed-effect) mean
    w = 1 / se[keep] ** 2
    out["meta_log2FC"] = (w * lfc[keep]).sum(axis=1) / w.where(lfc[keep].notna()).sum(axis=1)
    out["meta_p"], out["meta_fdr"] = p, MA.bh(p)
    sign = np.sign(lfc[keep])
    out["n_same_dir"] = np.where(out.meta_log2FC > 0, (sign > 0).sum(axis=1), (sign < 0).sum(axis=1))
    return out


def sig(m, fdr=0.01, lfc=1.0, min_same=3):
    s = m[(m.meta_fdr < fdr) & (m.meta_log2FC.abs() >= lfc) & (m.n_same_dir >= min_same)]
    return set(s.index[s.meta_log2FC > 0]), set(s.index[s.meta_log2FC < 0])


CNT = pd.read_csv(R / "function" / "promoter_motif_counts.tsv", sep="\t", index_col=0)
BG = [g for g in pd.read_csv(R / "meta.tsv.gz", sep="\t", index_col=0).index if g in CNT.index]
RBH = pd.read_csv(R.parent / "ref" / "rbh_wm_cu.tsv", sep="\t")
O2C = dict(zip(RBH.watermelon, RBH.cucumber))
CU = pd.read_csv(R / "de" / "CU_Leaf_PEG_2025.tsv", sep="\t", index_col=0).log2FoldChange
SALT = pd.read_csv(R / "de" / "WM_salt_LR_2025__leaf.tsv", sep="\t", index_col=0).log2FoldChange
COLD = pd.read_csv(R / "de" / "WM_Leaf_cold_2017.tsv", sep="\t", index_col=0).log2FoldChange


def downstream(up, dn):
    gs = [g for g in up if g in CNT.index]
    a = int((CNT.loc[gs, "ABRE"] >= 2).sum())
    c = int((CNT.loc[BG, "ABRE"] >= 2).sum()) - a
    OR, p = fisher_exact([[a, len(gs) - a], [c, len(BG) - len(gs) - c]])
    sgn = pd.Series({**{g: 1 for g in up}, **{g: -1 for g in dn}})

    def agree(series, mapper=None):
        vals = []
        for g, s in sgn.items():
            key = mapper.get(g) if mapper else g
            if key is not None and key in series.index and pd.notna(series[key]) and series[key] != 0:
                vals.append(np.sign(series[key]) == s)
        return np.mean(vals) if vals else np.nan, len(vals)
    cu, ncu = agree(CU, O2C)
    sa, _ = agree(SALT)
    co, _ = agree(COLD)
    return dict(ABRE_OR_up=OR, ABRE_p_up=p, cucumber_PEG_agreement=cu, cucumber_n=ncu,
                salt_leaf_agreement=sa, cold_agreement=co)


def main():
    lfc, se = load(MA.DISCOVERY)
    base = MA.meta(MA.DISCOVERY, min_k=3)
    up0, dn0 = sig(base)
    ref = up0 | dn0
    rows = []

    def record(name, group, up, dn):
        s = up | dn
        r = dict(variant=name, group=group, n_cds=len(s), n_up=len(up), n_down=len(dn),
                 retained_from_original=len(s & ref) / len(ref), jaccard=len(s & ref) / len(s | ref) if s else 0,
                 direction_conflicts=len((up & dn0) | (dn & up0)))
        r.update(downstream(up, dn))
        rows.append(r)
        print(f"  {name:45s} n={len(s):4d} retained={r['retained_from_original']:.2f} J={r['jaccard']:.2f}", flush=True)

    record("Original (DL, FDR<0.01, |lfc|>=1, >=3/4 same)", "original", up0, dn0)
    genes = {}
    for meth in ["FE", "DL", "REML", "PM", "HKSJ"]:
        m = combine(lfc, se, meth)
        genes[meth] = m.meta_fdr
        record(f"Estimator: {meth}", "estimator", *sig(m))
    st = stouffer(lfc, se)
    genes["Stouffer"] = st.meta_fdr
    record("Combination: Stouffer z (unweighted)", "estimator", *sig(st))
    pd.DataFrame(genes).to_csv(O / "meta_estimators_genes.tsv", sep="\t")
    for fdr in [0.05, 0.01, 0.001]:
        for l in [0.585, 1.0, 1.5]:
            for ms in [3, 4]:
                if (fdr, l, ms) == (0.01, 1.0, 3):
                    continue
                record(f"Thresholds: FDR<{fdr}, |lfc|>={l}, >={ms}/4 same", "threshold", *sig(base, fdr, l, ms))
    # study influence: drop each discovery study (3 remaining, all must agree)
    for s in MA.DISCOVERY:
        rest = [x for x in MA.DISCOVERY if x != s]
        m = MA.meta(rest, min_k=3)
        record(f"Leave out: {s}", "leave-one-study-out", *sig(m, min_same=3))
    # first-4M quantification of PRJNA454040 instead of random-8M
    alt = R / "robustness" / "de" / "WM_Leaf_M20Y34_2018__first4M.tsv"
    if alt.exists():
        l2, s2 = load(MA.DISCOVERY, {"WM_Leaf_M20Y34_2018": alt})
        m = combine(l2, s2, "DL")
        record("PRJNA454040 quantified from first 4 M reads", "subsampling", *sig(m))
    alt2 = O / "de" / "WM_Root_PEG_2016__random4M.tsv"
    if alt2.exists():
        l2, s2 = load(MA.DISCOVERY, {"WM_Root_PEG_2016": alt2})
        m = combine(l2, s2, "DL")
        record("PRJNA326331 quantified from random 4 M pairs", "subsampling", *sig(m))
    res = pd.DataFrame(rows)
    res.to_csv(O / "meta_sensitivity.tsv", sep="\t", index=False)
    print(res.drop(columns=["ABRE_p_up"]).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
