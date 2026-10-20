"""uncertainty of classification performance and dependence-aware
significance of directional concordance.

For every held-out watermelon study (leave-one-study-out signature) and every cucumber study (CDS orthologues):
  * AUC of the signature score with an exact (or Monte-Carlo) permutation P value obtained by permuting the
    condition labels within design strata, the smallest P value attainable with that design, and a stratified
    bootstrap 95% CI where n >= 12;
  * directional concordance computed from a simple stratified mean-difference estimator, with
      - a self-contained null: the same statistic after permuting condition labels within strata (keeps the
        full gene-gene correlation structure, so it does not assume independent genes), and
      - a competitive null: concordance of 1,000 random expressed non-signature gene sets of the same size whose
        predicted direction is taken from the same training meta-analysis;
  * results with and without the low-mapping cucumber study PRJNA1150737.
Outputs: results/extended/validation_competitive.tsv
"""
import os
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

import a2_meta as MA

R = Path(os.environ.get("PROJECT_DIR", ".") + "/results")
O = R / "extended"
RNG = np.random.default_rng(2026)
MAXPERM = 2000

S = pd.read_csv(R / "samples.tsv", sep="\t", index_col=0, dtype=str).fillna("")


def norm_log(counts):
    """DESeq-style median-of-ratios size factors, then log2(normalised count + 1)."""
    c = counts.loc[(counts > 0).all(axis=1)]
    lg = np.log(c)
    sf = np.exp((lg.sub(lg.mean(axis=1), axis=0)).median(axis=0))
    return np.log2(counts / sf + 1)


def strata(ss):
    key = ss.genotype.where(ss.genotype != "", "all")
    if ss.study.iloc[0] == "WM_Leaf_M20Y34_2018":
        key = ss.genotype + "_" + ss.timepoint
    return key


def label_perms(y, key):
    """All (or MAXPERM random) label arrangements that keep the number of stressed samples per stratum."""
    groups = {k: np.where(key.values == k)[0] for k in key.unique()}
    per = []
    for k, idx in groups.items():
        n1 = int(y[idx].sum())
        per.append([(idx, set(c)) for c in itertools.combinations(idx, n1)])
    total = int(np.prod([len(p) for p in per]))
    def make(choice):
        yy = np.zeros(len(y), int)
        for (idx, sel) in choice:
            yy[list(sel)] = 1
        return yy
    if total <= MAXPERM:
        return [make(ch) for ch in itertools.product(*per)], total, True
    out = [make([p[RNG.integers(len(p))] for p in per]) for _ in range(MAXPERM)]
    return out, total, False


def lfc_simple(x, y, key):
    """Stratified difference of mean log2 normalised expression (stress - control), averaged over strata."""
    d = []
    for k in key.unique():
        m = (key == k).values
        a, b = x.loc[:, m & (y == 1)], x.loc[:, m & (y == 0)]
        if a.shape[1] and b.shape[1]:
            d.append(a.mean(axis=1) - b.mean(axis=1))
    return pd.concat(d, axis=1).mean(axis=1)


def zs(x):
    sd = x.std(axis=1).replace(0, np.nan)
    return x.sub(x.mean(axis=1), axis=0).div(sd, axis=0)


def evaluate(name, kind, x, y, key, up, dn, pool_dir):
    """x: log expression (genes x samples) of one study; up/dn: predicted genes; pool_dir: predicted sign for
    background genes (Series indexed like x)."""
    up = [g for g in up if g in x.index]
    dn = [g for g in dn if g in x.index]
    genes = up + dn
    pred = pd.Series([1] * len(up) + [-1] * len(dn), index=genes)
    z = zs(x)
    score = z.loc[up].mean() - z.loc[dn].mean()
    auc = roc_auc_score(y, score)
    perms, total, exact = label_perms(y, key)
    auc_null = np.array([roc_auc_score(yy, score) for yy in perms])
    p_auc = np.mean(auc_null >= auc - 1e-12) if exact else (np.sum(auc_null >= auc - 1e-12) + 1) / (len(perms) + 1)
    min_p = 1 / total if exact else 1 / (len(perms) + 1)
    ci = (np.nan, np.nan)
    if len(y) >= 12:
        bs = []
        for _ in range(2000):
            idx = np.concatenate([RNG.choice(np.where(key.values == k)[0], (key == k).sum(), replace=True)
                                  for k in key.unique()])
            if len(set(y[idx])) == 2:
                bs.append(roc_auc_score(y[idx], score.iloc[idx]))
        ci = tuple(np.percentile(bs, [2.5, 97.5]))
    xx = x.loc[genes]
    lf = lfc_simple(xx, y, key)
    conc = float((np.sign(lf) == pred).mean())
    conc_perm = np.array([float((np.sign(lfc_simple(xx, yy, key)) == pred).mean()) for yy in perms])
    real = [i for i, yy in enumerate(perms) if not np.array_equal(yy, y)]
    cp = conc_perm[real]
    p_conc_perm = (np.sum(cp >= conc - 1e-12) + 1) / (len(cp) + 1)
    # competitive null: random non-signature genes, predicted direction from the same training meta-analysis
    bg = [g for g in pool_dir.index if g in x.index and g not in set(genes)]
    lf_bg = lfc_simple(x.loc[bg], y, key)
    agree_bg = (np.sign(lf_bg) == np.sign(pool_dir.loc[bg])).astype(float)
    comp = np.array([agree_bg.iloc[RNG.choice(len(bg), len(genes), replace=False)].mean() for _ in range(1000)])
    binom = stats.binomtest(int(round(conc * len(genes))), len(genes), 0.5, alternative="greater").pvalue
    return dict(validation=name, kind=kind, n_samples=len(y), n_stress=int(y.sum()), n_genes=len(genes),
                AUC=auc, AUC_perm_p=p_auc, AUC_min_attainable_p=min_p, n_label_arrangements=total,
                AUC_boot_lo=ci[0], AUC_boot_hi=ci[1],
                concordance=conc, concordance_perm_mean=cp.mean(), concordance_perm_max=cp.max(),
                concordance_perm_p=p_conc_perm, binomial_p_independence=binom,
                background_concordance_mean=comp.mean(), background_concordance_q99=np.quantile(comp, 0.99),
                competitive_p=(np.sum(comp >= conc) + 1) / (len(comp) + 1))


def main():
    rows = []
    cw = pd.read_csv(R / "counts_watermelon.tsv.gz", sep="\t", index_col=0)
    for held in MA.DISCOVERY:
        train = [d for d in MA.DISCOVERY if d != held]
        m = MA.meta(train, min_k=len(train))
        up, dn = MA.signature(m, min_same=len(train))
        ss = S[S.study == held]
        x = norm_log(cw[ss.index])
        x = x[(cw[ss.index] >= 10).sum(axis=1) >= 3]
        y = (ss.condition == "stress").astype(int).values
        r = evaluate(held, "watermelon held-out (LOSO)", x, y, strata(ss).reset_index(drop=True), up, dn,
                     np.sign(m.meta_log2FC))
        rows.append(r)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    cds = pd.read_csv(R / "cds.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(R / "meta.tsv.gz", sep="\t", index_col=0)
    rbh = pd.read_csv(R.parent / "ref" / "rbh_wm_cu.tsv", sep="\t")
    o = dict(zip(rbh.watermelon, rbh.cucumber))
    cc = pd.read_csv(R / "counts_cucumber.tsv.gz", sep="\t", index_col=0)
    up = [o[g] for g in cds.index[cds.direction == "up"] if g in o]
    dn = [o[g] for g in cds.index[cds.direction == "down"] if g in o]
    bgdir = pd.Series({o[g]: np.sign(v) for g, v in meta.meta_log2FC.items() if g in o})
    for st in ["CU_Leaf_PEG_2025", "CU_Leaf_PEG48h_2024", "CU_Leaf_graft_2022"]:
        ss = S[S.study == st]
        x = norm_log(cc[ss.index])
        x = x[(cc[ss.index] >= 10).sum(axis=1) >= 2]
        y = (ss.condition == "stress").astype(int).values
        r = evaluate(st, "cucumber (orthologues)", x, y, strata(ss).reset_index(drop=True), up, dn, bgdir)
        rows.append(r)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(O / "validation_competitive.tsv", sep="\t", index=False)
    print(res.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
