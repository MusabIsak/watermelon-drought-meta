"""uncertainty of classification performance and dependence-aware
tests of directional concordance.

For every validation study (four held-out watermelon drought studies, three cucumber drought studies):
  * AUC of the signature score with (i) an exact label-permutation P value (all relabelings that keep the
    group sizes within design strata), (ii) the 2.5-97.5% range of AUCs of 1,000 random gene sets, and
    (iii) a stratified bootstrap 95% CI where each group has >= 6 samples.
  * Directional concordance of signature genes with two dependence-aware tests:
      a. sample-label permutation (self-contained test): labels are permuted within design strata, the
         log2 fold change of every signature gene is recomputed from size-factor-normalised counts, and the
         fraction of genes agreeing with the predicted direction is compared with the observed fraction.
         Gene-gene correlation is preserved because whole samples are permuted.
      b. effective number of independent genes (M_eff, Li & Ji 2005) estimated from the eigenvalues of the
         gene-gene correlation matrix of the signature genes in the corresponding compendium; the binomial
         test is repeated with n = M_eff (a deliberately conservative bound).
Outputs: results/extended/validation_uncertainty.tsv
"""
import os
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

import a2_meta as MA
import a3_validate as V

R = Path(os.environ.get("PROJECT_DIR", ".") + "/results")
O = R / "extended"
O.mkdir(exist_ok=True)
RNG = np.random.default_rng(7)

STRATA = {"WM_Leaf_M20Y34_2018": ["genotype", "timepoint"], "CU_Leaf_PEG48h_2024": ["genotype"],
          "CU_Leaf_graft_2022": ["genotype"]}


def relabelings(meta, max_n=5000):
    """All (or up to max_n random) label vectors with the same number of stress samples within each stratum."""
    keys = STRATA.get(meta.study.iloc[0], [])
    groups = [meta.index] if not keys else [g.index for _, g in meta.groupby(keys)]
    per = []
    for g in groups:
        k = int((meta.loc[g, "condition"] == "stress").sum())
        per.append([set(c) for c in itertools.combinations(g, k)])
    total = np.prod([len(p) for p in per])
    combos = itertools.product(*per) if total <= max_n else (
        tuple(p[RNG.integers(len(p))] for p in per) for _ in range(max_n))
    for c in combos:
        s = set().union(*c)
        yield np.array([1 if x in s else 0 for x in meta.index]), total


def size_factors(c):
    lg = np.log(c.replace(0, np.nan))
    ref = lg.mean(axis=1)
    ok = ref.notna() & np.isfinite(ref)
    return np.exp((lg[ok].sub(ref[ok], axis=0)).median(axis=0))


def lfc_simple(logn, y, meta):
    keys = STRATA.get(meta.study.iloc[0], [])
    if not keys:
        return logn.loc[:, y == 1].mean(axis=1) - logn.loc[:, y == 0].mean(axis=1)
    parts = []
    for _, g in meta.groupby(keys):
        idx = [meta.index.get_loc(i) for i in g.index]
        yy = y[idx]
        if yy.sum() == 0 or yy.sum() == len(yy):
            continue
        cols = g.index
        parts.append(logn[cols[yy == 1]].mean(axis=1) - logn[cols[yy == 0]].mean(axis=1))
    return pd.concat(parts, axis=1).mean(axis=1)


def m_eff(x):
    """Li & Ji (2005) effective number of independent variables from a genes x samples matrix."""
    x = x.dropna()
    c = np.corrcoef(x.values)
    ev = np.abs(np.linalg.eigvalsh(np.nan_to_num(c)))
    return float(np.sum((ev >= 1).astype(float) + (ev - np.floor(ev))))


def analyse(study, meta, counts, z, up, dn, comp_z, expressed):
    y = (meta.condition == "stress").astype(int).values
    sc = V.score(z[meta.index], up, dn).values
    auc_obs = roc_auc_score(y, sc)
    perm_auc, perm_conc = [], []
    sf = size_factors(counts[meta.index])
    logn = np.log2(counts[meta.index] / sf + 1)
    genes = [g for g in up + dn if g in logn.index and logn.loc[g].sum() > 0]
    sign = pd.Series({g: (1 if g in up else -1) for g in genes})
    conc_obs = float((np.sign(lfc_simple(logn.loc[genes], y, meta)) == sign).mean())
    total = None
    for yy, total in relabelings(meta):
        perm_auc.append(roc_auc_score(yy, sc))
        perm_conc.append(float((np.sign(lfc_simple(logn.loc[genes], yy, meta)) == sign).mean()))
    perm_auc, perm_conc = np.array(perm_auc), np.array(perm_conc)
    null = V.null_auc(z[meta.index], pd.Series(y, index=meta.index), len(up), len(dn), expressed)
    ci = (np.nan, np.nan)
    if min(y.sum(), (1 - y).sum()) >= 6:
        boots = []
        pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
        for _ in range(2000):
            i = np.concatenate([RNG.choice(pos, len(pos)), RNG.choice(neg, len(neg))])
            boots.append(roc_auc_score(y[i], sc[i]))
        ci = tuple(np.percentile(boots, [2.5, 97.5]))
    meff = m_eff(comp_z.reindex(genes))
    k_eff = round(conc_obs * meff)
    return dict(study=study, n_control=int((y == 0).sum()), n_stress=int(y.sum()), n_relabelings=int(total),
                AUC=auc_obs, AUC_perm_p=float(np.mean(perm_auc >= auc_obs - 1e-12)),
                AUC_boot_lo=ci[0], AUC_boot_hi=ci[1],
                null_AUC_mean=float(np.nanmean(null)), null_AUC_q025=float(np.nanpercentile(null, 2.5)),
                null_AUC_q975=float(np.nanpercentile(null, 97.5)),
                null_AUC_p=float((np.sum(null >= auc_obs) + 1) / (len(null) + 1)),
                n_genes=len(genes), concordance=conc_obs,
                concordance_perm_max=float(perm_conc[np.abs(perm_conc - conc_obs) > 1e-12].max()) if (np.abs(perm_conc - conc_obs) > 1e-12).any() else np.nan,
                concordance_perm_p=float(np.mean(perm_conc >= conc_obs - 1e-12)),
                min_attainable_p=1 / total,
                M_eff=meff, binom_p_independent=stats.binomtest(round(conc_obs * len(genes)), len(genes), 0.5, "greater").pvalue,
                binom_p_Meff=stats.binomtest(int(k_eff), int(round(meff)), 0.5, "greater").pvalue)


def main():
    samples = pd.read_csv(R / "samples.tsv", sep="\t", index_col=0, dtype=str).fillna("")
    rows = []
    # watermelon held-out studies
    s = samples[samples.study.isin(MA.DISCOVERY)]
    tpm = pd.read_csv(R / "tpm_watermelon.tsv.gz", sep="\t", index_col=0)
    counts = pd.read_csv(R / "counts_watermelon.tsv.gz", sep="\t", index_col=0)
    z = V.zscore_within(tpm[s.index], s)
    expressed = z.dropna().index
    for held in MA.DISCOVERY:
        train = [d for d in MA.DISCOVERY if d != held]
        up, dn = MA.signature(MA.meta(train, min_k=len(train)), min_same=len(train))
        up = [g for g in up if g in expressed]
        dn = [g for g in dn if g in expressed]
        comp = z[s.index[s.study != held]]          # M_eff from the training compendium
        r = analyse(held, s[s.study == held], counts, z, up, dn, comp, expressed)
        rows.append(r)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    # cucumber
    cds = pd.read_csv(R / "cds.tsv", sep="\t", index_col=0)
    rbh = pd.read_csv(R.parent / "ref" / "rbh_wm_cu.tsv", sep="\t")
    o = dict(zip(rbh.watermelon, rbh.cucumber))
    sc_ = samples[samples.crop == "cucumber"]
    tpm_c = pd.read_csv(R / "tpm_cucumber.tsv.gz", sep="\t", index_col=0)
    cnt_c = pd.read_csv(R / "counts_cucumber.tsv.gz", sep="\t", index_col=0)
    zc = V.zscore_within(tpm_c[sc_.index], sc_)
    exp_c = zc.dropna().index
    up = [o[g] for g in cds.index[cds.direction == "up"] if g in o and o[g] in exp_c]
    dn = [o[g] for g in cds.index[cds.direction == "down"] if g in o and o[g] in exp_c]
    zw = V.zscore_within(tpm[s.index], s)
    inv = {v: k for k, v in o.items()}
    comp_w = zw.rename(index=o)                       # M_eff from the watermelon discovery compendium
    for st in sorted(sc_.study.unique()):
        r = analyse(st, sc_[sc_.study == st], cnt_c, zc, up, dn, comp_w, exp_c)
        rows.append(r)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    pd.DataFrame(rows).to_csv(O / "validation_uncertainty.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
