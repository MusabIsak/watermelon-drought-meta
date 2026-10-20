"""Validation analyses for the Core Drought Signature (CDS).

1. Leave-one-study-out (LOSO) within watermelon drought: signature rebuilt without the held-out study;
   score AUC vs 1,000 random gene sets of equal size (empirical p); ML benchmark under the same scheme.
2. Cross-species: CDS mapped to cucumber 1:1 orthologs; score AUC per cucumber study; correlation of
   meta log2FC with cucumber log2FC.
3. Specificity: CDS behaviour in salt and cold contrasts; drought-specific vs shared genes.
4. Time course: per-time-point log2FC of CDS genes (1, 6, 12, 24 h PEG).
5. Tolerance: CDS response in tolerant M20 vs sensitive Y34.
Outputs in results/validation/.
"""
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.svm import SVC
from xgboost import XGBClassifier

import a2_meta as MA

warnings.filterwarnings("ignore")
M = Path(os.environ.get("PROJECT_DIR", ".") + "")
R = M / "results"
V = R / "validation"
V.mkdir(exist_ok=True)
RNG = np.random.default_rng(2026)
N_NULL = 1000


def zscore_within(tpm, samples, genes=None):
    """log2(TPM+1), z-scored per gene within each study (no labels used)."""
    x = np.log2(tpm + 1)
    out = []
    for _, s in samples.groupby("study"):
        sub = x[s.index]
        sd = sub.std(axis=1).replace(0, np.nan)
        out.append(sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0))
    z = pd.concat(out, axis=1)[samples.index]
    return z if genes is None else z.reindex(genes)


def score(z, up, dn):
    u = z.reindex(up).mean(axis=0) if len(up) else 0
    d = z.reindex(dn).mean(axis=0) if len(dn) else 0
    return u - d


def auc(y, s):
    return roc_auc_score(y, s) if len(set(y)) == 2 else np.nan


def null_auc(z, y, n_up, n_dn, pool):
    vals = []
    pool = np.array(pool)
    for _ in range(N_NULL):
        g = RNG.choice(pool, n_up + n_dn, replace=False)
        vals.append(auc(y, score(z, list(g[:n_up]), list(g[n_up:]))))
    return np.array(vals)


def loso(samples, tpm):
    disc = MA.DISCOVERY
    s = samples[samples.study.isin(disc)]
    z = zscore_within(tpm[s.index], s)
    y_all = (s.condition == "stress").astype(int)
    expressed = z.dropna().index
    rows = []
    for held in disc:
        train = [d for d in disc if d != held]
        m = MA.meta(train, min_k=len(train))
        up, dn = MA.signature(m, min_same=len(train))       # all training studies agree
        up = [g for g in up if g in expressed]
        dn = [g for g in dn if g in expressed]
        t = s.study == held
        y = y_all[t]
        sc = score(z.loc[:, t], up, dn)
        a = auc(y, sc)
        nul = null_auc(z.loc[:, t], y, len(up), len(dn), expressed)
        # direction concordance in the held-out study (independent DESeq2 of that study only)
        de = pd.read_csv(R / "de" / f"{held}.tsv", sep="\t", index_col=0)
        sig_genes = [g for g in up + dn if g in de.index]
        exp_sign = pd.Series({g: (1 if g in up else -1) for g in sig_genes})
        obs = np.sign(de.loc[sig_genes, "log2FoldChange"])
        agree = int((obs == exp_sign).sum())
        binom_p = stats.binomtest(agree, len(sig_genes), 0.5, alternative="greater").pvalue
        sig_held = (de.loc[sig_genes, "padj"] < 0.05) & (obs == exp_sign)
        rho = stats.spearmanr(m.loc[sig_genes, "meta_log2FC"], de.loc[sig_genes, "log2FoldChange"]).correlation
        rows.append(dict(held_out=held, model="CDS score", n_up=len(up), n_down=len(dn), AUC=a,
                         null_mean=np.nanmean(nul), null_p=(np.sum(nul >= a) + 1) / (len(nul) + 1),
                         accuracy=np.mean((sc > 0).astype(int) == y), direction_agreement=agree / len(sig_genes),
                         n_genes_tested=len(sig_genes), binom_p=binom_p, frac_sig_same_dir=sig_held.mean(),
                         spearman_lfc=rho))
        # ML benchmark
        X = z.dropna().T
        Xtr, ytr, Xte = X[~t.values], y_all[~t], X[t.values]
        F, _ = f_classif(Xtr, ytr)
        feats = X.columns[np.argsort(-np.nan_to_num(F))[:500]]
        models = {
            "RandomForest": RandomForestClassifier(n_estimators=1000, random_state=1, n_jobs=-1),
            "XGBoost": XGBClassifier(n_estimators=400, max_depth=2, learning_rate=0.05, subsample=0.8,
                                     colsample_bytree=0.5, random_state=1, n_jobs=-1, eval_metric="logloss"),
            "LogReg_L1": LogisticRegression(penalty="l1", C=0.5, solver="liblinear", max_iter=5000),
            "SVM_linear": SVC(kernel="linear", C=0.1, probability=True, random_state=1),
        }
        for name, mdl in models.items():
            mdl.fit(Xtr[feats], ytr)
            p = mdl.predict_proba(Xte[feats])[:, 1]
            rows.append(dict(held_out=held, model=name, AUC=auc(y, p), accuracy=np.mean((p >= .5) == y)))
        r0 = rows[-5]
        print(f"  LOSO {held}: signature {len(up)}+{len(dn)} genes; direction agreement {r0['direction_agreement']:.1%} "
              f"(binomial p={r0['binom_p']:.1e}), significant & concordant {r0['frac_sig_same_dir']:.1%}, rho={r0['spearman_lfc']:.2f}")
        print(f"  LOSO {held}: CDS AUC {a:.3f} (null mean {np.nanmean(nul):.2f}, p={rows[-5]['null_p']:.3f}); "
              + ", ".join(f"{r['model']} {r['AUC']:.2f}" for r in rows[-4:]), flush=True)
    pd.DataFrame(rows).to_csv(V / "loso.tsv", sep="\t", index=False)


def cross_species(samples, cds, meta):
    tpm = pd.read_csv(R / "tpm_cucumber.tsv.gz", sep="\t", index_col=0)
    rbh = pd.read_csv(M / "ref" / "rbh_wm_cu.tsv", sep="\t")
    o = dict(zip(rbh.watermelon, rbh.cucumber))
    s = samples[samples.crop == "cucumber"]
    z = zscore_within(tpm[s.index], s)
    up = [o[g] for g in cds.index[cds.direction == "up"] if g in o and o[g] in z.index]
    dn = [o[g] for g in cds.index[cds.direction == "down"] if g in o and o[g] in z.index]
    rows = []
    pool = z.dropna().index
    allsc = []
    for st, ss in s.groupby("study"):
        y = (ss.condition == "stress").astype(int)
        sc = score(z[ss.index], up, dn)
        allsc.append(pd.DataFrame(dict(score=sc, condition=ss.condition, study=st)))
        a = auc(y, sc)
        nul = null_auc(z[ss.index], y, len(up), len(dn), pool)
        # log2FC concordance for orthologs
        de = pd.read_csv(R / "de" / f"{st}.tsv", sep="\t", index_col=0)
        pairs = [(g, o[g]) for g in meta.index if g in o and o[g] in de.index]
        mm = pd.DataFrame([(meta.loc[g, "meta_log2FC"], de.loc[c, "log2FoldChange"]) for g, c in pairs],
                          columns=["wm", "cu"]).dropna()
        rho_all = stats.spearmanr(mm.wm, mm.cu).correlation
        cpairs = [(cds.loc[g, "meta_log2FC"], de.loc[o[g], "log2FoldChange"]) for g in cds.index
                  if g in o and o[g] in de.index]
        cc = pd.DataFrame(cpairs, columns=["wm", "cu"]).dropna()
        agree = np.mean(np.sign(cc.wm) == np.sign(cc.cu))
        rows.append(dict(study=st, n=len(ss), n_up=len(up), n_down=len(dn), AUC=a, null_mean=np.nanmean(nul),
                         null_p=(np.sum(nul >= a) + 1) / (len(nul) + 1), rho_all_orthologs=rho_all,
                         cds_direction_agreement=agree, cds_orthologs=len(cc),
                         rho_cds=stats.spearmanr(cc.wm, cc.cu).correlation))
        print(f"  cucumber {st}: AUC {a:.2f} (null {np.nanmean(nul):.2f}, p={rows[-1]['null_p']:.3f}); "
              f"direction agreement {agree:.2f} over {len(cc)} CDS orthologs", flush=True)
    pd.DataFrame(rows).to_csv(V / "cross_species.tsv", sep="\t", index=False)
    pd.concat(allsc).to_csv(V / "cross_species_scores.tsv", sep="\t")


def specificity(cds):
    con = pd.read_csv(R / "contrasts.tsv", sep="\t")
    spec = con[con.role == "specificity"].contrast.tolist()
    tab = cds[["meta_log2FC", "direction"]].copy()
    for c in spec:
        d = pd.read_csv(R / "de" / f"{c}.tsv", sep="\t", index_col=0)
        tab[f"lfc_{c}"] = d.log2FoldChange.reindex(tab.index)
        tab[f"sig_{c}"] = ((d.padj < 0.05) & (d.log2FoldChange.abs() >= 1)).reindex(tab.index).fillna(False)
    same = []
    for g, r in tab.iterrows():
        n = 0
        for c in spec:
            if r[f"sig_{c}"] and np.sign(r[f"lfc_{c}"]) == np.sign(r.meta_log2FC):
                n += 1
        same.append(n)
    tab["n_other_stress_same_dir"] = same
    tab["class"] = np.where(tab.n_other_stress_same_dir == 0, "drought-specific", "shared abiotic")
    tab.to_csv(V / "specificity.tsv", sep="\t")
    print("  specificity:", tab["class"].value_counts().to_dict(), flush=True)
    # how much of the drought signature does each non-drought stress reproduce?
    rows = []
    for c in spec:
        lfc, sig = tab[f"lfc_{c}"], tab[f"sig_{c}"]
        ok = lfc.notna()
        sg = np.sign(tab.meta_log2FC[ok])
        rows.append(dict(contrast=c, cds_tested=int(ok.sum()),
                         direction_agreement=float((np.sign(lfc[ok]) == sg).mean()),
                         frac_sig_same_dir=float((sig[ok] & (np.sign(lfc[ok]) == sg)).mean()),
                         frac_sig_opposite=float((sig[ok] & (np.sign(lfc[ok]) != sg)).mean()),
                         spearman=stats.spearmanr(tab.meta_log2FC[ok], lfc[ok]).correlation))
    # the same statistics for the cucumber drought and held-out watermelon drought studies, for comparison
    pd.DataFrame(rows).to_csv(V / "specificity_by_contrast.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False), flush=True)


def time_course(cds):
    t = cds[["meta_log2FC", "direction"]].copy()
    for tp in ["1h", "6h", "12h", "24h"]:
        d = pd.read_csv(R / "de" / f"WM_Leaf_PEGtime_2024__{tp}.tsv", sep="\t", index_col=0)
        t[tp] = d.log2FoldChange.reindex(t.index)
        t[f"sig_{tp}"] = ((d.padj < 0.05) & (d.log2FoldChange.abs() >= 1)).reindex(t.index).fillna(False)

    def first(r):
        for tp in ["1h", "6h", "12h", "24h"]:
            if r[f"sig_{tp}"] and np.sign(r[tp]) == np.sign(r.meta_log2FC):
                return tp
        return "not significant"
    t["first_response"] = t.apply(first, axis=1)
    t.to_csv(V / "time_course.tsv", sep="\t")
    print("  first response:", t.first_response.value_counts().to_dict(), flush=True)


def tolerance(cds):
    out = cds[["meta_log2FC", "direction"]].copy()
    for g in ["M20", "Y34"]:
        d = pd.read_csv(R / "de" / f"WM_Leaf_M20Y34_2018__{g}.tsv", sep="\t", index_col=0)
        out[f"lfc_{g}"] = d.log2FoldChange.reindex(out.index)
    sgn = np.sign(out.meta_log2FC)
    out["stronger_in_tolerant"] = (out.lfc_M20 * sgn) > (out.lfc_Y34 * sgn)
    ok = out.lfc_M20.notna() & out.lfc_Y34.notna()
    w = stats.wilcoxon((out.lfc_M20 * sgn)[ok], (out.lfc_Y34 * sgn)[ok])
    out.to_csv(V / "tolerance.tsv", sep="\t")
    print(f"  tolerance: response stronger in tolerant M20 for {out.stronger_in_tolerant.mean():.0%} of CDS genes "
          f"(Wilcoxon p={w.pvalue:.2e})", flush=True)


def main(parts=("loso", "cross", "spec", "time", "tol")):
    samples = pd.read_csv(R / "samples.tsv", sep="\t", index_col=0, dtype=str).fillna("")
    cds = pd.read_csv(R / "cds.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(R / "meta.tsv.gz", sep="\t", index_col=0)
    if "loso" in parts:
        loso(samples, pd.read_csv(R / "tpm_watermelon.tsv.gz", sep="\t", index_col=0))
    if "cross" in parts:
        cross_species(samples, cds, meta)
    if "spec" in parts:
        specificity(cds)
    if "time" in parts:
        time_course(cds)
    if "tol" in parts:
        tolerance(cds)


if __name__ == "__main__":
    import sys
    main(tuple(sys.argv[1:]) or ("loso", "cross", "spec", "time", "tol"))
