"""robustness of the drought-induced and drought-repressed co-expression
modules and of their hub genes.

1. Per-study behaviour: the M2/M3 eigengene is recomputed within each discovery study separately and tested for
   association with drought (Welch t test) - a module that only reflects one study's design would fail here.
2. Per-study co-expression preservation: mean absolute pairwise correlation of the module genes within each
   study versus 1,000 random gene sets of the same size (z score and empirical P).
3. Leave-one-study-out and bootstrap reconstruction: the full network pipeline (same parameters, power re-chosen)
   is rebuilt (a) without each study and (b) on 20 bootstrap resamples of libraries within study x condition.
   For each rebuild we record the best-matching module of M3/M2 (Jaccard), and whether the top-25 hubs remain in
   the top 10% of module membership (kME) of the matched module.
Outputs: results/extended/network_per_study.tsv, network_rebuilds.tsv, network_hub_stability.tsv
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import average, fcluster
from scipy.spatial.distance import squareform

import a2_meta as MA
from a5_network import eigengene, scale_free_fit

R = Path(os.environ.get("PROJECT_DIR", ".") + "/results")
O = R / "extended"
RNG = np.random.default_rng(11)


def zexpr(tpm, s):
    x = np.log2(tpm[s.index] + 1)
    z = []
    for _, ss in s.groupby("study"):
        sub = x[ss.index]
        sd = sub.std(axis=1).replace(0, np.nan)
        z.append(sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0))
    return pd.concat(z, axis=1)


def build(z, genes):
    """Same pipeline as a5_network on a fixed gene set; returns {module: gene list} and kME table."""
    X = z.reindex(genes).dropna()
    g = X.index
    X = X.values
    S = (1 + np.corrcoef(X)) / 2
    fits = []
    for b in [8, 10, 12, 14, 16, 18, 20]:
        A = S ** b
        np.fill_diagonal(A, 0)
        fits.append((b, scale_free_fit(A.sum(axis=1))[0]))
    ok = [b for b, r in fits if r >= 0.8]
    beta = ok[0] if ok else max(fits, key=lambda t: t[1])[0]
    A = S ** beta
    np.fill_diagonal(A, 0)
    k = A.sum(axis=1)
    TOM = (A @ A + A) / (np.minimum.outer(k, k) + 1 - A)
    np.fill_diagonal(TOM, 1)
    D = 1 - TOM
    Z = average(squareform(D, checks=False))
    from dynamicTreeCut import cutreeHybrid
    lab = np.asarray(cutreeHybrid(Z, squareform(D, checks=False), minClusterSize=30, deepSplit=2, pamStage=False,
                                  verbose=0)["labels"])
    mods = {l: np.where(lab == l)[0] for l in set(lab) if l != 0 and (lab == l).sum() >= 30}
    MEs = {l: eigengene(X[i]) for l, i in mods.items()}
    changed = True
    while changed:
        changed = False
        ls = list(MEs)
        for a in range(len(ls)):
            for b in range(a + 1, len(ls)):
                if np.corrcoef(MEs[ls[a]], MEs[ls[b]])[0, 1] > 0.75:
                    mods[ls[a]] = np.concatenate([mods[ls[a]], mods.pop(ls[b])])
                    MEs[ls[a]] = eigengene(X[mods[ls[a]]])
                    MEs.pop(ls[b])
                    changed = True
                    break
            if changed:
                break
    out, kme = {}, {}
    for l, idx in mods.items():
        out[l] = list(g[idx])
        kme[l] = pd.Series([np.corrcoef(X[j], MEs[l])[0, 1] for j in idx], index=g[idx])
    return out, kme, beta


def main():
    samples = pd.read_csv(R / "samples.tsv", sep="\t", index_col=0, dtype=str).fillna("")
    s = samples[samples.study.isin(MA.DISCOVERY)]
    tpm = pd.read_csv(R / "tpm_watermelon.tsv.gz", sep="\t", index_col=0)
    expressed = (tpm[s.index] >= 1).mean(axis=1) >= 0.5
    z = zexpr(tpm[expressed], s).dropna()
    gm = pd.read_csv(R / "network" / "gene_modules.tsv", sep="\t")
    xl = np.log2(tpm.loc[z.index, s.index] + 1)
    xc = pd.concat([xl[ss.index].sub(xl[ss.index].mean(axis=1), axis=0) for _, ss in s.groupby("study")], axis=1)
    genes = list(xc.var(axis=1).sort_values(ascending=False).index[:5000])   # the 5,000 network genes (as a5)
    mt = pd.read_csv(R / "network" / "modules.tsv", sep="\t")
    IND, REP = mt.loc[mt.r_drought.idxmax(), "module"], mt.loc[mt.r_drought.idxmin(), "module"]
    mods = {m: list(gm.gene[gm.module == m]) for m in [IND, REP]}
    hubs = pd.read_csv(R / "network" / "hub_genes.tsv", sep="\t")
    top = {m: list(hubs[hubs.module == m].sort_values("kME", ascending=False).gene[:25]) for m in mods}
    # 1-2 per-study behaviour and preservation
    rows = []
    for st in MA.DISCOVERY:
        ss = s[s.study == st]
        x = np.log2(tpm.loc[genes, ss.index] + 1)
        x = x[x.std(axis=1) > 0]
        y = (ss.condition == "stress").values
        for m, gl in mods.items():
            gl = [g for g in gl if g in x.index]
            xm = x.loc[gl].sub(x.loc[gl].mean(axis=1), axis=0).div(x.loc[gl].std(axis=1), axis=0).values
            me = eigengene(xm)
            t, p = stats.ttest_ind(me[y], me[~y], equal_var=False)
            c = np.corrcoef(xm)
            dens = np.nanmean(np.abs(c[np.triu_indices(len(gl), 1)]))
            nul = []
            pool = np.array(x.index)
            for _ in range(1000):
                rg = RNG.choice(pool, len(gl), replace=False)
                xr = x.loc[rg].values
                cr = np.corrcoef(xr)
                nul.append(np.nanmean(np.abs(cr[np.triu_indices(len(rg), 1)])))
            nul = np.array(nul)
            sign = 1 if m == IND else -1
            rows.append(dict(study=st, module=m, n_genes=len(gl), eigengene_shift=float(me[y].mean() - me[~y].mean()),
                             expected_sign=sign, t=t, p=p, mean_abs_corr=dens, null_mean=nul.mean(),
                             preservation_z=(dens - nul.mean()) / nul.std(), preservation_p=(np.sum(nul >= dens) + 1) / 1001))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(O / "network_per_study.tsv", sep="\t", index=False)
    # 3 rebuilds
    rb, hs = [], []

    def compare(tag, zz):
        new, kme, beta = build(zz, genes)
        for m, gl in mods.items():
            best = max(new, key=lambda l: len(set(new[l]) & set(gl)) / len(set(new[l]) | set(gl)))
            J = len(set(new[best]) & set(gl)) / len(set(new[best]) | set(gl))
            frac_in = len(set(new[best]) & set(gl)) / len(gl)
            k = kme[best]
            thr = k.quantile(0.9)
            kept = [h for h in top[m] if h in k.index and k[h] >= thr]
            rb.append(dict(rebuild=tag, power=beta, module=m, matched_size=len(new[best]), jaccard=J,
                           fraction_of_original_in_match=frac_in, top25_hubs_in_top10pct=len(kept) / 25))
            for h in top[m]:
                hs.append(dict(rebuild=tag, module=m, gene=h, in_match=h in k.index, top10=h in kept))
        print(tag, rb[-2:], flush=True)

    for st in MA.DISCOVERY:
        keep = s.index[s.study != st]
        compare(f"without {st}", zexpr(tpm.loc[genes], s.loc[keep]).dropna())
    for b in range(20):
        idx = []
        for _, g in s.groupby(["study", "condition"]):
            idx += list(RNG.choice(g.index, len(g), replace=True))
        zz = zexpr(tpm.loc[genes], s.loc[s.index.isin(idx)])   # unique libraries of the bootstrap draw
        # duplicate the resampled libraries so that multiplicity is respected
        cols = []
        for i, c in enumerate(idx):
            col = zz[c].rename(f"{c}_{i}")
            cols.append(col)
        zb = pd.concat(cols, axis=1).dropna()
        compare(f"bootstrap {b + 1}", zb)
    pd.DataFrame(rb).to_csv(O / "network_rebuilds.tsv", sep="\t", index=False)
    H = pd.DataFrame(hs)
    H.groupby(["module", "gene"]).top10.mean().rename("frequency_top10pct").reset_index().to_csv(
        O / "network_hub_stability.tsv", sep="\t", index=False)
    R_ = pd.DataFrame(rb)
    print(R_.groupby(["module", R_.rebuild.str.startswith("bootstrap")])[["jaccard", "fraction_of_original_in_match",
                                                                         "top25_hubs_in_top10pct"]].describe().round(2).T)


if __name__ == "__main__":
    main()
