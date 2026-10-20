"""Co-expression network (WGCNA-style) on the 43 watermelon drought discovery samples.

Expression: log2(TPM+1), z-scored within study (removes study/tissue baselines), top 5,000 genes by
variance across samples. Signed network, soft threshold chosen as the lowest power with scale-free
R^2 >= 0.8 (or the best available), topological overlap, average-linkage clustering, dynamic tree cut
(hybrid), module merging at eigengene correlation 0.75. Module-trait correlation with drought;
intramodular connectivity (kME) defines hub genes.
Outputs in results/network/.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import average, fcluster
from scipy.spatial.distance import squareform

import a2_meta as MA

M = Path(os.environ.get("PROJECT_DIR", ".") + "")
R = M / "results"
N = R / "network"
N.mkdir(exist_ok=True)


def scale_free_fit(adj_k, nbins=10):
    k = adj_k[adj_k > 0]
    h, edges = np.histogram(k, bins=nbins)
    mid = (edges[1:] + edges[:-1]) / 2
    ok = h > 0
    x, y = np.log10(mid[ok]), np.log10(h[ok] / h.sum())
    slope, inter, r, *_ = stats.linregress(x, y)
    return -np.sign(slope) * r ** 2, slope


def eigengene(x):
    xc = x - x.mean(axis=1, keepdims=True)
    u, s, vt = np.linalg.svd(xc, full_matrices=False)
    me = vt[0]
    if np.corrcoef(me, x.mean(axis=0))[0, 1] < 0:
        me = -me
    return me


def main():
    samples = pd.read_csv(R / "samples.tsv", sep="\t", index_col=0, dtype=str).fillna("")
    s = samples[samples.study.isin(MA.DISCOVERY)]
    tpm = pd.read_csv(R / "tpm_watermelon.tsv.gz", sep="\t", index_col=0)[s.index]
    x = np.log2(tpm + 1)
    expressed = (tpm >= 1).mean(axis=1) >= 0.5
    x = x[expressed]
    z = []
    for _, ss in s.groupby("study"):
        sub = x[ss.index]
        sd = sub.std(axis=1).replace(0, np.nan)
        z.append(sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0))
    z = pd.concat(z, axis=1)[s.index].dropna()
    # Gene selection: within-study variance of log2 expression (study means removed, NOT standardised; the variance of
    # within-study z-scores is ~constant across genes and therefore cannot rank genes).
    xc = pd.concat([x.loc[z.index, ss.index].sub(x.loc[z.index, ss.index].mean(axis=1), axis=0)
                    for _, ss in s.groupby("study")], axis=1)
    z = z.loc[xc.var(axis=1).sort_values(ascending=False).index[:5000]]
    X = z.values
    C = np.corrcoef(X)
    S = (1 + C) / 2                                           # signed
    fits = []
    for b in [4, 6, 8, 10, 12, 14, 16, 18, 20]:
        A = S ** b
        np.fill_diagonal(A, 0)
        r2, slope = scale_free_fit(A.sum(axis=1))
        fits.append(dict(power=b, R2=r2, slope=slope, mean_k=A.sum(axis=1).mean()))
    fits = pd.DataFrame(fits)
    fits.to_csv(N / "soft_threshold.tsv", sep="\t", index=False)
    ok = fits[fits.R2 >= 0.8]
    beta = int(ok.power.iloc[0]) if len(ok) else int(fits.loc[fits.R2.idxmax(), "power"])
    print(f"  soft threshold power {beta} (R2={fits.set_index('power').R2[beta]:.2f})", flush=True)
    A = S ** beta
    np.fill_diagonal(A, 0)
    k = A.sum(axis=1)
    L = A @ A
    TOM = (L + A) / (np.minimum.outer(k, k) + 1 - A)
    np.fill_diagonal(TOM, 1)
    D = 1 - TOM
    Z = average(squareform(D, checks=False))
    try:
        from dynamicTreeCut import cutreeHybrid
        cut = cutreeHybrid(Z, squareform(D, checks=False), minClusterSize=30, deepSplit=2, pamStage=False, verbose=0)
        labels = np.asarray(cut["labels"])
        print(f"  dynamic hybrid tree cut: {len(set(labels)) - (0 in labels)} modules", flush=True)
    except Exception as e:
        print("  dynamic tree cut failed:", repr(e)[:200], flush=True)
        # static tree cut (WGCNA cutreeStatic): choose the height giving the most genes in modules >= 30
        best = None
        for q in np.linspace(0.50, 0.995, 100):
            h = np.quantile(Z[:, 2], q)
            lab = fcluster(Z, t=h, criterion="distance")
            sizes = pd.Series(lab).value_counts()
            n_mod = int((sizes >= 30).sum())
            cover = sizes[sizes >= 30].sum() / len(lab)
            if cover < 0.6:
                continue
            if best is None or n_mod > best[2]:
                best = (sizes[sizes >= 30].sum(), lab, n_mod, h)
        labels = best[1]
        print(f"  static cut height {best[3]:.3f}: {best[2]} modules covering {best[0]} genes", flush=True)
    genes = z.index
    y = (s.condition == "stress").astype(int).values
    mods = {}
    for lab in sorted(set(labels)):
        if lab == 0:
            continue
        idx = np.where(labels == lab)[0]
        if len(idx) < 30:
            continue
        mods[lab] = idx
    # merge modules with highly correlated eigengenes
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
    cds = pd.read_csv(R / "cds.tsv", sep="\t", index_col=0)
    rows, gene_rows = [], []
    for i, (l, idx) in enumerate(sorted(mods.items(), key=lambda t: -len(t[1])), 1):
        me = MEs[l]
        r, p = stats.pointbiserialr(y, me)
        name = f"M{i}"
        g = genes[idx]
        kme = np.array([np.corrcoef(X[j], me)[0, 1] for j in idx])
        n_cds = len(set(g) & set(cds.index))
        hyp_p = stats.hypergeom.sf(n_cds - 1, len(genes), len(set(cds.index) & set(genes)), len(g))
        rows.append(dict(module=name, size=len(idx), r_drought=r, p_drought=p, n_cds=n_cds, cds_enrich_p=hyp_p))
        for gg, kk in zip(g, kme):
            gene_rows.append(dict(gene=gg, module=name, kME=kk, in_cds=gg in cds.index))
    mt = pd.DataFrame(rows).sort_values("p_drought")
    gm = pd.DataFrame(gene_rows)
    mt.to_csv(N / "modules.tsv", sep="\t", index=False)
    gm.to_csv(N / "gene_modules.tsv", sep="\t", index=False)
    print(mt.round(4).to_string(index=False), flush=True)
    hubs = gm.merge(mt[["module", "r_drought", "p_drought"]], on="module")
    hubs = hubs[hubs.p_drought < 0.001].sort_values("kME", ascending=False).groupby("module").head(25)
    hubs.to_csv(N / "hub_genes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
