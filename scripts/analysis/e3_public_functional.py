"""public experimental evidence for the regulatory interpretation of the Core Drought
Signature (CDS), using only existing public data.

Watermelon genes are mapped to Arabidopsis by the DIAMOND best hits of a4_function.py (RBH-only as sensitivity).
For each watermelon gene set (CDS induced / repressed, co-expression modules M3 / M1) the unique Arabidopsis
orthologues are compared with the orthologues of all other genes tested in the meta-analysis (Fisher exact test)
and, as a stricter comparison, with orthologues of genes that changed in the same direction in at least one
discovery study (padj < 0.05) but are not in the CDS ("non-reproducible drought-responsive genes").

Resources (all public, GEO):
 1. GSE110079 (Watanabe et al., 2018): Col-0, ABA-deficient aba2-2 and aba3-1 seedlings, mock vs PEG,
    3 h (solution) and 8 h (plate), three biological replicates, Agilent arrays, RMA-normalised.
    Per gene and time point: PEG effect in each genotype and genotype x treatment interaction (OLS, pooled
    residual variance of the 18 arrays, df = 12; BH-FDR).
      * osmotic-responsive in WT: FDR < 0.05 in the set's direction at either time point;
      * ABA-dependent: responsive in WT and the response significantly smaller in aba2-2 (interaction
        FDR < 0.05 at a time point where the gene is responsive); replication in aba3-1 reported separately;
      * fraction of the WT response retained in aba2-2 (mutant effect / WT effect).
 2. GSE80565 (Song et al., 2016): ABA vs mock RNA-seq time course (1-60 h) of Col-0 seedlings;
    ABA-responsive = FDR < 0.05 in the set's direction at any time point.
 3. GSE80564 (Song et al., 2016): ChIP-seq of 21 ABA-responsive transcription factors after ABA treatment;
    a gene is bound if a peak overlaps -1,000 to +500 bp of the transcription start (primary) or
    -1,000 bp to the gene end (sensitivity).
Outputs: results/extended/public_functional_sets.tsv, public_functional_genes.tsv
"""
import os
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import a2_meta as MA

R = Path(os.environ.get("PROJECT_DIR", ".") + "/results")
P = Path(os.environ.get("PROJECT_DIR", ".") + "/pub")
REF = Path(os.environ.get("PROJECT_DIR", ".") + "/ref")
O = R / "extended"


def bh(p):
    return MA.bh(np.asarray(p, float))


# ---------------------------------------------------------------- 1. ABA-deficient mutants under PEG
def gse110079():
    lines = gzip.open(P / "gse110079/matrix.txt.gz", "rt").read().splitlines()
    titles = [l for l in lines if l.startswith("!Sample_title")][0].split("\t")[1:]
    titles = [t.strip('"') for t in titles]
    tab = [l.split("\t") for l in lines if not l.startswith("!") and l]
    df = pd.DataFrame([r[1:] for r in tab[1:]], index=[r[0].strip('"') for r in tab[1:]]).astype(float)
    df.columns = titles
    df = df[df.index.str.fullmatch(r"AT[1-5CM]G\d{5}")]
    df = df[df.mean(axis=1) > np.percentile(df.mean(axis=1), 20)]           # drop the lowest-intensity 20%
    out = {}
    for tp, tag in [("3 h", "3h"), ("8 h", "8h")]:
        cols = [c for c in df.columns if f"_{tp} " in c]
        x = df[cols]
        grp = {}
        for g in ["Col-0", "aba2-2", "aba3-1"]:
            for t in ["mock", "PEG"]:
                grp[(g, t)] = [c for c in cols if c.startswith(g) and f" {t} " in c]
        assert all(len(v) == 3 for v in grp.values()), grp
        m = {k: x[v].mean(axis=1) for k, v in grp.items()}
        rss = sum(((x[v].sub(m[k], axis=0)) ** 2).sum(axis=1) for k, v in grp.items())
        df_res = len(cols) - len(grp)
        s = np.sqrt(rss / df_res)
        se_d, se_i = s * np.sqrt(2 / 3), s * np.sqrt(4 / 3)
        eff = {g: m[(g, "PEG")] - m[(g, "mock")] for g in ["Col-0", "aba2-2", "aba3-1"]}
        r = pd.DataFrame(index=x.index)
        for g in eff:
            r[f"{g}_lfc_{tag}"] = eff[g]
            r[f"{g}_fdr_{tag}"] = bh(2 * stats.t.sf(np.abs(eff[g] / se_d), df_res))
        for g in ["aba2-2", "aba3-1"]:
            d = eff["Col-0"] - eff[g]
            r[f"int_{g}_{tag}"] = d
            r[f"int_{g}_fdr_{tag}"] = bh(2 * stats.t.sf(np.abs(d / se_i), df_res))
        out[tag] = r
    return pd.concat(out.values(), axis=1)


def aba_dependence(t, sign):
    """Per-gene flags for one direction (sign = +1 induced, -1 repressed)."""
    res = pd.DataFrame(index=t.index)
    resp, dep, dep3, ret = [], [], [], []
    for tag in ["3h", "8h"]:
        wt = t[f"Col-0_lfc_{tag}"] * sign
        r_ = (wt > 0) & (t[f"Col-0_fdr_{tag}"] < 0.05)
        d2 = r_ & (t[f"int_aba2-2_{tag}"] * sign > 0) & (t[f"int_aba2-2_fdr_{tag}"] < 0.05)
        d3 = r_ & (t[f"int_aba3-1_{tag}"] * sign > 0) & (t[f"int_aba3-1_fdr_{tag}"] < 0.05)
        resp.append(r_)
        dep.append(d2)
        dep3.append(d3)
        ret.append((t[f"aba2-2_lfc_{tag}"] * sign / wt).where(r_))
    res["osmotic_responsive"] = pd.concat(resp, axis=1).any(axis=1)
    res["aba2_dependent"] = pd.concat(dep, axis=1).any(axis=1)
    res["aba3_dependent"] = pd.concat(dep3, axis=1).any(axis=1)
    res["fraction_retained_aba2"] = pd.concat(ret, axis=1).mean(axis=1).clip(-1, 2)
    return res


# ---------------------------------------------------------------- 2. ABA time course
def gse80565():
    t = pd.read_csv(P / "song2016/GSE80565_RNAseq_ABAtimeSeries_rpkm_log2FC_FDR_filterCPM2n2.txt.gz", sep=r"[\t ]",
                    engine="python", index_col=0)    # first column is space-delimited in the GEO file
    tps = [c.split("_")[0].replace("AvsE", "") for c in t.columns if c.endswith("_logFC")]
    up = pd.concat([(t[f"AvsE{h}_logFC"] > 0) & (t[f"AvsE{h}_FDR"] < 0.05) for h in tps], axis=1).any(axis=1)
    dn = pd.concat([(t[f"AvsE{h}_logFC"] < 0) & (t[f"AvsE{h}_FDR"] < 0.05) for h in tps], axis=1).any(axis=1)
    return pd.DataFrame({"ABA_induced": up, "ABA_repressed": dn})


# ---------------------------------------------------------------- 3. ChIP-seq of 21 ABA-responsive TFs
def gse80564(window="promoter"):
    """window 'promoter': -1,000 to +500 bp of the TSS (primary); 'gene': -1,000 bp to the gene end."""
    g = pd.read_csv(REF / "at_genes.bed", sep="\t", header=None, names=["chr", "s", "e", "strand", "agi"])
    g = g[g.chr.isin(list("12345"))]
    if window == "promoter":
        g["lo"] = np.where(g.strand == "+", g.s - 1000, g.e - 500)
        g["hi"] = np.where(g.strand == "+", g.s + 500, g.e + 1000)
    else:
        g["lo"] = np.where(g.strand == "+", g.s - 1000, g.s)
        g["hi"] = np.where(g.strand == "+", g.e, g.e + 1000)
    bound = {}
    for f in sorted((P / "song2016").glob("GSE80564_*_ABA_optimal_narrowPeak_p16.bed.gz")):
        parts = f.name.split("_")
        tf = parts[2] if parts[2] != "ABA" else parts[1]
        pk = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "s", "e"])
        pk["chr"] = pk.chr.str.replace("chr", "")
        hits = set()
        for ch, gg in g.groupby("chr"):
            p = pk[pk.chr == ch].sort_values("s")
            if p.empty:
                continue
            ps, pe = p.s.values, p.e.values
            maxlen = (pe - ps).max()
            for lo, hi, a in zip(gg.lo.values, gg.hi.values, gg.agi.values):
                i = np.searchsorted(ps, hi, side="right")      # peaks starting before the gene end
                j = np.searchsorted(ps, lo - maxlen, side="left")
                if np.any(pe[j:i] >= lo):
                    hits.add(a)
        bound[tf] = hits
        print(tf, len(pk), "peaks", len(hits), "genes", flush=True)
    B = pd.DataFrame({tf: g.agi.isin(h).values for tf, h in bound.items()}, index=g.agi.values)
    B["ABF_any"] = B[[c for c in B if c in ("ABF1", "ABF3", "ABF4")]].any(axis=1)
    B["n_TFs_bound"] = B[[c for c in bound]].sum(axis=1)
    return B, list(bound)


def fisher(test, bg):
    """test, bg: boolean Series (unique AGIs)."""
    a, b = int(test.sum()), int((~test).sum())
    c, d = int(bg.sum()), int((~bg).sum())
    orr, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return a, a + b, a / max(1, a + b), c / max(1, c + d), orr, p


def main():
    orth = pd.read_csv(R / "function/arabidopsis_orthologs.tsv", sep="\t")
    meta = pd.read_csv(R / "meta.tsv.gz", sep="\t", index_col=0)
    cds = pd.read_csv(R / "cds.tsv", sep="\t", index_col=0)
    gm = pd.read_csv(R / "network/gene_modules.tsv", sep="\t")
    mt = pd.read_csv(R / "network/modules.tsv", sep="\t")
    IND, REP = mt.loc[mt.r_drought.idxmax(), "module"], mt.loc[mt.r_drought.idxmin(), "module"]
    # non-reproducible drought-responsive genes: significant in >= 1 discovery study, same direction, not CDS
    any_up, any_dn = set(), set()
    for st in MA.DISCOVERY:
        d = pd.read_csv(R / "de" / f"{st}.tsv", sep="\t", index_col=0)
        any_up |= set(d.index[(d.padj < 0.05) & (d.log2FoldChange > 0)])
        any_dn |= set(d.index[(d.padj < 0.05) & (d.log2FoldChange < 0)])
    sets = {"CDS induced": (set(cds.index[cds.direction == "up"]), +1),
            "CDS repressed": (set(cds.index[cds.direction == "down"]), -1),
            f"{IND} (induced module)": (set(gm.gene[gm.module == IND]), +1),
            f"{REP} (repressed module)": (set(gm.gene[gm.module == REP]), -1)}
    comp = {+1: any_up - set(cds.index), -1: any_dn - set(cds.index)}

    t110 = gse110079()
    dep = {s: aba_dependence(t110, s) for s in (+1, -1)}
    aba = gse80565()
    chip, tfs = gse80564("promoter")
    chip_g, _ = gse80564("gene")

    rows, genes = [], []
    for variant, o in [("best hit", orth), ("RBH only", orth[orth.RBH])]:
        w2a = dict(zip(o.watermelon, o.AGI))
        tested = {w2a[g] for g in meta.index if g in w2a}
        for name, (ws, sign) in sets.items():
            A = {w2a[g] for g in ws if g in w2a}
            C = {w2a[g] for g in comp[sign] if g in w2a} - A
            props = []
            dd = dep[sign]
            props.append(("osmotic-responsive in WT (GSE110079)", dd.osmotic_responsive))
            props.append(("ABA-dependent in aba2-2 (GSE110079)", dd.aba2_dependent))
            props.append(("ABA-dependent in aba3-1 (GSE110079)", dd.aba3_dependent))
            props.append(("ABA-dependent in aba2-2 among WT-responsive (GSE110079)",
                          dd.aba2_dependent[dd.osmotic_responsive]))
            props.append(("ABA-responsive (GSE80565)", aba.ABA_induced if sign > 0 else aba.ABA_repressed))
            props.append(("bound by ABF1/3/4 (GSE80564)", chip.ABF_any))
            props.append(("bound by ABF1/3/4, gene window (GSE80564)", chip_g.ABF_any))
            for tf in ["GBF3", "HB7", "NFYB2", "NFYC2", "RD26", "DREB2A"]:
                if tf in chip:
                    props.append((f"bound by {tf} (GSE80564)", chip[tf]))
            for pname, s in props:
                univ = set(s.index)
                a_ = pd.Series(s.reindex(sorted(A & univ)).values)
                bg = pd.Series(s.reindex(sorted((tested & univ) - A)).values)
                cc = pd.Series(s.reindex(sorted(C & univ)).values)
                k, n, fa, fb, orr, p = fisher(a_, bg)
                _, _, _, fc, orc, pc = fisher(a_, cc)
                rows.append(dict(orthology=variant, gene_set=name, property=pname, n_orthologues=n, n_with=k,
                                 pct_set=100 * fa, pct_background=100 * fb, OR_vs_background=orr, P_vs_background=p,
                                 pct_nonreproducible=100 * fc, n_nonreproducible=len(cc),
                                 OR_vs_nonreproducible=orc, P_vs_nonreproducible=pc))
            # quantitative retention of the osmotic response in aba2-2
            ret_a = dep[sign].fraction_retained_aba2.reindex(sorted(A)).dropna()
            ret_c = dep[sign].fraction_retained_aba2.reindex(sorted(C)).dropna()
            ret_b = dep[sign].fraction_retained_aba2.reindex(sorted(tested - A)).dropna()
            u = stats.mannwhitneyu(ret_a, ret_b, alternative="less")
            uc = stats.mannwhitneyu(ret_a, ret_c, alternative="less")
            rows.append(dict(orthology=variant, gene_set=name, property="median fraction of WT osmotic response retained in aba2-2",
                             n_orthologues=len(ret_a), pct_set=100 * ret_a.median(), pct_background=100 * ret_b.median(),
                             P_vs_background=u.pvalue, pct_nonreproducible=100 * ret_c.median(),
                             n_nonreproducible=len(ret_c), P_vs_nonreproducible=uc.pvalue))
            if variant == "best hit":
                for g in sorted(ws):
                    if g not in w2a:
                        continue
                    a = w2a[g]
                    rec = dict(gene=g, gene_set=name, AGI=a)
                    if a in dd.index:
                        rec.update(osmotic_responsive_WT=bool(dd.osmotic_responsive[a]),
                                   aba2_dependent=bool(dd.aba2_dependent[a]), aba3_dependent=bool(dd.aba3_dependent[a]),
                                   fraction_retained_aba2=dd.fraction_retained_aba2[a])
                    if a in aba.index:
                        rec["ABA_responsive"] = bool((aba.ABA_induced if sign > 0 else aba.ABA_repressed)[a])
                    if a in chip.index:
                        rec.update(ABF_bound=bool(chip.ABF_any[a]), n_TFs_bound=int(chip.n_TFs_bound[a]),
                                   TFs_bound=",".join(tf for tf in tfs if chip.loc[a, tf]))
                    genes.append(rec)
    S = pd.DataFrame(rows)
    S.to_csv(O / "public_functional_sets.tsv", sep="\t", index=False)
    pd.DataFrame(genes).to_csv(O / "public_functional_genes.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(S[S.orthology == "best hit"].drop(columns="orthology").round(4).to_string())
    print(S[S.orthology == "RBH only"].drop(columns="orthology")[["gene_set", "property", "n_orthologues", "pct_set", "pct_background", "OR_vs_background", "P_vs_background", "pct_nonreproducible", "P_vs_nonreproducible"]].round(4).to_string())


if __name__ == "__main__":
    main()
