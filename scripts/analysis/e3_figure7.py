"""Figure 7: public Arabidopsis experimental data for orthologues of Core Drought Signature genes."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import a2_meta as MA
from a6_figures import R, UP, DOWN, INK, INK2, GRID, panel, save
import e3_public_functional as PF

LIGHT_UP, LIGHT_DN, BG = "#f4a9a8", "#9cc3ef", "#c3c2b7"


def bars(ax, S, gset, props, labels, dark, light, title):
    s = S[(S.orthology == "best hit") & (S.gene_set == gset)].set_index("property").loc[props]
    y = np.arange(len(props))[::-1]
    h = 0.26
    ax.barh(y + h, s.pct_set, h, color=dark, label="Signature orthologues")
    ax.barh(y, s.pct_nonreproducible, h, color=light, label="Non-reproducible drought-responsive")
    ax.barh(y - h, s.pct_background, h, color=BG, label="All tested genes")
    for yy, v, p in zip(y, s.pct_set, s.P_vs_nonreproducible):
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        ax.text(v + 1.5, yy + h, f"{v:.0f}% {star}", va="center", fontsize=6.5, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(0, 105)
    ax.set_xlabel("Orthologues with property (%)")
    ax.xaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title, loc="left")


def main():
    S = pd.read_csv(R / "extended/public_functional_sets.tsv", sep="\t")
    expr = ["osmotic-responsive in WT (GSE110079)", "ABA-dependent in aba2-2 (GSE110079)",
            "ABA-dependent in aba3-1 (GSE110079)", "ABA-responsive (GSE80565)"]
    elab = ["Osmotic (PEG) response\nin wild type", "Response reduced\nin aba2-2", "Response reduced\nin aba3-1",
            "Responsive to\nABA treatment"]
    chip = ["bound by ABF1/3/4 (GSE80564)", "bound by GBF3 (GSE80564)", "bound by HB7 (GSE80564)",
            "bound by NFYC2 (GSE80564)", "bound by RD26 (GSE80564)", "bound by DREB2A (GSE80564)"]
    clab = ["ABF1/3/4", "GBF3", "HB7", "NF-YC2", "RD26", "DREB2A"]
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.12, top=0.9)
    gs2 = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.62, top=0.9)
    ax = fig.add_subplot(gs[0, 0])
    bars(ax, S, "CDS induced", expr, elab, UP, LIGHT_UP, "Induced signature genes")
    panel(ax, "a")
    hs, ls = ax.get_legend_handles_labels()
    fig.legend(hs, ["Signature orthologues", "Orthologues of non-reproducible drought-responsive genes",
                    "Orthologues of all tested genes"], loc="upper center", ncol=3, fontsize=7,
               bbox_to_anchor=(0.5, 1.0))
    ax = fig.add_subplot(gs[0, 1])
    bars(ax, S, "CDS repressed", expr, elab, DOWN, LIGHT_DN, "Repressed signature genes")
    ax.set_yticklabels([])
    panel(ax, "b")
    ax = fig.add_subplot(gs2[1, 0])
    bars(ax, S, "CDS induced", chip, clab, UP, LIGHT_UP, "Promoter binding after ABA (ChIP-seq)")
    panel(ax, "c")
    # d: retention of the osmotic response in aba2-2
    ax = fig.add_subplot(gs2[1, 1])
    t = PF.gse110079()
    dep = PF.aba_dependence(t, +1)
    orth = pd.read_csv(R / "function/arabidopsis_orthologs.tsv", sep="\t")
    w2a = dict(zip(orth.watermelon, orth.AGI))
    meta = pd.read_csv(R / "meta.tsv.gz", sep="\t", index_col=0)
    cds = pd.read_csv(R / "cds.tsv", sep="\t", index_col=0)
    anyup = set()
    for st in MA.DISCOVERY:
        d = pd.read_csv(R / "de" / f"{st}.tsv", sep="\t", index_col=0)
        anyup |= set(d.index[(d.padj < 0.05) & (d.log2FoldChange > 0)])
    A = {w2a[g] for g in cds.index[cds.direction == "up"] if g in w2a}
    C = {w2a[g] for g in anyup - set(cds.index) if g in w2a} - A
    B = {w2a[g] for g in meta.index if g in w2a} - A
    data = [dep.fraction_retained_aba2.reindex(sorted(x)).dropna().clip(-0.5, 1.5) * 100 for x in (A, C, B)]
    bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False,
                    medianprops=dict(color=INK, lw=1.2), whiskerprops=dict(color=INK2), capprops=dict(color=INK2))
    for p, c in zip(bp["boxes"], [UP, LIGHT_UP, BG]):
        p.set_facecolor(c)
        p.set_edgecolor(INK2)
    ax.axhline(100, color=INK2, lw=0.6, ls="--")
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels([f"Induced\nsignature\n(n = {len(data[0])})", f"Non-\nreproducible\n(n = {len(data[1])})",
                        f"All tested\n(n = {len(data[2])})"], fontsize=7)
    ax.set_ylabel("Wild-type PEG induction\nretained in aba2-2 (%)")
    ax.set_title("ABA dependence of PEG induction", loc="left")
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    panel(ax, "d")
    save(fig, "Figure7")
    print([np.median(x) for x in data])


if __name__ == "__main__":
    main()
