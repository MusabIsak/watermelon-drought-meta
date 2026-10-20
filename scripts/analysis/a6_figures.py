"""Publication figures (PNG 300 dpi + PDF) in results/figures/.

Palette (validated reference instance): diverging blue #2a78d6 <-> red #e34948 with grey midpoint;
categorical slots in fixed order blue, orange, aqua, yellow; ink #0b0b0b / #52514e; grid #e1e0d9.
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

M = Path(os.environ.get("PROJECT_DIR", ".") + "")
R = M / "results"
FIG = R / "figures"
FIG.mkdir(exist_ok=True)

UP, DOWN, GREY = "#e34948", "#2a78d6", "#898781"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID, AXIS = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7"
DIV = LinearSegmentedColormap.from_list("div", [DOWN, "#f0efec", UP])
LABEL = {"WM_Leaf_M20Y34_2018": "Leaf, M20/Y34, 4–8 d", "WM_Leaf_CS_2020": "Leaf, Crimson Sweet, 8 d",
         "WM_Root_PEG_2016": "Root, M08, PEG 6 h", "WM_Leaf_PEGtime_2024": "Leaf, PEG 1–24 h",
         "CU_Leaf_PEG_2025": "Cucumber leaf, PEG", "CU_Leaf_PEG48h_2024": "Cucumber, PEG 48 h (secondary)",
         "CU_Leaf_graft_2022": "Cucumber, grafted, WS",
         "WM_Leaf_salt7h_2020": "Salt, leaf 7 h", "WM_salt_LR_2025__leaf": "Salt, leaf",
         "WM_salt_LR_2025__root": "Salt, root", "WM_Leaf_cold_2017": "Cold, leaf 4 °C"}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": AXIS, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 9, "axes.titleweight": "bold",
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
                     "savefig.dpi": 300, "figure.dpi": 150, "pdf.fonttype": 42, "ps.fonttype": 42})


def panel(ax, letter):
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", color=INK, va="bottom")


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def fig2():
    con = pd.read_csv(R / "contrasts.tsv", sep="\t")
    meta = pd.read_csv(R / "meta.tsv.gz", sep="\t", index_col=0)
    cds = pd.read_csv(R / "cds.tsv", sep="\t", index_col=0)
    disc = ["WM_Leaf_M20Y34_2018", "WM_Leaf_CS_2020", "WM_Root_PEG_2016", "WM_Leaf_PEGtime_2024"]
    fig = plt.figure(figsize=(7.2, 6.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.25], width_ratios=[1, 1.15], hspace=0.45, wspace=0.35)
    # a: DEGs per discovery study
    ax = fig.add_subplot(gs[0, 0])
    c = con.set_index("contrast").loc[disc]
    y = np.arange(len(disc))[::-1]
    ax.barh(y + 0.18, c.up, 0.34, color=UP, label="Up")
    ax.barh(y - 0.18, c.down, 0.34, color=DOWN, label="Down")
    ax.set_yticks(y, [LABEL[d] for d in disc])
    ax.set_xlabel("DEGs (padj < 0.05, |log2FC| ≥ 1)")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right")
    ax.set_title("DEGs per study", loc="left")
    panel(ax, "a")
    # b: overlap -> consistency histogram
    ax = fig.add_subplot(gs[0, 1])
    sig = {}
    for d in disc:
        de = pd.read_csv(R / "de" / f"{d}.tsv", sep="\t", index_col=0)
        sig[d] = set(de.index[(de.padj < 0.05) & (de.log2FoldChange.abs() >= 1)])
    allg = set().union(*sig.values())
    n_in = pd.Series({g: sum(g in s for s in sig.values()) for g in allg}).value_counts().sort_index()
    ax.bar(n_in.index, n_in.values, color=[GREY, GREY, "#5598e7", INK][:len(n_in)], width=0.7)
    for xv, yv in zip(n_in.index, n_in.values):
        ax.text(xv, yv, f"{yv:,}", ha="center", va="bottom", fontsize=7, color=INK2)
    ax.set_xticks(range(1, 5))
    ax.set_xlabel("Number of studies in which the gene is a DEG")
    ax.set_ylabel("Genes")
    ax.set_title("Few genes respond in every study", loc="left")
    panel(ax, "b")
    # c: meta volcano
    ax = fig.add_subplot(gs[1, 0])
    yv = -np.log10(meta.meta_fdr.clip(lower=1e-30))
    ax.scatter(meta.meta_log2FC, yv, s=2, color=GRID, rasterized=True)
    for d, col in [("up", UP), ("down", DOWN)]:
        s = cds[cds.direction == d]
        ax.scatter(s.meta_log2FC, -np.log10(s.meta_fdr.clip(lower=1e-30)), s=5, color=col, rasterized=True,
                   label=f"CDS {d} ({len(s)})")
    ax.axhline(2, color=AXIS, lw=0.6, ls="--")
    ax.axvline(1, color=AXIS, lw=0.6, ls="--")
    ax.axvline(-1, color=AXIS, lw=0.6, ls="--")
    ax.set_xlabel("Meta-analysis log2 fold change (drought vs control)")
    ax.set_ylabel("−log10 meta FDR")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, markerscale=3, fontsize=6.5, borderaxespad=0)
    ax.set_title("Core Drought Signature", loc="left", pad=16)
    panel(ax, "c")
    # d: heatmap of CDS log2FC across discovery + specificity + cucumber
    ax = fig.add_subplot(gs[1, 1])
    cols = disc + [x for x in con[con.role == "specificity"].contrast]
    mat = pd.DataFrame({x: pd.read_csv(R / "de" / f"{x}.tsv", sep="\t", index_col=0).log2FoldChange for x in cols})
    rbh = pd.read_csv(M / "ref" / "rbh_wm_cu.tsv", sep="\t")
    o = dict(zip(rbh.watermelon, rbh.cucumber))
    cu = [x for x in con[con.crop == "cucumber"].contrast]
    for x in cu:
        d = pd.read_csv(R / "de" / f"{x}.tsv", sep="\t", index_col=0).log2FoldChange
        mat[x] = pd.Series({g: d.get(o.get(g)) for g in cds.index}, dtype=float)
    mat = mat.reindex(cds.sort_values("meta_log2FC", ascending=False).index)
    im = ax.imshow(mat.values, aspect="auto", cmap=DIV, vmin=-4, vmax=4, interpolation="nearest")
    SHORT = {"WM_Leaf_M20Y34_2018": "WM leaf M20/Y34", "WM_Leaf_CS_2020": "WM leaf CS", "WM_Root_PEG_2016": "WM root PEG",
             "WM_Leaf_PEGtime_2024": "WM leaf PEG 1-24 h", "WM_Leaf_cold_2017": "WM cold", "WM_Leaf_salt7h_2020": "WM salt 7 h",
             "WM_salt_LR_2025__leaf": "WM salt leaf", "WM_salt_LR_2025__root": "WM salt root", "CU_Leaf_PEG48h_2024": "CU PEG 48 h (2nd)",
             "CU_Leaf_PEG_2025": "CU PEG", "CU_Leaf_graft_2022": "CU grafted WS"}
    ax.set_xticks(range(mat.shape[1]), [SHORT.get(x, x) for x in mat.columns], rotation=55, ha="right", fontsize=7.5)
    ax.set_yticks([])
    ax.set_ylabel(f"CDS genes (n = {len(cds)}; induced on top)")
    for xpos in [len(disc) - 0.5, len(cols) - 0.5]:
        ax.axvline(xpos, color=INK, lw=0.8)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("log2FC", fontsize=7)
    for xpos, lab in [((len(disc) - 1) / 2, "Drought\n(discovery)"), ((len(disc) + len(cols) - 1) / 2, "Salt / cold"),
                      ((len(cols) + mat.shape[1] - 1) / 2, "Cucumber\ndrought")]:
        ax.text(xpos, 1.01, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.5, color=INK2)
    ax.set_title("CDS genes across all studies", loc="left", pad=22)
    panel(ax, "d")
    save(fig, "Fig2_meta_analysis")


def fig3():
    """concordance (primary) with dependence-aware evidence; AUC shown against random-gene-set nulls."""
    lo = pd.read_csv(R / "validation" / "loso.tsv", sep="\t")
    vu = pd.read_csv(R / "extended" / "validation_uncertainty.tsv", sep="\t").set_index("study")
    cs = pd.read_csv(R / "validation" / "cross_species.tsv", sep="\t").set_index("study")
    sc = pd.read_csv(R / "validation" / "cross_species_scores.tsv", sep="\t", index_col=0)
    order = ["WM_Leaf_M20Y34_2018", "WM_Leaf_CS_2020", "WM_Root_PEG_2016", "WM_Leaf_PEGtime_2024",
             "CU_Leaf_PEG_2025", "CU_Leaf_PEG48h_2024", "CU_Leaf_graft_2022"]
    short = {"WM_Leaf_M20Y34_2018": "WM leaf\nM20/Y34", "WM_Leaf_CS_2020": "WM leaf\nCS", "WM_Root_PEG_2016": "WM root\nPEG 6 h",
             "WM_Leaf_PEGtime_2024": "WM leaf\nPEG 1-24 h", "CU_Leaf_PEG_2025": "CU leaf\nPEG", "CU_Leaf_PEG48h_2024": "CU leaf PEG 48 h\n(secondary)",
             "CU_Leaf_graft_2022": "CU grafted\nWS"}
    # significant-and-concordant fractions
    sig = {r.held_out: r.frac_sig_same_dir for r in lo[lo.model == "CDS score"].itertuples()}
    rbh = pd.read_csv(M / "ref" / "rbh_wm_cu.tsv", sep="\t")
    o = dict(zip(rbh.watermelon, rbh.cucumber))
    cds = pd.read_csv(R / "cds.tsv", sep="\t", index_col=0)
    conc = {}
    for st in order:
        if st.startswith("CU"):
            de = pd.read_csv(R / "de" / f"{st}.tsv", sep="\t", index_col=0)
            v = [(np.sign(de.at[o[g], "log2FoldChange"]) == (1 if d == "up" else -1), de.at[o[g], "padj"] < 0.05)
                 for g, d in cds.direction.items() if g in o and o[g] in de.index and pd.notna(de.at[o[g], "log2FoldChange"])]
            v = np.array(v)
            conc[st] = (v[:, 0].mean(), (v[:, 0] & v[:, 1]).mean(), len(v))
        else:
            r = lo[(lo.model == "CDS score") & (lo.held_out == st)].iloc[0]
            conc[st] = (r.direction_agreement, sig[st], int(r.n_genes_tested))
    fig = plt.figure(figsize=(7.2, 6.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], width_ratios=[1.45, 1], hspace=0.55, wspace=0.32)
    ax = fig.add_subplot(gs[0, :])
    x = np.arange(len(order))
    col = [CAT[0] if st.startswith("WM") else CAT[1] for st in order]
    tot = [conc[st][0] * 100 for st in order]
    sg = [conc[st][1] * 100 for st in order]
    ax.bar(x, tot, 0.6, color=col, alpha=0.35, label="Same direction (any)")
    ax.bar(x, sg, 0.6, color=col, label="Same direction and padj < 0.05")
    ax.axhline(50, color=INK2, lw=0.8, ls="--")
    pass  # dashed line = chance (50%), explained in the legend
    for i, st in enumerate(order):
        me = vu.loc[st, "M_eff"]
        pe = stats.binomtest(int(round(conc[st][0] * me)), int(round(me)), 0.5, alternative="greater").pvalue
        ax.text(i, tot[i] + 1.5, f"{tot[i]:.0f}%\nn={conc[st][2]}\nP(M$_{{eff}}$)={pe:.0e}".replace("e-0", "e-"), ha="center",
                va="bottom", fontsize=5.6, color=INK)
    ax.set_xticks(x, [short[st] for st in order], fontsize=6.5)
    ax.set_ylim(0, 140)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Signature genes (%)")
    ax.axvline(3.5, color=AXIS, lw=0.8)
    ax.text(1.5, 131, "Held-out watermelon drought", ha="center", fontsize=7, color=INK2)
    ax.text(5, 131, "Cucumber drought (orthologues)", ha="center", fontsize=7, color=INK2)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=INK2, alpha=0.35, label="same direction"), Patch(color=INK2, label="same direction, padj < 0.05")],
              fontsize=6, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, frameon=False, borderaxespad=0.2)
    ax.set_title("Directional reproducibility", loc="left", pad=6)
    panel(ax, "a")
    # b: AUC vs random-gene-set null and permutation P
    ax = fig.add_subplot(gs[1, 0])
    for i, st in enumerate(order):
        r = vu.loc[st]
        ax.plot([i, i], [r.null_AUC_q025, r.null_AUC_q975], color=GRID, lw=6, solid_capstyle="butt", zorder=1)
        ax.scatter(i, r.null_AUC_mean, marker="_", s=80, color=GREY, zorder=2)
        ax.scatter(i, r.AUC, s=28, color=col[i], zorder=3, edgecolor="white", linewidth=0.8)
        ax.text(i, 1.06, f"{r.AUC_perm_p:.2g}", ha="center", fontsize=5.6, color=INK2)
    ax.text(-0.9, 1.06, "P$_{perm}$", ha="left", fontsize=5.6, color=INK2)
    ax.set_xticks(x, [short[st].replace("\n", " ") for st in order], fontsize=5.8, rotation=40, ha="right")
    ax.set_ylim(-0.03, 1.13)
    ax.set_ylabel("AUC (signature score)")
    ax.set_title("Sample classification vs random gene sets", loc="left", pad=10)
    panel(ax, "b")
    # c: cucumber sample scores
    ax = fig.add_subplot(gs[1, 1])
    st3 = ["CU_Leaf_PEG_2025", "CU_Leaf_PEG48h_2024", "CU_Leaf_graft_2022"]
    rng = np.random.default_rng(0)
    for i, stn in enumerate(st3):
        d = sc[sc.study == stn]
        for cond, cc, off in [("control", DOWN, -0.15), ("stress", UP, 0.15)]:
            v = d[d.condition == cond].score
            ax.scatter(np.full(len(v), i + off) + rng.uniform(-0.04, 0.04, len(v)), v, s=20, color=cc,
                       label=("Control" if cond == "control" else "Drought") if i == 0 else None, zorder=3)
    ax.axhline(0, color=AXIS, lw=0.6, ls="--")
    ax.set_xticks(range(3), [short[s_].replace("CU ", "") for s_ in st3], fontsize=6)
    ax.set_ylabel("Signature score (orthologues)")
    ax.legend(fontsize=6, loc="lower right")
    ax.set_title("Cucumber samples", loc="left", pad=10)
    panel(ax, "c")
    save(fig, "Fig3_validation")


def fig4():
    spec = pd.read_csv(R / "validation" / "specificity.tsv", sep="\t", index_col=0)
    tc = pd.read_csv(R / "validation" / "time_course.tsv", sep="\t", index_col=0)
    tol = pd.read_csv(R / "validation" / "tolerance.tsv", sep="\t", index_col=0)
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.8), gridspec_kw=dict(wspace=0.75, width_ratios=[1.15, 1, 1]))
    ax = axs[0]
    lo = pd.read_csv(R / "validation" / "loso.tsv", sep="\t")
    lo = lo[lo.model == "CDS score"]
    cs = pd.read_csv(R / "validation" / "cross_species.tsv", sep="\t")
    sb = pd.read_csv(R / "validation" / "specificity_by_contrast.tsv", sep="\t")
    rbh = pd.read_csv(M / "ref" / "rbh_wm_cu.tsv", sep="\t")
    o = dict(zip(rbh.watermelon, rbh.cucumber))
    cdsd = pd.read_csv(R / "cds.tsv", sep="\t", index_col=0).direction
    def cu_sig(st):
        de = pd.read_csv(R / "de" / f"{st}.tsv", sep="\t", index_col=0)
        v = [(np.sign(de.at[o[g], "log2FoldChange"]) == (1 if d == "up" else -1)) and de.at[o[g], "padj"] < 0.05
             for g, d in cdsd.items() if g in o and o[g] in de.index and pd.notna(de.at[o[g], "log2FoldChange"])]
        return np.mean(v)
    rows = [(LABEL[r.held_out].split(",")[0] + ", " + LABEL[r.held_out].split(", ")[1], r.direction_agreement, r.frac_sig_same_dir, CAT[0]) for r in lo.itertuples()]
    rows += [(LABEL[r.study], r.cds_direction_agreement, cu_sig(r.study), CAT[1]) for r in cs.itertuples()]
    rows += [(LABEL[r.contrast], r.direction_agreement, r.frac_sig_same_dir, CAT[2] if "salt" in r.contrast else CAT[3]) for r in sb.itertuples()]
    rows = rows[::-1]
    ax.barh(range(len(rows)), [r[1] * 100 for r in rows], color=[r[3] for r in rows], height=0.7, alpha=0.35)
    ax.barh(range(len(rows)), [r[2] * 100 for r in rows], color=[r[3] for r in rows], height=0.7)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=5.8)
    ax.axvline(50, color=INK2, lw=0.8, ls="--")
    ax.set_xlim(0, 100)
    ax.set_xlabel("CDS genes in the signature direction (%)\nlight: any change; dark: padj < 0.05")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=CAT[0], label="WM drought (held out)"), Patch(color=CAT[1], label="Cucumber drought"),
                       Patch(color=CAT[2], label="WM salt"), Patch(color=CAT[3], label="WM cold")],
              fontsize=5.5, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, borderaxespad=0.2)
    ax.set_title("Drought vs salt and cold", loc="left", pad=26)
    panel(ax, "a")
    ax = axs[1]
    tps = ["1h", "6h", "12h", "24h"]
    for d, col in [("up", UP), ("down", DOWN)]:
        s = tc[tc.direction == d][tps]
        med = s.median()
        q1, q3 = s.quantile(0.25), s.quantile(0.75)
        x = [1, 6, 12, 24]
        ax.plot(x, med, color=col, lw=2, marker="o", ms=4)
        ax.text(24.8, med.iloc[-1], f"CDS {d}", color=INK2, fontsize=6.5, va="center")
        ax.fill_between(x, q1, q3, color=col, alpha=0.15, lw=0)
    ax.axhline(0, color=AXIS, lw=0.6)
    ax.set_xticks([1, 6, 12, 24])
    ax.set_xlabel("Hours of PEG")
    ax.set_ylabel("log2FC vs 0 h (median, IQR)")
    ax.set_xlim(0, 30)
    ax.set_title("Activation dynamics", loc="left")
    panel(ax, "b")
    ax = axs[2]
    sgn = np.sign(tol.meta_log2FC)
    ax.scatter(tol.lfc_Y34 * sgn, tol.lfc_M20 * sgn, s=5, color=INK2, alpha=0.6)
    lim = np.nanmax(np.abs([tol.lfc_Y34 * sgn, tol.lfc_M20 * sgn])) * 1.05
    ax.plot([-2, lim], [-2, lim], color=AXIS, lw=0.8, ls="--")
    ax.set_xlim(-2, lim)
    ax.set_ylim(-2, lim)
    frac = ((tol.lfc_M20 * sgn) < (tol.lfc_Y34 * sgn))[tol.lfc_M20.notna() & tol.lfc_Y34.notna()].mean()
    ax.text(0.03, 0.97, f"{frac:.0%} below diagonal\n(weaker in M20)", transform=ax.transAxes, ha="left", va="top", fontsize=6, color=INK2)
    ax.set_xlabel("Response in susceptible Y34\n(signed log2FC)")
    ax.set_ylabel("Response in tolerant M20\n(signed log2FC)")
    ax.set_title("Genotype", loc="left")
    panel(ax, "c")
    save(fig, "Fig4_specificity_dynamics")


def fig5():
    go = pd.read_csv(R / "function" / "go_enrichment.tsv", sep="\t")
    ae = pd.read_csv(R / "function" / "arabidopsis_enrichment.tsv", sep="\t")
    pe = pd.read_csv(R / "function" / "promoter_enrichment.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.3, 1], hspace=0.55, wspace=0.75)
    ax = fig.add_subplot(gs[0, :])
    rows = []
    for s, col in [("CDS_up", UP), ("CDS_down", DOWN)]:
        g = go[(go.set == s) & (go.aspect == "biological_process") & (go.FDR < 0.05)]
        g = g[g.K <= 500].nsmallest(7, "p")
        for r in g.itertuples():
            rows.append((r.term[:60], -np.log10(r.FDR), col, r.k))
    rows = rows[::-1]
    ax.barh(range(len(rows)), [r[1] for r in rows], color=[r[2] for r in rows], height=0.7)
    ax.set_yticks(range(len(rows)), [f"{r[0]} ({r[3]})" for r in rows], fontsize=6.5)
    ax.axvline(-np.log10(0.05), color=AXIS, lw=0.6, ls="--")
    ax.set_xlabel("−log10 FDR (GO biological process; gene count in brackets)")
    ax.set_title("Enriched GO processes (red: induced, blue: repressed)", loc="left")
    panel(ax, "a")
    ax = fig.add_subplot(gs[1, 0])
    xs = ["CDS_up", "CDS_down", "CDS_all"]
    w = 0.38
    for j, (ev, lab, alpha) in enumerate([("at_drought_any", "Any evidence", 1.0),
                                          ("at_drought_functional", "Mutant/assay evidence", 0.45)]):
        a = ae[ae.evidence == ev].set_index("set")
        xpos = np.arange(3) + (j - 0.5) * w
        ax.bar(xpos, [a.loc[x, "set_rate"] * 100 for x in xs], w * 0.92, color=[UP, DOWN, INK2], alpha=alpha,
               edgecolor=[UP, DOWN, INK2], label=lab)
        ax.hlines(a.bg_rate.iloc[0] * 100, xpos[0] - w / 2, xpos[-1] + w / 2, color=GREY, lw=0.8, ls="--")
        for i, x in enumerate(xs):
            p = a.loc[x, "p"]
            ax.text(xpos[i], a.loc[x, "set_rate"] * 100 + 0.15, f"P={p:.1g}",
                    ha="center", va="bottom", fontsize=5.2, color=INK2, rotation=90 if j else 0)
    ax.set_xticks(range(3), ["Up", "Down", "All"])
    ax.set_ylim(0, 9.5)
    ax.set_ylabel("% genes with ortholog of an\nArabidopsis drought gene")
    ax.legend(fontsize=5.8, loc="upper right")
    ax.set_title("Orthologue annotation (indirect)", loc="left")
    panel(ax, "b")
    ax = fig.add_subplot(gs[1, 1])
    motifs = list(dict.fromkeys(pe.motif))
    for j, (s, col) in enumerate([("CDS_up", UP), ("CDS_down", DOWN)]):
        v = pe[pe.set == s].set_index("motif").reindex(motifs)
        ax.scatter(np.log2(v.odds_ratio.clip(lower=0.1)), np.arange(len(motifs)) + (0.12 if j == 0 else -0.12),
                   color=col, s=np.where(v.FDR < 0.05, 36, 12), label=s.replace("CDS_", "CDS "), zorder=3)
    ax.axvline(0, color=AXIS, lw=0.8)
    ax.set_yticks(range(len(motifs)), [m.replace("ABRE", "ABRE (≥2)") for m in motifs], fontsize=6.5)
    ax.set_xlabel("log2 odds ratio vs all genes\n(large dot: FDR < 0.05)")
    ax.legend(fontsize=6.5, loc="center right")
    ax.set_title("Promoter motifs (indirect)", loc="left")
    panel(ax, "c")
    save(fig, "Fig5_function_regulation")


def fig1():
    """Study design schematic."""
    from matplotlib.patches import FancyBboxPatch
    d = pd.read_csv(M / "census" / "design.tsv", sep="\t")
    n = d.groupby("study").size()
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 60)
    ax.axis("off")

    def box(x, y, w, h, title, lines, col):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", fc="white", ec=col, lw=1.4))
        ax.text(x + 1.2, y + h - 1.6, title, fontsize=7.5, fontweight="bold", color=INK, va="top")
        for i, l in enumerate(lines):
            ax.text(x + 1.2, y + h - 8.6 - i * 3.0, l, fontsize=6.1, color=INK2 if i else INK, va="top")

    groups = [("Discovery\nwatermelon drought", ["WM_Leaf_M20Y34_2018", "WM_Leaf_CS_2020", "WM_Root_PEG_2016",
                                                 "WM_Leaf_PEGtime_2024"], CAT[0]),
              ("Validation\ncucumber drought", ["CU_Leaf_PEG_2025", "CU_Leaf_PEG48h_2024", "CU_Leaf_graft_2022"], CAT[1]),
              ("Specificity\nwatermelon salt & cold", ["WM_Leaf_salt7h_2020", "WM_salt_LR_2025", "WM_Leaf_cold_2017"], CAT[2])]
    lab = dict(LABEL)
    lab["WM_salt_LR_2025"] = "Salt, leaf/root, 2 genotypes"
    for k, (t, st, col) in enumerate(groups):
        x = 1 + k * 33.3
        tot = int(sum(n[s] for s in st))
        box(x, 33, 31, 25, t, [f"{tot} libraries"] + [f"• {lab[s]} (n = {n[s]})" for s in st], col)
    steps = [("Uniform\nprocessing", "ENA reads, Salmon\n4–8 M reads/library"),
             ("Per-study\nDE", "PyDESeq2\nwith covariates"),
             ("Random-effects\nmeta-analysis", "DerSimonian–Laird\nCore Drought\nSignature"),
             ("Validation", "leave-one-study-out\ncucumber\northologues"),
             ("Biology", "specificity, timing\ngenotype, GO\npromoters, network")]
    for i, (t, sub) in enumerate(steps):
        x = 1 + i * 20
        ax.add_patch(FancyBboxPatch((x, 4), 17.5, 18, boxstyle="round,pad=0.3,rounding_size=1.2", fc="#f0efec", ec=AXIS, lw=0.8))
        ax.text(x + 8.75, 19, t, fontsize=6.6, fontweight="bold", ha="center", va="top", color=INK, wrap=True)
        ax.text(x + 8.75, 12.2, sub, fontsize=5.9, ha="center", va="top", color=INK2)
        if i < 4:
            ax.annotate("", xy=(x + 19.6, 13), xytext=(x + 18.0, 13), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.8))
    for k in range(3):
        ax.plot([16.5 + k * 33.3] * 2, [32.4, 27.5], color=CAT[k], lw=1.0)
    ax.plot([9.75, 83.1], [27.5, 27.5], color=INK2, lw=0.8)
    ax.annotate("", xy=(9.75, 22.6), xytext=(9.75, 27.5), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.8))
    ax.text(50, 24.6, "91 public RNA-seq libraries, 10 studies (ENA/SRA/DDBJ)", ha="center", fontsize=6.3, color=INK2)
    save(fig, "Fig1_study_design")


def fig6():
    mt = pd.read_csv(R / "network" / "modules.tsv", sep="\t").sort_values("r_drought")
    IND, REP = mt.iloc[-1].module, mt.iloc[0].module
    hubs = pd.read_csv(R / "network" / "hub_genes.tsv", sep="\t")
    o = pd.read_csv(R / "function" / "arabidopsis_orthologs.tsv", sep="\t", index_col=0)
    ps = pd.read_csv(R / "extended" / "network_per_study.tsv", sep="\t")
    fig = plt.figure(figsize=(7.4, 6.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], hspace=0.6, wspace=0.95)
    ax = fig.add_subplot(gs[0, 0])
    cols = [UP if (r > 0 and p < 0.001) else DOWN if (r < 0 and p < 0.001) else GRID
            for r, p in zip(mt.r_drought, mt.p_drought)]
    ax.barh(range(len(mt)), mt.r_drought, color=cols, height=0.7)
    ax.set_yticks(range(len(mt)), [f"{m} ({s_}; {c} CDS)" for m, s_, c in zip(mt.module, mt["size"], mt.n_cds)], fontsize=6)
    ax.axvline(0, color=AXIS, lw=0.8)
    ax.set_xlabel("Eigengene correlation\nwith drought (r)")
    ax.set_title("Modules", loc="left")
    panel(ax, "a")
    for j, (mod, col, letter) in enumerate([(IND, UP, "b"), (REP, DOWN, "c")]):
        ax = fig.add_subplot(gs[0, j + 1])
        h = hubs[hubs.module == mod].nlargest(15, "kME").iloc[::-1]
        names = [o.at[g, "at_symbol"] if g in o.index and isinstance(o.at[g, "at_symbol"], str) else g.replace("Cla97", "")
                 for g in h.gene]
        names = ["MFT" if n == "E12A11" else n.replace("AtATG18a", "ATG18a") for n in names]
        ax.scatter(h.kME, range(len(h)), s=np.where(h.in_cds, 30, 12), color=col, zorder=3)
        ax.hlines(range(len(h)), 0.8, h.kME, color=col, lw=0.8)
        ax.set_yticks(range(len(h)), names, fontsize=6)
        ax.set_xlim(0.8, 1.0)
        ax.set_xlabel("Module membership (kME)\nlarge dot: CDS gene")
        r = mt.set_index("module").loc[mod, "r_drought"]
        ax.set_title(f"{mod}: most connected genes\n(r = {r:+.2f}; candidates)", loc="left", fontsize=8)
        panel(ax, letter)
    ax = fig.add_subplot(gs[1, :])
    studies = ["WM_Leaf_M20Y34_2018", "WM_Leaf_CS_2020", "WM_Root_PEG_2016", "WM_Leaf_PEGtime_2024"]
    w = 0.36
    for k, (mod, col) in enumerate([(IND, UP), (REP, DOWN)]):
        d = ps[ps.module == mod].set_index("study").reindex(studies)
        xs = np.arange(4) + (k - 0.5) * w
        ax.bar(xs, d.eigengene_shift, w * 0.9, color=col, label=f"{mod} ({'induced' if k == 0 else 'repressed'})")
        for x_, v, p_, z_ in zip(xs, d.eigengene_shift, d.p, d.preservation_z):
            ax.text(x_, v + (0.03 if v > 0 else -0.03), f"P={p_:.1g}\nz={z_:.0f}", ha="center", va="bottom" if v > 0 else "top",
                    fontsize=5.5, color=INK2)
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.set_xticks(range(4), [LABEL[s_] for s_ in studies], fontsize=6.5)
    ax.set_ylim(-1.2, 1.2)
    ax.set_ylabel("Eigengene shift\n(drought − control)")
    ax.legend(fontsize=6.5, loc="lower left", ncol=2)
    ax.set_title("Module behaviour within each study (P: Welch t test; z: co-expression preservation)", loc="left", fontsize=8)
    panel(ax, "d")
    save(fig, "Fig6_network")


def figS1():
    t = pd.read_csv(R / "network" / "soft_threshold.tsv", sep="\t")
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.7), gridspec_kw=dict(wspace=0.45))
    sel = int(open(R / "network" / "selected_power.txt").read()) if (R / "network" / "selected_power.txt").exists() else 16
    ax = axs[0]
    ax.plot(t.power, t.R2, marker="o", color=INK, ms=4)
    ax.axhline(0.8, color=UP, lw=0.8, ls="--")
    ax.text(t.power.min(), 0.81, "R² = 0.80 criterion", fontsize=6.5, color=UP, va="bottom")
    r = t.set_index("power").R2[sel]
    ax.scatter([sel], [r], s=90, facecolor="none", edgecolor=UP, lw=1.5, zorder=3)
    ax.annotate(f"selected: {sel}\n(R² = {r:.2f})", (sel, r), xytext=(sel - 7, r - 0.35), fontsize=6.5,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    ax.set_xlabel("Soft-threshold power (β)")
    ax.set_ylabel("Scale-free fit (signed R²)")
    ax.set_title("Scale-free topology", loc="left")
    panel(ax, "a")
    ax = axs[1]
    ax.plot(t.power, t.mean_k, marker="o", color=INK, ms=4)
    ax.axvline(sel, color=UP, lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_xlabel("Soft-threshold power (β)")
    ax.set_ylabel("Mean connectivity")
    ax.set_title("Connectivity", loc="left")
    panel(ax, "b")
    save(fig, "FigS1_soft_threshold")


def figS2():
    """Sensitivity of the CDS."""
    d = pd.read_csv(R / "extended" / "meta_sensitivity.tsv", sep="\t")
    d = d[d.n_cds > 0].copy()
    lab = d.variant.str.replace("Thresholds: ", "").str.replace("Estimator: ", "").str.replace("Combination: ", "")
    lab = lab.str.replace("Leave out: ", "without ").str.replace("Original (DL, FDR<0.01, |lfc|>=1, >=3/4 same)", "Original CDS (DL)")
    lab = lab.str.replace("_", " ")
    d["label"] = lab
    grp_col = {"original": INK, "estimator": CAT[0], "threshold": GREY, "leave-one-study-out": CAT[1], "subsampling": CAT[2]}
    d = d.iloc[::-1]
    fig, axs = plt.subplots(1, 4, figsize=(7.4, 7.6), sharey=True, gridspec_kw=dict(wspace=0.12))
    y = np.arange(len(d))
    cols = [grp_col[g] for g in d.group]
    specs = [("n_cds", "Signature size (genes)", None), ("retained_from_original", "Original CDS genes\nretained (%)", 100),
             ("ABRE_OR_up", "Paired-ABRE odds ratio\n(induced genes)", None),
             ("cucumber_PEG_agreement", "Concordance (%)\ncucumber PEG / salt / cold", 100)]
    for ax, (c, t, sc) in zip(axs, specs):
        v = d[c] * (sc or 1)
        if c == "cucumber_PEG_agreement":
            ax.scatter(d.salt_leaf_agreement * 100, y, s=10, color=CAT[2], label="salt leaf", zorder=3)
            ax.scatter(d.cold_agreement * 100, y, s=10, color=CAT[3], label="cold", zorder=3)
            ax.scatter(v, y, s=10, color=CAT[1], label="cucumber PEG", zorder=3)
            ax.axvline(50, color=INK2, lw=0.6, ls="--")
            ax.set_xlim(40, 100)
            ax.legend(fontsize=5.5, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=1, frameon=False)
        else:
            ax.barh(y, v, color=cols, height=0.7)
            if c == "ABRE_OR_up":
                ax.axvline(1, color=INK2, lw=0.6, ls="--")
        ax.set_xlabel(t, fontsize=6.5)
        ax.tick_params(axis="x", labelsize=6)
        ax.grid(axis="x", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
    axs[0].set_yticks(y, d.label, fontsize=5.6)
    panel(axs[0], "a")
    save(fig, "FigS2_meta_sensitivity")


def figS3():
    """Read composition of cucumber libraries and subsampling robustness."""
    lm = pd.read_csv(R / "extended" / "lowmap_classification.tsv", sep="\t")
    lab = {"SRR30555221": "PRJNA1150737\nDT", "SRR30555222": "PRJNA1150737\nDS", "SRR30555223": "PRJNA1150737\nCT",
           "SRR30555224": "PRJNA1150737\nCS", "DRR784147": "PRJDB37115\nCTRL2", "DRR784148": "PRJDB37115\nCTRL3",
           "SRR14206075": "PRJNA720243\nP-WW"}
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 2, hspace=0.75, wspace=0.35)
    ax = fig.add_subplot(gs[0, :])
    cats = [("mRNA", "Nuclear mRNA (cDNA)", CAT[0]), ("genome_other", "Genome, other", "#9ec5f0"),
            ("rRNA", "rRNA", CAT[1]), ("chloroplast", "Chloroplast", CAT[2]), ("unmapped", "Unmapped", GRID)]
    left = np.zeros(len(lm))
    yy = np.arange(len(lm))[::-1]
    for c, n, colr in cats:
        ax.barh(yy, lm[c], left=left, color=colr, height=0.65, label=n)
        left += lm[c].values
    ax.set_yticks(yy, [lab[r].replace("\n", " ") for r in lm.run], fontsize=6.5)
    ax.set_xlabel("Reads (%)")
    ax.set_xlim(0, 100)
    ax.legend(fontsize=6, ncol=5, loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False)
    ax.set_title("Composition of cucumber read-1 sequences (200,000 per library)", loc="left", pad=16)
    panel(ax, "a")
    for k, (f0, f1, ttl, name) in enumerate([
            (R / "de" / "WM_Leaf_M20Y34_2018.tsv", R / "robustness" / "de" / "WM_Leaf_M20Y34_2018__first4M.tsv",
             "PRJNA454040: random 8 M vs first 4 M", ("random 8 M", "first 4 M")),
            (R / "de" / "WM_Root_PEG_2016.tsv", R / "extended" / "de" / "WM_Root_PEG_2016__random4M.tsv",
             "PRJNA326331: first 4 M vs random 4 M", ("first 4 M", "random 4 M"))]):
        ax = fig.add_subplot(gs[1, k])
        if not f1.exists():
            ax.text(0.5, 0.5, "pending", ha="center", transform=ax.transAxes)
            continue
        a = pd.read_csv(f0, sep="\t", index_col=0).log2FoldChange
        b = pd.read_csv(f1, sep="\t", index_col=0).log2FoldChange
        i = a.index.intersection(b.index)
        a, b = a[i], b[i]
        ok = a.notna() & b.notna()
        cds = pd.read_csv(R / "cds.tsv", sep="\t", index_col=0)
        isc = i[ok].isin(cds.index)
        ax.scatter(a[ok][~isc], b[ok][~isc], s=1.5, color=GRID, rasterized=True)
        ax.scatter(a[ok][isc], b[ok][isc], s=3, color=INK2, rasterized=True, label="CDS genes")
        lim = 10
        ax.plot([-lim, lim], [-lim, lim], color=UP, lw=0.6, ls="--")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        r = np.corrcoef(a[ok], b[ok])[0, 1]
        ax.text(0.03, 0.97, f"r = {r:.3f}\nn = {ok.sum():,} genes", transform=ax.transAxes, va="top", fontsize=6.5)
        ax.set_xlabel(f"log2FC ({name[0]})")
        ax.set_ylabel(f"log2FC ({name[1]})")
        ax.set_title(ttl, loc="left", fontsize=7.5)
        ax.legend(fontsize=6, loc="lower right", markerscale=3)
        panel(ax, "bc"[k])
    save(fig, "FigS3_reads_subsampling")


if __name__ == "__main__":
    import sys
    for f in (sys.argv[1:] or ["fig1", "fig2", "fig3", "fig4", "fig5", "fig6", "figS1", "figS2", "figS3"]):
        globals()[f]()
        print("made", f)
