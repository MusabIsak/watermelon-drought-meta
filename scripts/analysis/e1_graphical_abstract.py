"""Graphical abstract: wording avoids presenting predicted regulation as mechanism."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

import a6_figures as F

fig, ax = plt.subplots(figsize=(7.0, 3.1))
ax.set_xlim(0, 110)
ax.set_ylim(-2, 44)
ax.axis("off")


def box(x, y, w, h, title, lines, col, fc="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", fc=fc, ec=col, lw=1.6))
    ax.text(x + w / 2, y + h - 2.2, title, ha="center", va="top", fontsize=8.4, fontweight="bold", color=F.INK)
    for i, l in enumerate(lines):
        ax.text(x + w / 2, y + h - 8.4 - i * 4.5, l, ha="center", va="top", fontsize=6.8, color=F.INK2)


box(1, 7, 25, 31, "Public data only", ["91 RNA-seq libraries", "10 studies, 2 species", "watermelon drought (4)",
                                       "cucumber drought (3)", "salt and cold (3)"], F.CAT[0])
box(33, 7, 33, 31, "Core drought signature", ["random-effects meta-analysis", "845 genes (299 up, 546 down)",
                                              "conclusions robust to model,", "thresholds and subsampling",
                                              "ABA-associated genes; ABRE-rich"], F.INK, fc="#f0efec")
for x0, x1 in [(26.8, 32.2), (66.8, 72.2)]:
    ax.annotate("", xy=(x1, 22.5), xytext=(x0, 22.5), arrowprops=dict(arrowstyle="-|>", color=F.INK2, lw=1.2))
vals = [("Held-out watermelon drought", 0.65, 0.93, "65–93%", F.CAT[0]), ("Cucumber drought (2 studies)", 0.70, 0.91, "70–91%", F.CAT[1]),
        ("Watermelon salt", 0.64, 0.89, "64–89%", F.CAT[2]), ("Watermelon cold", 0.51, 0.51, "51%", F.CAT[3])]
x0, L = 74, 26
ax.text(x0 + L / 2 + 2, 38.5, "Genes changing as predicted", ha="center", va="top", fontsize=8.4, fontweight="bold", color=F.INK)
ax.plot([x0 + L * 0.5] * 2, [6.5, 30.5], color=F.GREY, lw=0.7, ls="--", zorder=0)
for i, (lab, lo, hi, txt, c) in enumerate(vals):
    y = 29 - i * 6.4
    ax.text(x0, y + 1.9, lab, fontsize=6.4, color=F.INK2, va="bottom")
    ax.add_patch(FancyBboxPatch((x0, y - 1.3), L * lo, 2.8, boxstyle="round,pad=0,rounding_size=0.5", fc=c, ec="none"))
    if hi > lo:
        ax.add_patch(FancyBboxPatch((x0 + L * lo, y - 1.3), L * (hi - lo), 2.8, boxstyle="round,pad=0,rounding_size=0.5",
                                    fc=c, ec="none", alpha=0.4))
    ax.text(x0 + L * hi + 0.8, y + 0.1, txt, fontsize=6.6, va="center", color=F.INK)
ax.text(x0 + L * 0.5, 4.6, "chance (50%)", fontsize=5.9, ha="center", color=F.INK2)
ax.text(55, 0.5, "Computational associations and predictions; candidate regulators require experimental validation",
        ha="center", fontsize=6.2, color=F.INK2, style="italic")
for ext, kw in [("tiff", dict(dpi=330, pil_kwargs={"compression": "tiff_lzw"})), ("pdf", {}), ("png", dict(dpi=200))]:
    fig.savefig(F.FIG / f"Graphical_abstract.{ext}", bbox_inches="tight", **kw)
print("graphical abstract written")
