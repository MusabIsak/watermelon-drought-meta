"""effect sizes and detectability for the 195 CDS genes not detected as responsive under
the available salt and cold conditions. For every gene x salt/cold contrast: expressed/tested or not, log2FC with
95% CI, and classification:
  not_tested      - filtered out (too low expression in that study)
  same_dir_trend  - same direction as drought, not significant
  opposite_trend  - opposite direction, not significant
  small_effect    - 95% CI entirely within (-1, 1): a response of drought-like size (|log2FC|>=1) is excluded
(categories overlap: small_effect is reported separately). Outputs: results/extended/specificity_detectability.tsv
and specificity_detectability_summary.tsv
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(os.environ.get("PROJECT_DIR", ".") + "/results")
O = R / "extended"
sp = pd.read_csv(R / "validation" / "specificity.tsv", sep="\t", index_col=0)
con = pd.read_csv(R / "contrasts.tsv", sep="\t")
spec = con[con.role == "specificity"].contrast.tolist()
rows, summ = [], []
for cls in ["drought-specific", "shared abiotic"]:
    genes = sp.index[sp["class"] == cls]
    for c in spec:
        d = pd.read_csv(R / "de" / f"{c}.tsv", sep="\t", index_col=0).reindex(genes)
        sgn = np.sign(sp.loc[genes, "meta_log2FC"])
        tested = d.log2FoldChange.notna()
        lo, hi = d.log2FoldChange - 1.96 * d.lfcSE, d.log2FoldChange + 1.96 * d.lfcSE
        small = tested & (lo > -1) & (hi < 1)
        same = tested & (np.sign(d.log2FoldChange) == sgn)
        summ.append(dict(class_=cls, contrast=c, n=len(genes), tested=int(tested.sum()),
                         same_direction=int(same.sum()), same_direction_pct=same.sum() / tested.sum() * 100,
                         small_effect_CI_within_1=int(small.sum()), small_effect_pct=small.sum() / tested.sum() * 100,
                         median_abs_lfc=float(d.log2FoldChange.abs().median()),
                         median_drought_abs_lfc=float(sp.loc[genes, "meta_log2FC"].abs().median()),
                         median_lfcSE=float(d.lfcSE.median())))
        if cls == "drought-specific":
            for g in genes:
                rows.append(dict(gene=g, contrast=c, drought_meta_log2FC=sp.at[g, "meta_log2FC"],
                                 log2FC=d.at[g, "log2FoldChange"], lfcSE=d.at[g, "lfcSE"], padj=d.at[g, "padj"],
                                 baseMean=d.at[g, "baseMean"], small_effect=bool(small[g])))
S = pd.DataFrame(summ)
S.to_csv(O / "specificity_detectability_summary.tsv", sep="\t", index=False)
G = pd.DataFrame(rows)
G.to_csv(O / "specificity_detectability.tsv", sep="\t", index=False)
print(S.round(2).to_string(index=False))
# per gene: how many contrasts exclude a drought-sized response
k = G.groupby("gene").small_effect.sum()
print("genes with drought-sized response excluded in >=1 / >=2 / all 4 contrasts:",
      (k >= 1).sum(), (k >= 2).sum(), (k == 4).sum(), "of", len(k))
t = G.groupby("gene").log2FC.apply(lambda x: x.notna().sum())
print("genes tested in all 4 contrasts:", (t == 4).sum())
