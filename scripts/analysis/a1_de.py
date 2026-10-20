"""Gene matrices + per-study / per-contrast DESeq2 for all 10 studies.

Outputs (results/):
  counts_<crop>.tsv.gz, tpm_<crop>.tsv.gz, samples.tsv, qc.tsv
  de/<contrast>.tsv  (baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
  contrasts.tsv      (one row per contrast: study, role, stress, crop, tissue, n_ctrl, n_stress, design)
"""
import os
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

M = Path(os.environ.get("PROJECT_DIR", ".") + "")
Q = M / "quant"
R = M / "results"
(R / "de").mkdir(parents=True, exist_ok=True)


def tx2gene_cucumber():
    t = pd.read_csv(M / "ref" / "cu_tx2gene.tsv", sep="\t", header=None, names=["tx", "gene"])
    return dict(zip(t.tx, t.gene))


def build_matrices(design):
    cu_map = tx2gene_cucumber()
    qc = []
    for crop, d in design.groupby("crop"):
        counts, tpm = {}, {}
        for run in d.run:
            f = Q / run / "quant.sf"
            if not f.exists():
                continue
            q = pd.read_csv(f, sep="\t", index_col=0)
            if crop == "cucumber":
                g = q.index.map(lambda x: cu_map.get(x, x))
            else:
                g = q.index.map(lambda x: re.sub(r"\.\d+$", "", x))
            counts[run] = q.NumReads.groupby(g).sum()
            tpm[run] = q.TPM.groupby(g).sum()
            info = json.loads((Q / run / "aux_info" / "meta_info.json").read_text())
            qc.append(dict(run=run, crop=crop, processed=info["num_processed"], mapped=info["num_mapped"],
                           pct_mapped=info["percent_mapped"]))
        pd.DataFrame(counts).fillna(0).round().astype(int).to_csv(R / f"counts_{crop}.tsv.gz", sep="\t")
        pd.DataFrame(tpm).fillna(0).to_csv(R / f"tpm_{crop}.tsv.gz", sep="\t")
    pd.DataFrame(qc).to_csv(R / "qc.tsv", sep="\t", index=False)


def run_deseq(counts, meta, design_terms, name, info):
    inference = DefaultInference(n_cpus=1)
    c = counts[meta.index].T
    keep = (c >= 10).sum(axis=0) >= max(2, min(meta.condition.value_counts()))
    c = c.loc[:, keep]
    dds = DeseqDataSet(counts=c, metadata=meta, design="~" + " + ".join(design_terms),
                       refit_cooks=True, inference=inference, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "stress", "control"], inference=inference, quiet=True)
    st.summary()
    res = st.results_df
    res.to_csv(R / "de" / f"{name}.tsv", sep="\t")
    sig = res[(res.padj < 0.05) & (res.log2FoldChange.abs() >= 1)]
    info.update(contrast=name, n_ctrl=int((meta.condition == "control").sum()),
                n_stress=int((meta.condition == "stress").sum()), design="~" + " + ".join(design_terms),
                genes_tested=len(res), up=int((sig.log2FoldChange > 0).sum()), down=int((sig.log2FoldChange < 0).sum()))
    print(f"  {name}: {info['n_ctrl']}C/{info['n_stress']}S, {len(res)} genes, up {info['up']} down {info['down']}",
          flush=True)
    return info


def covariates(m):
    terms = []
    for cov in ["genotype", "timepoint", "tissue"]:
        v = m[cov].fillna("")
        if v.nunique() > 1 and (v != "").all():
            # only adjust for covariates that are crossed with condition
            if m.groupby(cov).condition.nunique().min() == 2:
                terms.append(cov)
    return terms + ["condition"]


def main(only=None):
    design = pd.read_csv(M / "census" / "design.tsv", sep="\t", dtype=str).fillna("")
    design = design[[ (Q / r / "quant.sf").exists() for r in design.run ]]
    build_matrices(design)
    design.set_index("run").to_csv(R / "samples.tsv", sep="\t")
    rows = []
    for (crop, study), m in design.groupby(["crop", "study"]):
        if only and study not in only:
            continue
        if m.condition.nunique() < 2:
            continue
        counts = pd.read_csv(R / f"counts_{crop}.tsv.gz", sep="\t", index_col=0)
        m = m.set_index("run")
        base = dict(study=study, crop=crop, role=m.role.iloc[0], stress=m.stress.iloc[0],
                    tissue=",".join(sorted(m.tissue.unique())))
        print(study, flush=True)
        if study == "WM_Leaf_PEGtime_2024":
            # pooled (all PEG time points vs 0 h) for the meta-analysis, plus each time point
            rows.append(run_deseq(counts, m, ["condition"], study, dict(base)))
            for tp in ["1h", "6h", "12h", "24h"]:
                mm = m[(m.condition == "control") | (m.timepoint == tp)]
                rows.append(run_deseq(counts, mm, ["condition"], f"{study}__{tp}", dict(base, role="time")))
        elif study == "WM_salt_LR_2025":
            for t in ["leaf", "root"]:
                mm = m[m.tissue == t]
                rows.append(run_deseq(counts, mm, covariates(mm), f"{study}__{t}", dict(base, tissue=t)))
        else:
            rows.append(run_deseq(counts, m, covariates(m), study, dict(base)))
            if study == "WM_Leaf_M20Y34_2018":
                for g in ["M20", "Y34"]:
                    mm = m[m.genotype == g]
                    rows.append(run_deseq(counts, mm, covariates(mm), f"{study}__{g}", dict(base, role="genotype")))
    out = R / "contrasts.tsv"
    new = pd.DataFrame(rows)
    if out.exists() and only:
        old = pd.read_csv(out, sep="\t")
        new = pd.concat([old[~old.contrast.isin(new.contrast)], new])
    new.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main(sys.argv[1:] or None)
