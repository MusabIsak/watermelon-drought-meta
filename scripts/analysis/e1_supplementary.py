"""Supplementary Tables S1-S14 in one workbook."""
import os
from pathlib import Path

import pandas as pd

M = Path(os.environ.get("PROJECT_DIR", "."))
R = M / "results"
V = R / "extended"
out = R / "Supplementary_Tables_S1-S14.xlsx"

d = pd.read_csv(M / "census" / "design.tsv", sep="\t").merge(pd.read_csv(R / "qc.tsv", sep="\t"), on=["run", "crop"], how="left")
notes = {"SRR30555221": "PRJNA1150737 metadata conflict: treatment attribute/alias = drought (DT); experiment title = 'Drought tolerant treated'",
         "SRR30555222": "PRJNA1150737 metadata conflict: treatment attribute/alias = drought (DS); experiment title = 'Drought tolerant control'",
         "SRR30555223": "PRJNA1150737 metadata conflict: treatment attribute/alias = control (CT); experiment title = 'Drought Susceptible treated'",
         "SRR30555224": "PRJNA1150737 metadata conflict: treatment attribute/alias = control (CS); experiment title = 'Drought Susceptible control'",
         "SRR30555221_": ""}
d["metadata_note"] = d.run.map(notes).fillna("")
d.loc[d.run.isin(["SRR30555221", "SRR30555224"]), "metadata_note"] += "; deposited as three files, mate files used"
scr = pd.read_csv(M / "census" / "all_projects.tsv", sep="\t")
cds = pd.read_csv(R / "cds_annotated.tsv", sep="\t", index_col=0)
sp = pd.read_csv(R / "validation" / "specificity.tsv", sep="\t", index_col=0)
tc = pd.read_csv(R / "validation" / "time_course.tsv", sep="\t", index_col=0)
tol = pd.read_csv(R / "validation" / "tolerance.tsv", sep="\t", index_col=0)
s3 = sp.join(tc[["1h", "6h", "12h", "24h", "first_response"]]).join(tol[["lfc_M20", "lfc_Y34"]])
s3["class"] = s3["class"].replace({"drought-specific": "not detected under salt/cold", "shared abiotic": "shared with salt/cold"})
s1b = pd.read_csv(V / "S1b_root_summary.tsv", sep="\t", header=None, names=["metric", "value"]) if (V / "S1b_root_summary.tsv").exists() else pd.DataFrame()
s1a = pd.read_csv(R / "robustness" / "S1_summary.tsv", sep="\t", header=None, names=["metric", "value"])


def s7a():
    """validation_uncertainty.tsv plus the DESeq2-based concordance and its M_eff-adjusted binomial P value
    (the values reported in the text and Fig. 3a)."""
    from scipy import stats as _st
    vu = pd.read_csv(V / "validation_uncertainty.tsv", sep="\t")
    lo = pd.read_csv(R / "validation" / "loso.tsv", sep="\t")
    lo = lo[lo.model == "CDS score"].set_index("held_out")
    cs = pd.read_csv(R / "validation" / "cross_species.tsv", sep="\t").set_index("study")
    c, p = [], []
    for r in vu.itertuples():
        conc = lo.at[r.study, "direction_agreement"] if r.study in lo.index else cs.at[r.study, "cds_direction_agreement"]
        c.append(conc)
        p.append(_st.binomtest(int(round(conc * r.M_eff)), int(round(r.M_eff)), 0.5, alternative="greater").pvalue)
    vu.insert(vu.columns.get_loc("concordance"), "concordance_DESeq2_primary", c)
    vu["binom_p_Meff_DESeq2_primary"] = p
    vu = vu.rename(columns={"concordance": "concordance_simple_estimator", "binom_p_Meff": "binom_p_Meff_simple_estimator"})
    return vu

sheets = {
    "S1a_samples_QC": (d.drop(columns=["fastq_ftp", "fastq_md5", "fastq_bytes"], errors="ignore"), False,
                       "Table S1a. Sample sheet, group assignment, metadata notes and Salmon QC"),
    "S1b_projects_screened": (scr, False, "Table S1b. All cucurbit BioProjects retrieved for screening"),
    "S2_CDS": (cds, True, "Table S2. Core Drought Signature (845 genes): study-level and meta statistics, annotation"),
    "S3_salt_cold_time_genotype": (s3, True, "Table S3. Salt/cold classification, PEG time course, M20 vs Y34"),
    "S4a_LOSO": (pd.read_csv(R / "validation" / "loso.tsv", sep="\t"), False, "Table S4a. Leave-one-study-out validation"),
    "S4b_cross_species": (pd.read_csv(R / "validation" / "cross_species.tsv", sep="\t"), False, "Table S4b. Cucumber validation"),
    "S4c_salt_cold_by_contrast": (pd.read_csv(R / "validation" / "specificity_by_contrast.tsv", sep="\t"), False,
                                  "Table S4c. Agreement with each salt/cold contrast"),
    "S5a_GO": (pd.read_csv(R / "function" / "go_enrichment.tsv", sep="\t"), False, "Table S5a. GO enrichment"),
    "S5b_Arabidopsis": (pd.read_csv(R / "function" / "arabidopsis_enrichment.tsv", sep="\t"), False,
                        "Table S5b. Arabidopsis water-deprivation orthologue enrichment"),
    "S5c_promoters": (pd.read_csv(R / "function" / "promoter_enrichment.tsv", sep="\t"), False, "Table S5c. Promoter motifs"),
    "S5d_TFs": (pd.read_csv(R / "function" / "cds_transcription_factors.tsv", sep="\t"), False, "Table S5d. CDS transcription factors"),
    "S5e_modules": (pd.read_csv(R / "network" / "modules.tsv", sep="\t"), False, "Table S5e. Co-expression modules"),
    "S5f_module_membership": (pd.read_csv(R / "network" / "gene_modules.tsv", sep="\t"), False, "Table S5f. Module membership (kME)"),
    "S5g_contrasts": (pd.read_csv(R / "contrasts.tsv", sep="\t"), False, "Table S5g. All DE contrasts"),
    "S6_meta_sensitivity": (pd.read_csv(V / "meta_sensitivity.tsv", sep="\t"), False,
                            "Table S6. CDS sensitivity to estimator, thresholds, study exclusion and subsampling"),
    "S7a_validation_uncertainty": (s7a(), False,
                                  "Table S7a. AUC uncertainty and dependence-aware concordance tests"),
    "S7b_competitive_null": (pd.read_csv(V / "validation_competitive.tsv", sep="\t"), False,
                             "Table S7b. Label-permutation and competitive (random gene set) nulls for directional concordance"),
    "S8_read_composition": (pd.read_csv(V / "lowmap_classification.tsv", sep="\t"), False,
                            "Table S8. Read composition of cucumber libraries (%)"),
    "S8b_PRJNA1150737_genotype": (pd.read_csv(V / "prjna1150737_genotype_check.tsv", sep="\t"), False,
                                  "Table S8b. Genotype identity of PRJNA1150737 libraries (shared homozygous variants)"),
    "S8c_PRJNA1150737_sites": (pd.read_csv(V / "prjna1150737_genotype_sites.tsv", sep="\t"), False,
                               "Table S8c. Homozygous-differential sites used for S8b"),
    "S9a_detectability_summary": (pd.read_csv(V / "specificity_detectability_summary.tsv", sep="\t"), False,
                                  "Table S9a. Effect sizes/detectability summary"),
    "S9b_detectability_genes": (pd.read_csv(V / "specificity_detectability.tsv", sep="\t"), False,
                                "Table S9b. The 195 genes not detected under salt/cold, per contrast"),
    "S10a_network_per_study": (pd.read_csv(V / "network_per_study.tsv", sep="\t"), False, "Table S10a. Module behaviour per study"),
    "S10b_network_rebuilds": (pd.read_csv(V / "network_rebuilds.tsv", sep="\t"), False, "Table S10b. Network rebuilds"),
    "S10c_hub_stability": (pd.read_csv(V / "network_hub_stability.tsv", sep="\t"), False, "Table S10c. Stability of highly connected genes"),
    "S11a_subsample_PRJNA454040": (s1a, False, "Table S11a. First-4M vs random-8M, PRJNA454040"),
    "S11b_subsample_PRJNA326331": (s1b, False, "Table S11b. First-4M vs random-4M, PRJNA326331"),
    "S12_candidates": (pd.read_csv(V / "candidate_priority.tsv", sep="\t"), False,
                       "Table S12. Prioritised candidate genes for functional validation (computational evidence only)"),
    "S13_candidate_evidence": (pd.read_csv(V / "candidate_evidence.tsv", sep="\t"), False,
                               "Table S13. Evidence profile of all 845 signature genes in three levels: L1 transcriptional "
                               "reproducibility (this study); L2 predicted regulatory association (watermelon promoter "
                               "motifs, co-expression module) and public Arabidopsis experiments on the orthologue "
                               "(GSE110079 ABA-deficient mutants, GSE80565 ABA time course, GSE80564 ChIP-seq); L3 "
                               "published genetic evidence for the Arabidopsis orthologue or a cucurbit homologue. "
                               "No gene has been functionally tested in watermelon."),
    "S14a_public_sets": (pd.read_csv(V / "public_functional_sets.tsv", sep="\t"), False,
                         "Table S14a. Public Arabidopsis experimental data for orthologues of signature genes and modules "
                         "versus orthologues of all tested genes and of non-reproducible drought-responsive genes "
                         "(Fisher exact / Mann-Whitney tests; best-hit and RBH-only orthology)"),
    "S14b_public_genes": (pd.read_csv(V / "public_functional_genes.tsv", sep="\t"), False,
                          "Table S14b. Per-gene public Arabidopsis data for orthologues of signature and module genes"),
}
with pd.ExcelWriter(out) as w:
    pd.DataFrame({"Sheet": list(sheets), "Content": [v[2] for v in sheets.values()]}).to_excel(w, sheet_name="README", index=False)
    for k, (df, idx, _) in sheets.items():
        df = df.replace({"drought-specific": "not detected under salt/cold", "shared abiotic": "shared with salt/cold"})
        df.to_excel(w, sheet_name=k[:31], index=idx)
print("written", out)
