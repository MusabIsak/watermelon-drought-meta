"""evidence profile of prioritised candidate genes, separating three levels of evidence:
  Level 1 - transcriptional reproducibility (this study: discovery studies, held-out cucumber, salt/cold);
  Level 2 - predicted regulatory association (watermelon promoter motifs, co-expression module) and, as public
            experimental support for that prediction, Arabidopsis orthologue data (ABA dependence in aba2-2/aba3-1,
            ABA responsiveness, ChIP-seq binding by ABA-responsive TFs; e3_public_functional.py);
  Level 3 - experimentally established function, from published genetic studies (Europe PMC screen of the 40
            highest-priority candidates plus the ABA-signalling and hub genes discussed in the text) of the Arabidopsis orthologue or
            a cucurbit homologue (curated below; every entry was checked against the Europe PMC record).
No candidate has been functionally tested in watermelon.
Output: results/extended/candidate_evidence.tsv (Table S13; main-text Table 3 is a subset)
"""
import os
from pathlib import Path

import pandas as pd

R = Path(os.environ.get("PROJECT_DIR", ".") + "/results")
O = R / "extended"

# AGI -> (published genetic evidence, reference(s), direction of effect on drought/ABA response)
LIT = {
    "AT3G14440": ("Knockout reduced ABA accumulation and drought tolerance; overexpression reduced transpiration and "
                  "improved drought tolerance", "Iuchi et al., 2001", "positive"),
    "AT4G26080": ("Loss-of-function alleles hypersensitive to ABA; negative regulator of ABA signalling",
                  "Gosti et al., 1999", "negative (feedback)"),
    "AT1G72770": ("Overexpression ABA-insensitive, loss of function ABA-hypersensitive; abi1 hab1 double mutant more "
                  "ABA-sensitive with reduced water consumption", "Saez et al., 2004, 2006", "negative (feedback)"),
    "AT5G51760": ("ahg1 mutant hypersensitive to ABA and osmotica during germination; no adult phenotype",
                  "Nishimura et al., 2007", "negative (seed)"),
    "AT3G11410": ("Overexpression ABA-insensitive (germination, stomatal closure); knockout ABA-hypersensitive",
                  "Kuhn et al., 2006", "negative (feedback)"),
    "AT4G34000": ("Overexpression ABA-hypersensitive with reduced transpiration and improved drought tolerance; "
                  "areb1 areb2 abf3 triple mutant less drought tolerant", "Kang et al., 2002; Yoshida et al., 2010",
                  "positive"),
    "AT2G46270": ("Overexpression improved osmotic, salt and drought tolerance; loss of function increased drought "
                  "sensitivity", "Ramegowda et al., 2017", "positive"),
    "AT2G46680": ("Induced by water deficit and ABA; with ATHB12 activates PP2C and represses PYL5/PYL8 genes "
                  "(ChIP, mutants): negative feedback on ABA signalling", "Söderman et al., 1996; Valdés et al., 2012",
                  "negative (feedback)"),
    "AT3G15500": ("Overexpression (with ANAC019/072) improved drought tolerance", "Tran et al., 2004", "positive"),
    "AT1G01720": ("Arabidopsis: knockout improved drought recovery, overexpression improved drought tolerance "
                  "(conflicting); cucumber CsATAF1 (Csa4M361820, the RBH orthologue): RNAi reduced and "
                  "overexpression increased drought tolerance, ABA-dependent",
                  "Lu et al., 2007; Wu et al., 2009; Wang et al., 2018", "context-dependent (positive in cucumber)"),
    "AT3G63060": ("Overexpression and knock-down: positive regulator of ABA responses; regulates anthocyanin "
                  "accumulation under drought", "Koops et al., 2011", "positive (ABA signalling)"),
    "AT4G37790": ("abig1 (HAT22) mutants: reduced ABA growth inhibition and fewer senesced leaves under drought; "
                  "not required for stomatal closure", "Liu et al., 2016", "growth restraint / senescence"),
    "AT1G54100": ("Overexpression conferred osmotic and oxidative stress tolerance; knockout more sensitive to "
                  "dehydration and salt", "Kotchoni et al., 2006", "positive"),
    "AT1G02205": ("Arabidopsis: controls very-long-chain alkane synthesis; overexpression reduced cuticle permeability "
                  "and soil water deficit susceptibility; cucumber CsCER1: altered expression changed alkanes, "
                  "cuticle permeability and drought resistance", "Bourdenx et al., 2011; Wang et al., 2015",
                  "positive"),
    "AT4G34220": ("RDK1: ABA-induced receptor-like kinase; rdk1 mutants ABA-resistant and drought-hypersensitive, "
                  "overexpression ABA-hypersensitive (note: the watermelon orthologue is repressed by drought)",
                  "Kumar et al., 2017", "positive (ABA signalling)"),
    "AT2G40340": ("Overexpression: dehydration-sensitive but heat- and freezing-tolerant; interacts with ABF2/3/4",
                  "Lee et al., 2010", "heat/freezing (not drought)"),
}


def main():
    c = pd.read_csv(O / "candidate_priority.tsv", sep="\t")
    pf = pd.read_csv(O / "public_functional_genes.tsv", sep="\t")
    pf = pf[pf.gene_set.isin(["CDS induced", "CDS repressed"])].drop_duplicates("gene").set_index("gene")
    gm = pd.read_csv(R / "network/gene_modules.tsv", sep="\t")
    gm["kME_rank"] = gm.groupby("module").kME.rank(ascending=False).astype(int)
    gm["module_size"] = gm.groupby("module").gene.transform("size")
    gm = gm.set_index("gene")
    rows = []
    for _, r in c.iterrows():
        g = r.Name
        p = pf.loc[g] if g in pf.index else None
        lit = LIT.get(r.AGI)
        rows.append({
            "watermelon_gene": g, "Arabidopsis_orthologue": r.AGI, "symbol": {"AT4G34220": "RDK1"}.get(r.AGI, r.At), "direction": r.direction,
            "priority_rank": _ + 1,
            # Level 1
            "L1_discovery_studies_same_direction": r.n_same_dir, "L1_meta_log2FC": round(r.meta_log2FC, 2),
            "L1_cucumber_PEG_log2FC": r.cucumber_PEG_log2FC, "L1_salt_cold_class": r.salt_cold_class,
            "L1_first_response_PEG": r.first_response_PEG,
            # Level 2 (watermelon predictions)
            "L2_paired_ABRE_promoter": (r.ABRE_n >= 2) if pd.notna(r.ABRE_n) else None,
            "L2_module": r.module,
            "L2_kME_rank": f"{gm.kME_rank[g]}/{gm.module_size[g]}" if g in gm.index and gm.module[g] in ("M1", "M3") else "",
            # Level 2 (public Arabidopsis experiments on the orthologue)
            "L2_At_osmotic_responsive": None if p is None else p.get("osmotic_responsive_WT"),
            "L2_At_ABA_dependent_aba2": None if p is None else p.get("aba2_dependent"),
            "L2_At_ABA_dependent_aba3": None if p is None else p.get("aba3_dependent"),
            "L2_At_ABA_responsive": None if p is None else p.get("ABA_responsive"),
            "L2_At_bound_ABF": None if p is None else p.get("ABF_bound"),
            "L2_At_n_TFs_bound_of_21": None if p is None else p.get("n_TFs_bound"),
            # Level 3
            "L3_published_function": lit[0] if lit else "",
            "L3_effect": lit[2] if lit else "",
            "L3_reference": lit[1] if lit else "",
            "tested_in_watermelon": "no",
        })
    E = pd.DataFrame(rows)
    E.to_csv(O / "candidate_evidence.tsv", sep="\t", index=False)
    lit_set = E[E.L3_reference != ""]
    print(len(E), "candidates;", len(lit_set), "with published genetic evidence")
    pd.set_option("display.width", 250)
    print(lit_set[["priority_rank", "watermelon_gene", "symbol", "L1_discovery_studies_same_direction",
                   "L1_cucumber_PEG_log2FC", "L2_paired_ABRE_promoter", "L2_kME_rank", "L2_At_ABA_dependent_aba2",
                   "L2_At_bound_ABF", "L3_reference"]].to_string())


if __name__ == "__main__":
    main()
