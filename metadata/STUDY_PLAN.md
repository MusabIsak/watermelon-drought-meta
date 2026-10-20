# Study plan (fixed before analysis of validation data)

Working title: **A conserved, drought-specific transcriptional signature in watermelon revealed by
multi-study meta-analysis and validated across tissues, time and species**

Date fixed: 12 October 2026. Validation datasets (cucumber) are processed but not inspected until the
signature is frozen.

## Data (91 public RNA-seq libraries, NCBI SRA / ENA)

| Role | Study label | Accession | Samples (C/S) | Design |
|---|---|---|---|---|
| Discovery | WM_Leaf_M20Y34_2018 | PRJNA454040 | 8/8 | leaf; tolerant M20 vs sensitive Y34; 4 and 8 d water withholding |
| Discovery | WM_Leaf_CS_2020 | PRJNA604984 | 3/3 | leaf; Crimson Sweet; 8 d drought |
| Discovery | WM_Root_PEG_2016 | PRJNA326331 | 3/3 | root; M08; 6 h 20% PEG-6000 (duplicate upload SRR3706796 removed) |
| Discovery | WM_Leaf_PEGtime_2024 | PRJNA1055702 | 3/12 | leaf; PEG time course 0, 1, 6, 12, 24 h |
| Validation | CU_Leaf_PEG_2025 | PRJDB37115 | 3/3 | cucumber leaf; 20% PEG |
| Validation | CU_Leaf_PEG48h_2024 | PRJNA1150737 | 2/2 | cucumber leaf; tolerant & sensitive lines; PEG 48 h |
| Validation | CU_Leaf_graft_2022 | PRJNA720243 | 2/2 | cucumber leaf; self- vs pumpkin-grafted; water stress |
| Specificity | WM_Leaf_salt7h_2020 | PRJNA609260 | 3/3 | leaf; 7 h salt |
| Specificity | WM_salt_LR_2025 | PRJNA844416 | 12/12 | leaf + root; 2 genotypes; 150 mM NaCl |
| Specificity | WM_Leaf_cold_2017 | PRJNA328189 | 2/2 | leaf; 4 °C |

Excluded after inspection: studies without controls (PRJNA1429590), unreadable group labels
(PRJNA770012, PRJNA1050590, PRJNA1047768), single-sample or unreplicated designs without controls
(PRJNA246706, PRJNA219371, PRJNA952866), non-drought topics returned by keyword search.

## Processing
Reads streamed from ENA, randomly subsampled to 8 M reads (fixed seed, identical draw for both mates),
quantified with Salmon 1.10 (selective alignment, GC/sequence bias correction) against 97103 v2.5
transcripts (CuGenDB); gene-level aggregation.

## Primary analysis (watermelon drought, discovery)
1. Per-study DESeq2 (covariates: genotype, time point where present; time course: each time vs 0 h and pooled).
2. Random-effects meta-analysis (DerSimonian–Laird) of per-study log2FC; heterogeneity (I², Cochran Q);
   BH-FDR.
3. **Core Drought Signature (CDS)**: meta-FDR < 0.01, |meta log2FC| ≥ 1, same direction in ≥ 3 of 4 studies.
4. Leave-one-study-out: signature rebuilt without the held-out study; signature score
   (mean z of up genes − mean z of down genes) evaluated by AUC; null = 1,000 random gene sets of equal size.
5. Benchmark: ML classifiers (RF, XGBoost, L1-logistic, SVM) under the same LOSO scheme.

## Secondary analyses
6. **Cross-species validation** (cucumber, untouched until signature frozen): orthologs by reciprocal best
   DIAMOND hits; signature score AUC; correlation of ortholog log2FC.
7. **Specificity**: score in salt and cold datasets; split CDS into drought-specific vs shared abiotic genes.
8. **Temporal dynamics**: activation of CDS genes across 1–24 h PEG (early vs late responders).
9. **Tolerance**: genotype × drought interaction (M20 vs Y34) for CDS genes.
10. **Regulatory layer**: co-expression network (WGCNA) on discovery data; hub genes; transcription-factor
    families; promoter cis-element enrichment (ABRE, DRE/CRT, MYB, MYC, W-box) if genome files available.
11. **Functional evidence**: Arabidopsis orthologs of CDS genes cross-referenced with known drought genes.
12. GO/KEGG enrichment with propagated annotations; permutation-based checks.

## Robustness
Sensitivity to subsampling depth, to removal of each study, to signature thresholds; all code and
processed data deposited publicly (GitHub + Zenodo).
