# Conserved core drought transcriptome in watermelon — code and processed data

Code and processed results for:

> Isak, M.A., Şimşek, Ö. A conserved core drought transcriptome in watermelon: multi-study meta-analysis validated
> across tissues, time, stresses and species. Manuscript submitted.

The study uses public RNA-seq data only (91 libraries from 10 BioProjects; see `metadata/design.tsv`). No new
sequencing data were generated.

## Contents

| Folder | Content |
|---|---|
| `metadata/design.tsv` | Sample sheet: run, study, crop, stress, role (discovery / validation / specificity), condition, tissue, genotype, time point, layout, read count, ENA FASTQ URLs, MD5 checksums and sizes. Group labels were assigned from BioSample attributes before analysis. |
| `metadata/all_projects_screened.tsv` | All cucurbit BioProjects retrieved from NCBI for screening (12 October 2026). |
| `metadata/STUDY_PLAN.md` | Pre-specified analysis plan (signature definition and validation tests). |
| `scripts/census/` | BioProject census, run metadata and sample-sheet construction. |
| `scripts/processing/` | Download, depth standardisation and Salmon quantification per run. |
| `scripts/analysis/` | Differential expression, meta-analysis, validation, functional analysis, network, figures and the sensitivity, validation and orthologue-evidence analyses. |
| `results/` | Gene-level count and TPM matrices, per-study DE tables, meta-analysis, Core Drought Signature, validation, function, network and extended results. |

## Software

Python 3.13.16 with the packages pinned in `requirements.txt`; Salmon 1.10.0; DIAMOND 2.1.10; curl, gzip, mawk.
Dynamic Tree Cut is used through its Python port, vendored and patched by `bash scripts/patch_dynamicTreeCut.sh`.

```bash
pip install -r requirements.txt   # dynamicTreeCut line may fail to build: it is vendored by the next command
bash scripts/patch_dynamicTreeCut.sh
export PROJECT_DIR=$(pwd)         # all scripts read and write relative to PROJECT_DIR
export SALMON_BIN=/path/to/salmon-1.10.0/bin
```

## References and indices

| File (in `reference/`) | Source |
|---|---|
| `97103_mRNA_v2.5.fa.gz`, `97103_AHRD_v2.5.txt.gz`, `97103_GO_v2.5.gaf.txt.gz` | CuGenDB v2, watermelon 97103 v2.5 |
| `Citrullus_lanatus.Cla97_v1.dna.toplevel.fa.gz`, `...Cla97_v1.63.chr.gff3.gz`, `...pep.all.fa.gz` | Ensembl Plants (promoters, orthology) |
| `Cucumis_sativus.ASM407v2.cdna.all.fa.gz`, `...pep.all.fa.gz`, `...dna.toplevel.fa.gz`, `...63.gff3.gz`, `...ncrna.fa.gz` | Ensembl Plants |
| `Arabidopsis_thaliana.TAIR10.pep.all.fa.gz` | Ensembl Plants |
| cucumber chloroplast DQ119058; mitochondrion HQ860792–HQ860794 | ENA |
| `at_genes.bed` (TAIR10 gene coordinates) | Ensembl Plants release 63 (`fetch_arabidopsis_public.sh`) |
| `pub/` GSE110079, GSE80564, GSE80565 | NCBI GEO (`fetch_arabidopsis_public.sh`) |

```bash
salmon index -t reference/97103_mRNA_v2.5.fa.gz -i reference/salmon_index -k 25
salmon index -t reference/Cucumis_sativus.ASM407v2.cdna.all.fa.gz -i reference/cu_salmon_index -k 25
```

## Run order

1. **Quantification** (per run; `jobs_*.txt` lines are `RUN READCOUNT URL1;URL2 MD5_1;MD5_2 INDEX BYTES1;BYTES2`):
   `cat jobs.txt | xargs -P 4 -L 1 scripts/processing/proc_fast.sh` → `quant/<run>/quant.sf`.
   Files < 1.2 GB are downloaded completely (MD5-checked) and subsampled to 8 M reads; for larger files the first
   4 M reads / read pairs are used. Robustness runs: `proc_s1.sh` (first 4 M reads of PRJNA454040) and `proc_s2.sh`
   (complete PRJNA326331 files subsampled to 4 M pairs by deterministic record-level selection applied identically
   to both mates).
2. `scripts/analysis/a1_de.py` — gene matrices and PyDESeq2 per study/contrast (`results/de/`, `contrasts.tsv`).
3. `a2_meta.py` — DerSimonian–Laird meta-analysis and Core Drought Signature (`meta.tsv.gz`, `cds.tsv`).
4. `a3_validate.py` — leave-one-study-out, cross-species, specificity, time course, genotype (`results/validation/`).
5. `a4_function.py` — GO, Arabidopsis orthology, transcription factors, promoter motifs (`results/function/`).
6. `a5_network.py` — signed co-expression network, hybrid Dynamic Tree Cut (`results/network/`).
7. Sensitivity, extended validation and orthologue-evidence analyses (`results/extended/`, `results/robustness/`):
   `s1_robust.py`, `e1_subsample_root.py` (subsampling), `e1_meta_sensitivity.py` (estimators, thresholds, study
   exclusion), `e1_lowmap.py` (read composition of low-mapping cucumber libraries), `e1_validation_stats.py`
   (AUC permutation/bootstrap, M_eff, label permutation), `e1b_validation_competitive.py` (competitive null),
   `e1_specificity.py` (salt/cold effect sizes and detectability), `e1_network_stability.py` (per-study, leave-one-out
   and bootstrap network robustness), `e1_candidates.py` (prioritised candidates), `e2_genotype_check.py` (genotype identity of the PRJNA1150737 libraries), `e3_public_functional.py` (public Arabidopsis
   ABA-mutant, ABA time-course and ChIP-seq data for orthologues; run `scripts/census/fetch_arabidopsis_public.sh` first)
   and `e3_candidate_evidence.py` (three-level evidence profile of each signature gene, curated published genetic evidence).
8. `a6_figures.py`, `e3_figure7.py`, `e1_graphical_abstract.py` — Figures 1–7, S1–S3 and the graphical abstract.
9. `e1_supplementary.py` — Supplementary Tables S1–S14 (Excel).

Random seeds are fixed in every script; all numbers in the manuscript are generated from the files in `results/`.

## Mapping of results to the manuscript

| Manuscript | File(s) |
|---|---|
| Table 1, Table S1 | `metadata/design.tsv`, `results/qc.tsv` |
| Fig. 2, Table S2 | `results/de/`, `results/meta.tsv.gz`, `results/cds_annotated.tsv` |
| Fig. 3, Tables S4, S7 | `results/validation/loso.tsv`, `cross_species*.tsv`, `results/extended/validation_uncertainty.tsv`, `validation_competitive.tsv` |
| Fig. 4, Tables S3, S9 | `results/validation/specificity.tsv`, `time_course.tsv`, `tolerance.tsv`, `results/extended/specificity_detectability*.tsv` |
| Fig. 5, Table S5 | `results/function/` |
| Fig. 6, Fig. S1, Tables S5, S10 | `results/network/`, `results/extended/network_*.tsv` |
| Fig. S2, Table S6 | `results/extended/meta_sensitivity.tsv` |
| Fig. S3, Tables S8, S11 | `results/robustness/`, `results/extended/S1b_root_*.tsv`, `lowmap_classification.tsv`, `prjna1150737_genotype_*.tsv` |
| Table S12 | `results/extended/candidate_priority.tsv` |
| Fig. 7, Table S14 | `results/extended/public_functional_sets.tsv`, `public_functional_genes.tsv` |
| Tables 2–3, Table S13 | `results/extended/candidate_evidence.tsv` (literature entries are curated in `e3_candidate_evidence.py`) |

## Notes on the data

* The large tables (`results/tpm_*.tsv.gz`, `results/counts_*.tsv.gz`, `results/meta.tsv.gz`) are gzip-compressed;
  pandas reads them directly, and `gunzip` gives the plain tab-separated files.

* PRJNA1150737 (cucumber) consists mainly of rRNA and chloroplast reads (random-primed total RNA) and its metadata
  are internally inconsistent; groups were taken from the BioSample alias/treatment/cultivar fields, which agree with each other and with shared sequence variants (`e2_genotype_check.py`); the dataset is used as secondary evidence only.
* Two runs of PRJNA1150737 were deposited as three files; only the two mate files were used.
* SRR3706796 duplicates SRR3707382 and was excluded.

## Licence and citation

Code: MIT licence (see `LICENSE`). Processed data: CC BY 4.0. Please cite the article and the archived dataset
(`CITATION.cff`). Raw data belong to the original submitters and are cited in the article.
