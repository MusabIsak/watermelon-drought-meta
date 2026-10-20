#!/usr/bin/env bash
# Public Arabidopsis data used by e3_public_functional.py (orthologue-level evidence, Section 2.8).
set -euo pipefail
cd "${PROJECT_DIR:-.}"
mkdir -p pub/gse110079 pub/song2016 ref
G=https://ftp.ncbi.nlm.nih.gov/geo/series
# ABA-deficient mutants aba2-2/aba3-1 vs Col-0, mock vs PEG (Watanabe et al., 2018)
curl -sS -o pub/gse110079/matrix.txt.gz $G/GSE110nnn/GSE110079/matrix/GSE110079_series_matrix.txt.gz
# ABA time course RNA-seq (Song et al., 2016)
curl -sS -o pub/song2016/GSE80565_RNAseq_ABAtimeSeries_rpkm_log2FC_FDR_filterCPM2n2.txt.gz \
  $G/GSE80nnn/GSE80565/suppl/GSE80565_RNAseq_ABAtimeSeries_rpkm_log2FC_FDR_filterCPM2n2.txt.gz
# ChIP-seq peaks of 21 ABA-responsive TFs, ABA condition (Song et al., 2016)
curl -sS "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80564&targ=self&form=text&view=brief" \
  | grep Series_supplementary | grep _ABA_ | sed 's/.*= ftp:/https:/' | tr -d '\r' \
  | while read -r u; do (cd pub/song2016 && curl -sS -O "$u"); done
# TAIR10 gene coordinates (Ensembl Plants release 63)
curl -sS https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-63/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3.gz \
  | gzip -dc | awk -F'\t' '$3=="gene"{split($9,a,";"); print $1"\t"$4"\t"$5"\t"$7"\t"a[1]}' | sed 's/ID=gene://' > ref/at_genes.bed
