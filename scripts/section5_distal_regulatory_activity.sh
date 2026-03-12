#!/bin/bash
# Section 5: Distal regulatory activity
# Author: Alon Shavit

## Task 1: Create folder
mkdir -p regulatory_elements

## Task 2: Identify candidate regulatory elements
# Intersect distal ATAC-seq peaks with H3K27ac AND H3K4me1

# Stomach: 12,740 candidate enhancers
bedtools intersect \
  -a ATAC-seq/analyses/peaks.analysis/ENCFF762IFP.outside.gene.body.bed \
  -b ChIP-seq/data/bed.files/ENCFF736SWI.bed \
  | bedtools intersect -a stdin \
  -b ChIP-seq/data/bed.files/ENCFF193RHT.bed \
  > regulatory_elements/stomach.regulatory_elements.bed

# Sigmoid colon: 16,439 candidate enhancers
bedtools intersect \
  -a ATAC-seq/analyses/peaks.analysis/ENCFF287UHP.outside.gene.body.bed \
  -b ChIP-seq/data/bed.files/ENCFF843KKF.bed \
  | bedtools intersect -a stdin \
  -b ChIP-seq/data/bed.files/ENCFF442KWF.bed \
  > regulatory_elements/sigmoid_colon.regulatory_elements.bed

## Task 3: Extract chr1 regulatory element starts
# Stomach: 1,545 elements on chr1
awk '$1=="chr1"' regulatory_elements/stomach.regulatory_elements.bed \
  | awk 'BEGIN{FS=OFS="\t"}{print $4, $2}' \
  > regulatory_elements/stomach.regulatory.elements.starts.tsv

# Sigmoid colon: 1,734 elements on chr1
awk '$1=="chr1"' regulatory_elements/sigmoid_colon.regulatory_elements.bed \
  | awk 'BEGIN{FS=OFS="\t"}{print $4, $2}' \
  > regulatory_elements/sigmoid_colon.regulatory.elements.starts.tsv

## Task 4: Extract chr1 gene starts (strand-aware)
# 2,047 protein-coding genes on chr1
# For plus strand: start = column 2 (5' end)
# For minus strand: start = column 3 (3' coordinate = 5' end of gene)
awk 'BEGIN{FS=OFS="\t"}$1=="chr1"{if ($6=="+"){start=$2} else {start=$3}; print $4, start}' \
  ATAC-seq/annotation/gencode.v24.protein.coding.gene.body.bed \
  > regulatory_elements/gene.starts.tsv

## Task 5: Test get.distance.py
# Expected output: ENSG00000187642.9    982093    2093
python3 bin/get.distance.py --input regulatory_elements/gene.starts.tsv --start 980000

## Task 6: Compute distances for all regulatory elements
# Stomach (1,545 elements)
cat regulatory_elements/stomach.regulatory.elements.starts.tsv | while read element start; do
  python3 bin/get.distance.py --input regulatory_elements/gene.starts.tsv --start "$start"
done > regulatory_elements/stomach.regulatoryElements.genes.distances.tsv

# Sigmoid colon (1,734 elements)
cat regulatory_elements/sigmoid_colon.regulatory.elements.starts.tsv | while read element start; do
  python3 bin/get.distance.py --input regulatory_elements/gene.starts.tsv --start "$start"
done > regulatory_elements/sigmoid_colon.regulatoryElements.genes.distances.tsv

## Task 7: Compute mean and median distances using R
# Stomach - mean: 47985.17, median: 27808
cat regulatory_elements/stomach.regulatoryElements.genes.distances.tsv | awk '{print $3}' | \
  Rscript -e 'd=scan("stdin",quiet=TRUE); cat("Stomach mean:", mean(d), "median:", median(d), "\n")'

# Sigmoid colon - mean: 73019.11, median: 35365
cat regulatory_elements/sigmoid_colon.regulatoryElements.genes.distances.tsv | awk '{print $3}' | \
  Rscript -e 'd=scan("stdin",quiet=TRUE); cat("Sigmoid_colon mean:", mean(d), "median:", median(d), "\n")'
