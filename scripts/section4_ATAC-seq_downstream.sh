#!/bin/bash
# Section 4: EN-TEx ATAC-seq downstream analyses
# Author: Alon Shavit

## Task 1: Create folder structure
mkdir -p ATAC-seq/analyses/peaks.analysis
mkdir -p ATAC-seq/data/bigBed.files
mkdir -p ATAC-seq/data/bed.files

## Task 2: Download ATAC-seq peaks from ENCODE
# Peaks were downloaded as bigBed files and converted to BED format
# ENCFF287UHP - sigmoid_colon (110,999 peaks)
# ENCFF762IFP - stomach (103,609 peaks)

## Task 3: Intersection analysis with BEDTools

# Peaks overlapping promoter/TSS regions
# Sigmoid colon: 47,871
bedtools intersect -a ATAC-seq/data/bed.files/ENCFF287UHP.bed \
  -b ATAC-seq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed \
  -u | wc -l

# Stomach: 44,749
bedtools intersect -a ATAC-seq/data/bed.files/ENCFF762IFP.bed \
  -b ATAC-seq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed \
  -u | wc -l

# Peaks outside gene body coordinates
# Sigmoid colon: 37,035
bedtools intersect -a ATAC-seq/data/bed.files/ENCFF287UHP.bed \
  -b ATAC-seq/annotation/gencode.v24.protein.coding.gene.body.bed \
  -v > ATAC-seq/analyses/peaks.analysis/ENCFF287UHP.outside.gene.body.bed

# Stomach: 34,537
bedtools intersect -a ATAC-seq/data/bed.files/ENCFF762IFP.bed \
  -b ATAC-seq/annotation/gencode.v24.protein.coding.gene.body.bed \
  -v > ATAC-seq/analyses/peaks.analysis/ENCFF762IFP.outside.gene.body.bed
