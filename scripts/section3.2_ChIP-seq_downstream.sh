#!/bin/bash
# Section 3.2: EN-TEx ChIP-seq downstream analyses (Exercise 4)
# Author: Alon Shavit

#############################
## Data Setup
#############################
# ChIP-seq peaks from ENCODE (donor ENCDO451RUA, GRCh38):
# ENCFF668DFB - sigmoid_colon H3K4me3
# ENCFF5021JF - stomach H3K4me3
# ENCFF477DXO - sigmoid_colon POLR2A
# ENCFF979WAV - stomach POLR2A

#############################
## Aggregation Plots
#############################
# Generate signal matrices over TSS of most/least expressed genes
# using deepTools computeMatrix and plotProfile

# Stomach H3K4me3
computeMatrix reference-point -S ChIP-seq/data/bigWig.files/*.stomach.H3K4me3.bigWig \
  -R ChIP-seq/analyses/stomach.1000.most.expressed.genes.txt \
     ChIP-seq/analyses/stomach.1000.least.expressed.genes.txt \
  --referencePoint TSS -a 2000 -b 2000 \
  -o ChIP-seq/analyses/aggregation.plot/stomach.H3K4me3.matrix.gz

plotProfile -m ChIP-seq/analyses/aggregation.plot/stomach.H3K4me3.matrix.gz \
  -o ChIP-seq/analyses/aggregation.plot/aggregation.plot.stomach.pdf \
  --perGroup

# Sigmoid colon H3K4me3
computeMatrix reference-point -S ChIP-seq/data/bigWig.files/*.sigmoid_colon.H3K4me3.bigWig \
  -R ChIP-seq/analyses/sigmoid_colon.1000.most.expressed.genes.txt \
     ChIP-seq/analyses/sigmoid_colon.1000.least.expressed.genes.txt \
  --referencePoint TSS -a 2000 -b 2000 \
  -o ChIP-seq/analyses/aggregation.plot/sigmoid_colon.H3K4me3.matrix.gz

plotProfile -m ChIP-seq/analyses/aggregation.plot/sigmoid_colon.H3K4me3.matrix.gz \
  -o ChIP-seq/analyses/aggregation.plot/aggregation.plot.sigmoid_colon.pdf \
  --perGroup

#############################
## Scatterplot Correlation
#############################
# Pearson and Spearman correlation of H3K4me3 signal between tissues
# using deepTools multiBigwigSummary and plotCorrelation

multiBigwigSummary bins \
  -b ChIP-seq/data/bigWig.files/*.H3K4me3.bigWig \
  -o ChIP-seq/analyses/H3K4me3.matrix.tsv

plotCorrelation -in ChIP-seq/analyses/H3K4me3.matrix.tsv \
  --corMethod spearman --whatToPlot scatterplot \
  -o ChIP-seq/analyses/scatterplot.correlation/scatterplot.correlation.stomach.pdf

plotCorrelation -in ChIP-seq/analyses/H3K4me3.matrix.tsv \
  --corMethod spearman --whatToPlot scatterplot \
  -o ChIP-seq/analyses/scatterplot.correlation/scatterplot.correlation.sigmoid_colon.pdf

#############################
## Peak Analysis - Gene Lists
#############################
# Identify genes with tissue-specific H3K4me3 peaks

# Genes with H3K4me3 peaks in stomach
bedtools intersect -a ChIP-seq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed \
  -b ChIP-seq/data/bed.files/ENCFF5021JF.bed -u \
  | awk '{print $4}' | sort -u > ChIP-seq/analyses/peaks.analysis/genes.with.peaks.stomach.H3K4me3.txt

# Genes with H3K4me3 peaks in sigmoid colon
bedtools intersect -a ChIP-seq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed \
  -b ChIP-seq/data/bed.files/ENCFF668DFB.bed -u \
  | awk '{print $4}' | sort -u > ChIP-seq/analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.H3K4me3.txt

# Genes with POLR2A peaks in stomach
bedtools intersect -a ChIP-seq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed \
  -b ChIP-seq/data/bed.files/ENCFF979WAV.bed -u \
  | awk '{print $4}' | sort -u > ChIP-seq/analyses/peaks.analysis/genes.with.peaks.stomach.POLR2A.txt

# Genes with POLR2A peaks in sigmoid colon
bedtools intersect -a ChIP-seq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed \
  -b ChIP-seq/data/bed.files/ENCFF477DXO.bed -u \
  | awk '{print $4}' | sort -u > ChIP-seq/analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.POLR2A.txt

# Genes marked in both tissues
comm -12 \
  ChIP-seq/analyses/peaks.analysis/genes.with.peaks.stomach.H3K4me3.txt \
  ChIP-seq/analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.H3K4me3.txt \
  > ChIP-seq/analyses/peaks.analysis/genes.marked.both.tissues.H3K4me3.txt

# Genes specific to sigmoid colon
comm -13 \
  ChIP-seq/analyses/peaks.analysis/genes.with.peaks.stomach.H3K4me3.txt \
  ChIP-seq/analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.H3K4me3.txt \
  > ChIP-seq/analyses/peaks.analysis/genes.specific.sigmoid_colon.H3K4me3.txt

# Genes not marked in either tissue
comm -23 \
  ChIP-seq/analyses/peaks.analysis/universe.genes.txt \
  <(cat ChIP-seq/analyses/peaks.analysis/genes.with.peaks.stomach.H3K4me3.txt \
        ChIP-seq/analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.H3K4me3.txt \
    | sort -u) \
  > ChIP-seq/analyses/peaks.analysis/genes.not.marked.H3K4me3.txt

#############################
## Venn Diagram
#############################
# Overlap between H3K4me3 and POLR2A peaks
# Generated using R VennDiagram package
# Output: ChIP-seq/analyses/peaks.analysis/Venn.Diagram.H3K4me3.POLR2A.png

#############################
## Boxplot Expression
#############################
# Compare expression levels of genes with vs without H3K4me3 peaks
# Output: ChIP-seq/analyses/peaks.analysis/boxplot.expression.pdf

#############################
## Metascape GO Enrichment
#############################
# Tissue-specific gene lists submitted to Metascape (https://metascape.org)
# for Gene Ontology enrichment analysis
# Results: stomach-specific genes enriched for gastric acid production (WP2596)
