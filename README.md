# Epigenomics Pipeline - EN-TEx Data Analysis
**Author:** Alon Shavit  
**Course:** Master in Omics Data Analysis - Epigenomics  
**Instructor:** Beatrice Borsari

## 1. Overview
This repository contains the analysis pipelines and results for Exercises 4 and 5 of the Epigenomics module. The project processes EN-TEx ChIP-seq and ATAC-seq data to characterize the epigenomic landscape of two human tissues (Stomach vs. Sigmoid Colon) using ENCODE data from donor ENCDO451RUA.

## 2. Repository Structure
```
epigenomics_uvic/
├── scripts/                          # Executable Bash pipelines
│   ├── section3.2_ChIP-seq_downstream.sh   # Exercise 4
│   ├── section4_ATAC-seq_downstream.sh     # ATAC-seq analysis
│   └── section5_distal_regulatory_activity.sh  # Exercise 5
├── ChIP-seq/                         # Exercise 4 data & results
│   ├── analyses/
│   │   ├── aggregation.plot/         # H3K4me3 signal over TSS
│   │   ├── scatterplot.correlation/  # Tissue correlation plots
│   │   └── peaks.analysis/           # Gene lists, Venn diagrams
│   └── data/bed.files/               # ENCODE peak files
├── ATAC-seq/                         # ATAC-seq downstream analyses
│   ├── analyses/peaks.analysis/      # Peaks outside gene bodies
│   └── annotation/                   # Gencode v24 annotations
├── regulatory_elements/              # Exercise 5 results
│   ├── *.regulatory_elements.bed     # Candidate enhancers
│   └── *.distances.tsv               # Gene distance calculations
├── bin/                              # Analysis scripts
│   └── get.distance.py               # Nearest-gene calculator
└── README.md
```

---

## 3. Exercise 4: ChIP-seq Downstream Analyses (Section 3.2)

### 3.1 Background
H3K4me3 is a histone modification strongly associated with active promoters and transcription initiation. By integrating ChIP-seq data with RNA-seq expression data, we can confirm this relationship and identify tissue-specific regulatory programs. POLR2A (RNA Polymerase II) ChIP-seq provides complementary evidence of active transcription.

### 3.2 Data
Pre-computed peaks from ENCODE for donor ENCDO451RUA:

| Mark | Stomach | Sigmoid Colon |
|------|---------|---------------|
| H3K4me3 | ENCFF5021JF | ENCFF668DFB |
| POLR2A | ENCFF979WAV | ENCFF477DXO |

### 3.3 Aggregation Plots
H3K4me3 ChIP-seq signal was profiled over the TSS (±2 kb) of the 1,000 most and least expressed genes using deepTools `computeMatrix` and `plotProfile`.

**Stomach:**

![Aggregation Plot Stomach](ChIP-seq/analyses/aggregation.plot/aggregation.plot.stomach.pdf)

**Sigmoid Colon:**

![Aggregation Plot Sigmoid Colon](ChIP-seq/analyses/aggregation.plot/aggregation.plot.sigmoid_colon.pdf)

**Interpretation:** The H3K4me3 signal shows a sharp, bimodal peak directly over the TSS of highly expressed genes in both tissues, with a characteristic dip at the nucleosome-free region at the exact TSS position. The signal is markedly reduced at the TSS of lowly expressed genes. This confirms H3K4me3 as a mark of active promoters — its presence correlates with transcriptional activity. The bimodal pattern reflects the positioned +1 and -1 nucleosomes flanking the TSS.

### 3.4 Scatterplot Correlation
Spearman correlation of H3K4me3 signal between stomach and sigmoid colon is consistently higher than Pearson correlation, reflecting the non-linear, rank-based relationship between epigenomic signals across tissues. Most promoters share similar H3K4me3 patterns, with tissue-specific differences at a subset of loci.

### 3.5 Tissue-Specific Peak Analysis & Gene Ontology
Genes were classified by H3K4me3 peak presence at their TSS:

| Category | Description |
|----------|-------------|
| **Both tissues** | Genes with H3K4me3 in stomach AND sigmoid colon — housekeeping genes |
| **Stomach only** | Genes with H3K4me3 only in stomach — tissue-specific functions |
| **Sigmoid colon only** | Genes with H3K4me3 only in sigmoid colon — tissue-specific functions |
| **Neither** | Genes without H3K4me3 in either tissue — silenced/lowly expressed |

**GO Enrichment Results (Metascape):**
- **Stomach-specific genes** were enriched for gastric acid production pathways (WP2596), validating that H3K4me3 marks promoters of tissue-appropriate functional genes.
- **Sigmoid colon-specific genes** showed enrichment for intestinal absorption and epithelial differentiation pathways.

### 3.6 H3K4me3 and POLR2A Overlap (Venn Diagram)

![Venn Diagram](ChIP-seq/analyses/peaks.analysis/Venn.Diagram.H3K4me3.POLR2A.png)

**Interpretation:** The Venn diagram reveals substantial overlap between H3K4me3 and POLR2A peaks at promoters. Genes marked by both H3K4me3 AND POLR2A represent actively transcribing loci — H3K4me3 maintains the open chromatin state while POLR2A is physically engaged in transcription. Genes with H3K4me3 but without POLR2A may represent poised promoters — accessible and marked but not currently being transcribed. Notably, very few genes show POLR2A without H3K4me3, confirming that H3K4me3 is a prerequisite for stable polymerase recruitment at most promoters.

---

## 4. Exercise 5: Distal Regulatory Activity (Section 5)

### 4.1 Background
Not all open chromatin regions are enhancers. To identify high-confidence active enhancers, we use the classic epigenomic enhancer signature: co-occurrence of H3K27ac (active regulatory elements) and H3K4me1 (enhancer mark). Regions that are simultaneously open (ATAC-seq), H3K27ac-positive, and H3K4me1-positive represent high-confidence active distal regulatory elements.

### 4.2 Histone Mark Data

| Mark | Tissue | Accession | Peaks |
|------|--------|-----------|-------|
| H3K27ac | Stomach | ENCFF736SWI | 57,121 |
| H3K27ac | Sigmoid Colon | ENCFF843KKF | 58,062 |
| H3K4me1 | Stomach | ENCFF193RHT | 152,394 |
| H3K4me1 | Sigmoid Colon | ENCFF442KWF | 95,892 |

### 4.3 ATAC-seq Peak Classification

| Tissue | Total Peaks | Overlapping Promoters | Outside Gene Body |
|--------|------------|----------------------|-------------------|
| Stomach | 103,609 | 44,749 (43.2%) | 34,537 (33.3%) |
| Sigmoid Colon | 110,999 | 47,871 (43.1%) | 37,035 (33.4%) |

### 4.4 Candidate Enhancers (Triple Intersection)

| Tissue | Distal ATAC Peaks | After H3K27ac + H3K4me1 Filter |
|--------|-------------------|-------------------------------|
| Stomach | 34,537 | 12,740 (36.9%) |
| Sigmoid Colon | 37,035 | 16,439 (44.4%) |

### 4.5 Distance to Nearest Gene (Chromosome 1)

| Tissue | Elements (chr1) | Mean Distance | Median Distance |
|--------|----------------|---------------|-----------------|
| Stomach | 1,545 | 47,985 bp | 27,808 bp |
| Sigmoid Colon | 1,734 | 73,019 bp | 35,365 bp |

### 4.6 Biological Interpretation
The median distance is significantly lower than the mean in both tissues, indicating a right-skewed distribution: most enhancers are within ~50 kb of their target genes, but a subset acts over >100 kb distances. This is consistent with enhancer-promoter interactions occurring primarily within topologically associating domains (TADs).

Sigmoid colon shows more regulatory elements (16,439 vs. 12,740) at greater distances from genes (mean: 73 kb vs. 48 kb), suggesting a more distributed regulatory network. This may reflect the higher cellular heterogeneity of colonic tissue, requiring more diverse distal regulatory inputs to control cell-type-specific gene expression programs.

The biological significance of distal regulatory elements lies in their role as enhancers — they physically loop to contact gene promoters (often via cohesin/CTCF-mediated chromatin interactions) and amplify transcription of their target genes. Understanding enhancer-gene distances is critical for interpreting GWAS variants, as many disease-associated SNPs fall in non-coding enhancer regions far from the genes they regulate.

---

## 5. Tools Used
| Tool | Purpose |
|------|---------|
| **BEDTools** | Genomic interval intersection and filtering |
| **deepTools** | computeMatrix, plotProfile for aggregation plots |
| **AWK** | Text processing and coordinate extraction |
| **Python 3** | get.distance.py for nearest-gene computation |
| **R** | Statistical summaries, Venn diagrams, boxplots |
| **Metascape** | Gene Ontology enrichment analysis |
| **ENCODE Portal** | Source of all epigenomic data (GRCh38) |

## 6. Reproducibility
All analysis performed on macOS using pre-computed ENCODE peak files. Files retrievable from https://www.encodeproject.org using accession numbers listed above. Required tools: BEDTools, deepTools, Python 3, R.
