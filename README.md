# Epigenomics Pipeline - EN-TEx Data Analysis
**Author:** Alon Shavit
**Course:** Master in Omics Data Analysis - Epigenomics
**Instructor:** Beatrice Borsari

## 1. Overview
This repository contains the analysis pipelines and results for Exercises 4 and 5 of the Epigenomics module. The pipeline processes EN-TEx ATAC-seq and ChIP-seq data to identify distal regulatory elements (putative enhancers) in human tissues (Stomach vs. Sigmoid Colon) using ENCODE data from donor ENCDO451RUA.

## 2. Repository Structure
* `ATAC-seq/` — Section 4: ATAC-seq downstream analyses
  * `data/bed.files/` — ENCODE peak files (ENCFF287UHP, ENCFF762IFP)
  * `data/bigBed.files/` — Original bigBed format files with md5sum verification
  * `analyses/peaks.analysis/` — BEDTools intersection results and answers
  * `annotation/` — Gencode v24 gene body and TSS annotations
* `regulatory_elements/` — Section 5: Distal regulatory activity
  * Candidate enhancer BED files for both tissues
  * Gene start coordinates and distance calculations (chr1)
* `bin/` — Scripts including completed `get.distance.py`

## 3. Section 4: EN-TEx ATAC-seq Downstream Analyses

### 3.1 Background
ATAC-seq (Assay for Transposase-Accessible Chromatin with sequencing) identifies regions of open chromatin across the genome. Open chromatin regions correspond to active regulatory elements — promoters, enhancers, and insulators — where transcription factors can bind. By comparing ATAC-seq peaks against gene annotations, we can classify these regions by their regulatory function.

### 3.2 Data
ATAC-seq pseudoreplicated peaks (bigBed narrow, GRCh38) were downloaded from ENCODE for stomach (ENCFF762IFP: 103,609 peaks) and sigmoid colon (ENCFF287UHP: 110,999 peaks). Files were verified using md5sum checksums against ENCODE-provided values.

### 3.3 Intersection Analysis
Using BEDTools intersect, peaks were classified by their genomic context relative to Gencode v24 protein-coding gene annotations:

| Tissue | Total Peaks | Overlapping Promoters (TSS) | Outside Gene Body |
|--------|------------|----------------------------|-------------------|
| Stomach | 103,609 | 44,749 (43.2%) | 34,537 (33.3%) |
| Sigmoid Colon | 110,999 | 47,871 (43.1%) | 37,035 (33.4%) |

### 3.4 Biological Interpretation
The distribution of ATAC-seq peaks across genomic features is remarkably consistent between the two tissues: approximately 43% of peaks overlap promoter regions (TSS), while about 33% fall outside gene bodies entirely. The remaining ~24% overlap gene bodies but not promoters, likely representing intronic enhancers or other intragenic regulatory elements.

The high proportion of peaks at promoters is expected — promoters of actively transcribed genes must maintain open chromatin to allow transcription factor binding and RNA Polymerase II recruitment. The significant fraction of peaks outside gene bodies (~34,000-37,000 per tissue) represents candidate distal regulatory elements, primarily enhancers, which are the focus of Section 5.

Sigmoid colon shows slightly more ATAC-seq peaks overall (110,999 vs. 103,609), suggesting a somewhat more complex regulatory landscape, which could reflect the higher cellular heterogeneity of colonic tissue compared to gastric epithelium.

## 4. Section 5: Distal Regulatory Activity

### 4.1 Background
Not all open chromatin regions are enhancers. To identify high-confidence active enhancers, we use the classic epigenomic enhancer signature: the co-occurrence of H3K27ac (marks active regulatory elements) and H3K4me1 (marks both active and poised enhancers). Regions that are simultaneously open (ATAC-seq), H3K27ac-positive, and H3K4me1-positive represent high-confidence active distal regulatory elements.

### 4.2 Histone Mark Data
H3K27ac and H3K4me1 pseudoreplicated peaks were downloaded from ENCODE for the same donor (ENCDO451RUA):

| Mark | Tissue | Accession | Experiment | Peaks |
|------|--------|-----------|------------|-------|
| H3K27ac | Stomach | ENCFF736SWI | ENCSR204OJS | 57,121 |
| H3K27ac | Sigmoid Colon | ENCFF843KKF | ENCSR937EVN | 58,062 |
| H3K4me1 | Stomach | ENCFF193RHT | ENCSR158WBG | 152,394 |
| H3K4me1 | Sigmoid Colon | ENCFF442KWF | ENCSR775LGE | 95,892 |

### 4.3 Identifying Candidate Enhancers
Distal ATAC-seq peaks (outside gene body) were sequentially intersected with H3K27ac and H3K4me1 peaks using BEDTools:

| Tissue | Distal ATAC Peaks | After H3K27ac + H3K4me1 Filter |
|--------|-------------------|-------------------------------|
| Stomach | 34,537 | 12,740 (36.9%) |
| Sigmoid Colon | 37,035 | 16,439 (44.4%) |

Approximately 37-44% of distal open chromatin regions carry the active enhancer signature, while the remainder may represent poised enhancers (H3K4me1 only), primed regions, or other non-enhancer regulatory elements such as insulators.

### 4.4 Distance to Nearest Gene (Chromosome 1)
The `get.distance.py` script was completed to compute the distance from each regulatory element to its closest protein-coding gene on chromosome 1 (2,047 genes).

| Tissue | Regulatory Elements (chr1) | Mean Distance | Median Distance |
|--------|---------------------------|---------------|-----------------|
| Stomach | 1,545 | 47,985 bp | 27,808 bp |
| Sigmoid Colon | 1,734 | 73,019 bp | 35,365 bp |

### 4.5 Biological Interpretation
The median distance is significantly lower than the mean in both tissues, indicating a right-skewed distribution: most regulatory elements are relatively close to their target genes (<50 kb), but a subset of distal enhancers acts over much larger genomic distances (>100 kb). This is consistent with the known genomic architecture where most enhancer-promoter interactions occur within topologically associating domains (TADs).

Sigmoid colon regulatory elements are on average farther from genes than stomach (mean: 73 kb vs. 48 kb), which could reflect differences in chromatin organization and enhancer-gene regulatory architecture between these tissues. The higher number of regulatory elements in sigmoid colon (16,439 vs. 12,740) combined with their greater distance from genes suggests a more distributed regulatory network in colonic tissue.

## 5. Tools Used
* **BEDTools** — genomic interval intersection and filtering
* **AWK** — text processing and coordinate extraction
* **Python 3** — get.distance.py script for nearest-gene computation
* **R** — statistical summary (mean/median calculations)
* **ENCODE Portal** — source of all epigenomic data (GRCh38 assembly)

## 6. Reproducibility
All analysis was performed on macOS. ENCODE files can be retrieved from https://www.encodeproject.org using the accession numbers listed above. BEDTools must be installed and available in PATH.
