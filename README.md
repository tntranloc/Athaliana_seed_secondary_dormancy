# Master Project Scripts Repository  
**Seed Secondary Dormancy & Associated Manuscript**

**Author:** Nhu L. T. Tran  
**Email:** ntran5@uni-koeln.de  
**Last updated:** 02 September 2024  

---

## Overview
This repository contains scripts developed during my Master’s project on **seed secondary dormancy**, including analyses used in the associated manuscript.  

The workflows focus on:
- Genomic data processing and GWAS  
- Statistical modeling with kinship correction  
- Species distribution (ecological niche) modeling  
- Data visualization and exploratory analysis  

---

## Main Workflows

### Genomics & GWAS
- **`complete_mapping_from_fastq_to_vcf.sh`**  
  End-to-end pipeline for mapping FASTQ files to a reference genome, followed by variant calling and GWAS.

- **`runGWAS_with_statgenGWAS.R`**  
  GWAS pipeline in R using the `statgenGWAS` package.

---

### Statistical Modeling
- **`run_GLM_with_kinshipcovariate.R`**  
  Generalized linear model (GLM) with kinship as a covariate using `lme4` and `lme4qtl`.

---

### Species Distribution Modeling
- **`species_distribution_model_with_biomod2.R`**  
  Workflow for species distribution (ecological niche) modeling using the `biomod2` package.

---

## Supporting Scripts

### Data Exploration & Visualization
- **`multicollinearity_check_and_visulisation.md`**  
  Assess multicollinearity among predictors using Variance Inflation Factor (VIF).

- **`draw_PCA_for_plink_output.R`**  
  Generate PCA plots from `.eigenvec` and `.eigenval` outputs (PLINK).

- **`drawing_nice_sampling_map_ggplot.R`**  
  Create geographical sampling maps with `ggplot2`.

- **`drawing_nice_single_and_multilple_histograms.R`**  
  Functions for plotting single and overlaid histograms.

---

### Data Processing & Utilities
- **`correlation_sliding_window.txt`**  
  R and Python scripts to compute correlations between two variables across sliding windows (user-defined window size and step).

- **`get_bioclimatic_data.R`**  
  Retrieve bioclimatic variables from climate databases.

---

## Associated Publication
This repository accompanies a manuscript on **seed secondary dormancy**.  
Tran NLT, Ali T, Schmitz G, de Meaux J. Heat-Induced Secondary Dormancy Contributes to Local Adaptation in Arabidopsis thaliana. Mol Ecol. 2025 Oct;34(19):e70086. doi: 10.1111/mec.70086. Epub 2025 Aug 26. PMID: 40856109; PMCID: PMC12456118.


---

## Requirements
Workflows rely on:
- **R** (e.g., `lme4`, `lme4qtl`, `biomod2`, `statgenGWAS`, `ggplot2`)  
- **Python** (for selected data processing tasks)  
- **Bash** (for genomic pipelines)  

External tools:
- PLINK  
- Standard bioinformatics tools for mapping and variant calling  

---

## Notes
- Scripts are modular and can be adapted to different species or datasets.  
- Some workflows assume familiarity with genomic data formats (FASTQ, VCF, PLINK outputs).  
- Visualization scripts are designed for publication-quality figures.  
