Authors: Nhu L.T.Tran

Email: ntran5@uni-koeln.de

Last updated: 02.09.2024

Collected scripts for my Master project and the manuscript thereof in Seed secondary dormancy

Main scripts:
- complete_mapping_from_fastq_to_vcf.sh shows how to map fastq files to a reference genome, followed by variant calling and GWAS.
- run_GLM_with_kinshipcovariate.R shows how to run GLM with kinship as covariate using lme4 and lme4qtl packages
- species_distribution_model_with_biomod2.R shows how to make species distribution model, or aka ecological niche model, using biomod2 package

Supporting scripts:
- multicollinearity_check_and_visulisation.md shows how to check for multicollinearity of predictor variables using Variance Inflation Factor (VIF) values.
- draw_PCA_for_plink_output.R shows how to make nice PCA from .eigenvec and .eigenval output from PLINK.
- correlation_sliding_window.txt uses R and Python to compute correlation between 2 variables through many sliding windows, one can specify window size and step of sliding.
- drawing_nice_sampling_map_ggplot.R shows how to draw geographical map with data points of interest
- drawing_nice_single_and_multilple_histograms.R have functions to draw single or overlayed histograms
- get_bioclimatic_data.R shows how to get bioclimatic variables from climate database
- runGWAS_with_statgenGWAS.R shows how to run GWAS in R usingf statgenGWAS package
