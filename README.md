# CSE283_Heritability_calculations-using causal TWAS

------ THIS README IS IN THE PROCESS OF BEING UPDATED ------

This repo contains scripts that wrap causal ctwas functions from - Zhao, S., Crouse, W., Qian, S. et al. Adjusting for genetic confounders in transcriptome-wide association studies improves discovery of risk genes of complex traits. Nat Genet 56, 336–347 (2024). https://doi.org/10.1038/s41588-023-01648-9

**The goal of this project is to critically compare two cTWAS modalities for calculating heritability**

1) Single tissue cTWAS
   In single cTWAS, each gene contributes one predicted expression valriable. The model estimates a single effect size BETA-gene and corresponding PIP for that gene using
   the association between GWAA signal and predicted expression / molecular trait for that tissue.
   
2) Multi tissue cTWAS
    In multi tissue model, multiple predicted expression variables are calculated for all processed tissues. These will be correlated if genes across tissues share the 
    same / similar QTLs. cTWAS will therefore jointly estimate effects of expression variables from all tissues whilst accounting for correlation structure stored in a covariance matrix. cTWAS will shrink the effect sizes towards a shared prior and allocates posterior probability (measure of causality) across tissues.

An **expression variable** is a measure of the **expression level of a gene that can be exaplined by genetics**.

For gene \( g \):

\[
\hat{E}_g = \sum_{i=1}^{m} w_{gi} G_i
\]

where

- \( G_i \) = genotype of SNP \( i \) (often standardized)  
- \( w_{gi} \) = prediction weight learned from an eQTL model (e.g., PrediXcan / mashR)  
- \( m \) = number of SNPs in the prediction model  

This predicted value \( \hat{E}_g \) represents the **component of gene expression determined by genetics**.

In **cTWAS / TWAS**, this predicted expression acts as a **gene-level predictor variable** in the association model, allowing the method to test whether **genetically regulated expression of a gene is associated with the trait**.

#############################
#############################
**Heritability as a readout**
To be added 
— this section will describe how cTWAS estimates the proportion of trait variance explained by genetically predicted gene expression (gene-level PVE) relative to total trait heritability, providing a quantitative measure of how much disease signal can be attributed to regulatory mechanisms captured by expression prediction models.
#############################
#############################

**Biological Background**

**CSE283**
In CSE 283 we have studied population genetics, linkage disequilibrium (LD), and statistical assocaition models to assess how genetic variants (mainly signle nucleotide polymorphisms) are linked to complex traits and disease.

**The Problem**
One missing peice of GWAS is that the association testing does not give insight into the causal tissue/s or cell type/s through which detected variants are acting.
Such information is critical for understanding disease pathology and etiology.

The mechanisms through which genetic variants act can be broadly categorised into -
1) variants falling into coding regions and alter RNA code (and protein structure for protein coding genes) 
2) variants falling iside non-coding regions of the genome, impacting regulatory elements.

It is well known that over 90% of GWAS hits fall in non-coding regions of the genome (Maurano et al., 2012; ENCODE Project Consortium, 2012), implying that most signals likely act through gene regulation rather than protein-coding changes. 
Regulatory activity is highly tissue- and cell-type-specific (GTEx Consortium, 2015; 2020).
Therefore, without integrating tissue-resolved functional data, GWAS loci cannot directly reveal where the biological mechanism operates.

**QTL Studies**
To start narrowing down causality of gwas variants, we can measure molcular phenotypes (e.g. mRNA levels, splice junctions, protein abundance) and stratify molecular phenotype measurements by genotype to identify genotype–molecular phenotype associations in what are commonly known as quantitative trait loci (QTL) studies.

**TWAS Studies**
These associations can be used to construct genetic prediction models of molecular traits, which form the basis of transcriptome-wide association studies (TWAS) and related methods that link genetically predicted molecular variation to complex disease risk.

CTWAS - ...
(Ignore. Methods such as TWAS (Gusev et al., 2016), colocalization (Giambartolomei et al., 2014), and cTWAS.)

Heritability - 


**------ Set Up Instructions ------**

This file consists of wrapper scripts to process using cTWAS tool.
Three inputs are required - 
1) GWAS DATA -
   (Must be in vcf, tsv or csv format)
   For readily availabe GWAS refer to ieu open gwas project - 
   https://opengwas.io/

2) Pretrained weights from Prediction models e.g. PrediXcan (or QTL summary statistics)
   In .db format as prepared by PrediXcan.
   For readily available prediction models please refer to:
   GTEx v8 models on eQTL and sQTL 
   https://predictdb.org/post/2021/07/21/gtex-v8-models-on-eqtl-and-sqtl/

3) LD matrices
   Please note LD structure is critical for cTWAS. Due to imputation of GWAS weights from
   prediction model statistics, in addition to causality being tested per habloblock,
   it is important to match LD information to GWAS data.
   (GWAS, model and LD should be matched on a pop level.)
   For EUR LD matrix please refer to -
   https://uchicago.app.box.com/s/jqocacd2fulskmhoqnasrknbt59x3xkn

**------ Environment Set Up: ------**

If not already please install mamba or conda onto your local / cluster account:
https://github.com/conda-forge/miniforge

subsequently run the following command to generate SnakeMake and ctwas environments from yaml files.

mamba env create -f environment.yml

Environment yamls can be found in the envs folder.

**------ CTWAS Wrapper scripts: ------**

These scripts can be run in command line 

**1_prepare_referece**
Constructs the reference data required for cTWAS by extracting variants from the LD reference panel, 
harmonizing variant identifiers, and computing or subsetting LD matrices for the variants used in the analysis. 
It ensures that SNPIDs, allele orientation, and genomic coordinates are consistent across reference, weights, and GWAS inputs.

**Script arguments**

- **`--region_path`**  
  Path to the genomic region definition file used to partition the genome into LD regions for cTWAS analysis.

- **`--ld_dir`**  
  Directory containing LD reference data (e.g., SNP covariance/correlation matrices) used to compute variant correlation structure.

- **`--genome_version`** *(default: `b38`)*  
  Genome build used for variant coordinates and region definitions (e.g., `b37` or `b38`).

- **`--chromosome`** *(optional)*  
  Integer specifying a chromosome to subset the analysis to a single chromosome instead of processing all regions.

- **`--outpath`**  
  Output directory where processed reference objects (typically `.RDS` files) will be written for downstream cTWAS steps.

**2_Process_GWAS**
Pre-precess GWAS summary statistics to produce / extract GWAS Z-score from GWAS summary statistics.
This step harmonises GWAS SNPIDs with the LD reference panel, checks allele orientation, removes problematic variants
and formats the GWAS summary statistics for running cTWAS.

**Script arguments**

- **`--GWAS_Path`**  
  Path to the GWAS summary statistics file that will be processed and harmonized for use in the cTWAS pipeline.

- **`--SNP_map_path`**  
  Path to the SNP reference map used to align GWAS variants with the reference panel (e.g., mapping rsIDs or variant IDs to reference coordinates).

- **`--genome_version`** *(default: `b38`)*  
  Genome build used for variant coordinates (e.g., `b37` or `b38`) to ensure consistency between GWAS data and the reference panel.

- **`--chromosome`** *(optional)*  
  Integer specifying a chromosome to subset the analysis to a single chromosome rather than processing the full GWAS dataset.

- **`--outpath`**  
  Directory where the processed GWAS objects (typically `.RDS` files) will be written for downstream cTWAS steps.

- **`--disease`** *(optional)*  
  Label or identifier for the trait/disease being analyzed, typically used for naming output files or organizing results.


**3_process_weights**
Pre-process prediction model weights (e.g. PrediXcan - elastic net or mashR models) by harmonising SNPIDs and ensure that 
LD information is available for each SNP. SNPs are then harmonised with LD panels and GWAS.

**Script arguments**

- **`--weight_files`**  
  Path to gene expression prediction model weight files (e.g., PrediXcan/mashR models) used to construct genetically predicted expression variables.

- **`--region_path`**  
  Path to the genomic region definition file used to assign variants and gene models to LD regions.

- **`--chromosome`** *(optional)*  
  Integer specifying a chromosome to subset the analysis to a single chromosome instead of processing all regions.

- **`--snp_map_path`**  
  Path to the SNP reference map used to harmonize variant identifiers between weight models, GWAS data, and the LD reference.

- **`--gwas_snpid_path`**  
  Path to the file containing SNP identifiers from the processed GWAS dataset, used to ensure consistency between GWAS variants and prediction model SNPs.

- **`--prediction_model`** *(optional)*  
  Identifier for the prediction model being used (e.g., specific tissue or model set), typically used for labeling outputs.

- **`--type`** *(optional)*  
  QTL type associated with the prediction model (e.g., `eQTL`, `sQTL`, `pQTL`).

- **`--context`** *(optional)*  
  Biological context for the prediction model, commonly indicating the tissue or cell type.

- **`--LD_map_path`** *(optional)*  
  Path to a table containing filenames for LD matrices and associated SNP information for each genomic region. Required when LD is loaded externally rather than from PredictDB resources.

- **`--outpath`**  
  Output directory where processed weight objects and associated metadata (typically `.RDS` files) will be written.

- **`--ncores`** *(default: `1`)*  
  Number of CPU cores used for parallel processing during weight preparation.


**4_cTWAS_Runner**
Runs the core cTWAS analysis, integrating GWAS summary statistics, LD reference data, and gene expression prediction weights and applying a fine-mapping utilising a Bayesian statistical approach. 
**This method jointly models SNP and gene effects** to estimate posterior inclusion probabilities (PIPs) and effect size distributions for both variant and gene features.


**Script arguments**

- **`--harmonised_gwas_score`**  
  Path to the harmonized GWAS summary statistics containing standardized association statistics (e.g., Z-scores) aligned to the reference SNP set.

- **`--harmonised_weights`**  
  One or more paths to processed prediction model weight files that have been harmonized with the reference SNP map. Multiple files can be supplied (e.g., for multiple tissues or conditions).

- **`--region_path`**  
  Path to the genomic region definition file used to assign variants and gene models to LD regions for cTWAS analysis.

- **`--snp_map`**  
  Path to the SNP reference map used to align variant identifiers between GWAS data, prediction models, and LD reference panels.

- **`--LD_map`**  
  Path to a table containing filenames and metadata for LD matrices corresponding to each genomic region.

- **`--maxSNP`** *(default: `20000`)*  
  Maximum number of SNPs allowed within a genomic region during analysis to control computational complexity.

- **`--min_group_size`** *(default: `100`)*  
  Minimum number of variants required within a group or region for it to be included in the cTWAS analysis.

- **`--ncore`** *(default: `6`)*  
  Number of CPU cores used for parallel processing.

- **`--genome_version`** *(default: `b38`)*  
  Genome build used for variant coordinates (e.g., `b37` or `b38`).

- **`--chromosome`** *(optional)*  
  Integer specifying a chromosome to subset the analysis to a single chromosome.

- **`--outpath`**  
  Directory where cTWAS result objects (typically `.RDS` files) will be written.

- **`--fname`**  
  Optional filename prefix used when saving output results.


**5_summarise**
Aggregates and formats cTWAS outputs across regions or tissues into analysis-ready tables.
This step typically combines PIPs, effect size summaries, and variance explained metrics (e.g., PVE), producing final result files for downstream interpretation, comparison, and visualization.

**Script arguments**

- **`--finemap_path`**  
  One or more paths to cTWAS fine-mapping result files (e.g., SuSiE output tables) that will be aggregated during the summarization step.

- **`--param_path`**  
  One or more paths to parameter files containing estimated model parameters from the cTWAS runs (e.g., priors, variance estimates).

- **`--snp_map_path`** *(optional)*  
  Path to the SNP reference map used to annotate or harmonize SNP identifiers when generating the final summary tables.

- **`--gwas_n`** *(optional)*  
  GWAS sample size used to compute or rescale downstream metrics (e.g., variance explained or effect size summaries).

- **`--outpath`**  
  Directory where the aggregated summary outputs (e.g., combined tables or `.RDS` files) will be written.




**------ Running SnakeFile ------**

Snakefile  - for SLURM server - 
configs


**On going work**

