# CSE283_Heritability_calculations-using causal TWAS

------ THIS README IS IN THE PROCESS OF BEING UPDATED ------

Methods such as TWAS (Gusev et al., 2016), colocalization (Giambartolomei et al., 2014), and cTWAS.

This repo contains scripts that wrap causal ctwas functions from - Zhao, S., Crouse, W., Qian, S. et al. Adjusting for genetic confounders in transcriptome-wide association studies improves discovery of risk genes of complex traits. Nat Genet 56, 336–347 (2024). https://doi.org/10.1038/s41588-023-01648-9

**The goal of this project is to critically compare two cTWAS modalities for calculating heritability**
1) Single tissue cTWAS
   
2) Multi tissue cTWAS


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

Heritability - 


------ Set Up Instructions ------
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

------ Environment Set Up: ------

If not already please install mamba or conda onto your local / cluster account:
https://github.com/conda-forge/miniforge

subsequently run the following command to generate SnakeMake and ctwas environments from yaml files.

mamba env create -f environment.yml

Environment yamls can be found in the envs folder.

------ CTWAS Wrapper scripts: ------
These scripts can be run in command line 

**1_prepare_referece**
Constructs the reference data required for cTWAS by extracting variants from the LD reference panel, 
harmonizing variant identifiers, and computing or subsetting LD matrices for the variants used in the analysis. 
It ensures that SNPIDs, allele orientation, and genomic coordinates are consistent across reference, weights, and GWAS inputs.

args - 

**2_Process_GWAS**
Pre-precess GWAS summary statistics to produce / extract GWAS Z-score from GWAS summary statistics.
This step harmonises GWAS SNPIDs with the LD reference panel, checks allele orientation, removes problematic variants
and formats the GWAS summary statistics for running cTWAS.

args - 


**3_process_weights**
Pre-process prediction model weights (e.g. PrediXcan - elastic net or mashR models) by harmonising SNPIDs 

args - 


**4_cTWAS_Runner**
Runs the core cTWAS analysis, integrating GWAS summary statistics, LD reference data, and gene expression prediction weights within a Bayesian fine-mapping framework. 
The method jointly models SNP and gene effects to estimate posterior inclusion probabilities (PIPs) and effect size distributions for both variant and gene features.

args - 

**5_summarise**
Aggregates and formats cTWAS outputs across regions or tissues into analysis-ready tables.
This step typically combines PIPs, effect size summaries, and variance explained metrics (e.g., PVE), producing final result files for downstream interpretation, comparison, and visualization.
args - 

------ Helper Script ------


------ Running SnakeFile ------

Snakefile  - for SLURM server - 
configs


**On going work**

