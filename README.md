# CSE283_Heritability_calculations-using causal TWAS

------ THIS README IS NOT COMPLETE YET ------

Methods such as TWAS (Gusev et al., 2016), colocalization (Giambartolomei et al., 2014), and cTWAS.

This repo contains scripts that wrap causal ctwas functions from - Zhao, S., Crouse, W., Qian, S. et al. Adjusting for genetic confounders in transcriptome-wide association studies improves discovery of risk genes of complex traits. Nat Genet 56, 336–347 (2024). https://doi.org/10.1038/s41588-023-01648-9

------ Set Up Instructions ------

------ Environment Set Up: ------

SnakeMake env

ctwas env 

------ CTWAS Wrapper scripts: ------

1_prepare_referece
2_Process_GWAS
3_process_weights
4_cTWAS_Runner
5_summarise

------ Running SnakeFile ------

Snakefile  - for SLURM server - 
configs




Background 

In CSE 283 we have studied population genetics, linkage disequilibrium (LD), and statistical assocaition models to assess how genetic variants (mainly signle nucleotide polymorphisms) are linked to complex traits and disease.

One fundamental limitation of GWAS is that the association testing does not give insight into the causal tissue/s or cell type/s through which detected variants are acting. Such information is critical in understanding pathology and disease etiology.

The mechanisms through which genetic variants act can be broadly categorised into 1) variants falling into coding regions and alter RNA code (and protein structure for protein coding genes) 2) variants falling iside non-coding regions of the genome, impacting regulatory elements.

It is well known that over 90% of GWAS hits fall in non-coding regions of the genome (Maurano et al., 2012; ENCODE Project Consortium, 2012), implying that most signals likely act through gene regulation rather than protein-coding changes. Regulatory activity is highly tissue- and cell-type-specific (GTEx Consortium, 2015; 2020). Therefore, without integrating tissue-resolved functional data, GWAS loci cannot directly reveal where the biological mechanism operates.

