 Causal Transcriptome Wide Association Study to calculate heritability of GWAS variants associated with gene expression (option 2)

 Johannes Ottensmeier
 
 Project Description and Methodology – This project will assess two modalities of transcriptome-wide association using the causal-TWAS (cTWAS) tool. https://doi.org/10.1038/s41588-023-01648-9
 
 1.	Single-tissue cTWAS, where gene expression prediction models from one tissue are used to infer gene-level associations
 2.	Multi-tissue cTWAS, where models across multiple tissues are jointly modeled to estimate shared and tissue-specific effects
    
    The output of these approaches calculates GWAS heritability on both a variant-level and a gene-level. The methodologies perform a combined transcriptome wide-association study followed by fine-mapping to calculate gene and variant level posterior inclusion probabilities (PIPs) and proportion of variance explained (PVE).

 1. Run cTWAS using GTEx v8 mashr-based prediction models trained on expression data.
 2. Perform cTWAS on both single-tissue and multi-tissue models (~10 tissues) to 2-3 GWAS summary statistics producing Gene-level PIPs (from SuSiE fine-mapping, gene-level PVE and total PVE for heritability to subsequently calculate proportion heritability explained.
 
 (code written and implemented into snake-make system for single tissue – need to finish for multi tissue).

 3. Determine concordance between tissue and modeling approaches by –
 
 a)	Correlation of of gene PIPs.

 b)	Correlation of gene-level PVE / total PVE,

 c)	Overlap among top 1%, 5% and 10% of genes.

 d)	Rank consistency of top prioritized genes.

 **Dataset Description**

 The primary dataset will consist of UK Biobank GWAS summary statistics for a complex trait (e.g., LDL - ukb-d-30780_irnt, IBD - ebi-a-GCST004131 and SBP - ukb-a-360) – at hand. 
 (I note that these data are in hg37. However, harmonisation is conducted between all datasets to ensure SNPID consistency (on RSID) considering REF-ALT also. Please refer to this tutorial)
 Gene expression prediction models have been obtained from PredictDB GTEx v8 mash-r based expression models modelling transcriptomic data by eQTL weights aligned to hg38 – at hand.
 Linkage disequilibrium (LD) reference matrices from European-ancestry reference panels aligned to hg38 genome version – at hand.
 