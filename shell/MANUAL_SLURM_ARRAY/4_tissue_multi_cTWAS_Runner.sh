#!/bin/bash

#SBATCH --job-name=cTWAS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=3:00:00
#SBATCH --output=/logging/ctwas_runner_%A.out
#SBATCH --error=/logging/ctwas_runner_%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=###############

###############################################################
#                                                             #
#                      Main Program                           #
#                                                             #
###############################################################
Adipose="../../intermediary/harmonised_weights/ukb-d-30780_irnt/predixcan_mashr_eqtl_Adipose.rds"
harmonised_weights=($Adipose)
####################################################

Rscript "../../scripts/5_cTWAS_Runner.R" \
  --harmonised_gwas_score "/../../intermediary/harmonised_gwas_zsnps/ukb-d-30780_irnt.RDS" \
  --harmonised_weights ${harmonised_weights[@]} \
  --region_path "../../data/INPUT/HAPLOBLOCK_REF/ukb_b38_0.1_chrom_all.rds" \
  --snp_map "/../../intermediary/snp_map/ukb_b38_0.1_chrom_all.rds" \
  --LD_map "../../LD_map/ukb_b38_0.1_chrom_all.rds" \
  --genome_version "b38" \
  --outpath "../../data/OUT/multi/raw/ukb-d-30780_irnt" \
  --fname "Liver_Adipose_predixcan_mashr_eqtl"

########################################################
