################################################################################
################################################################################

extract_GWAS_z_snp <- function(path) {
  gwas <- VariantAnnotation::readVcf(path)
  gwas <- as.data.frame(gwasvcf::vcf_to_tibble(gwas))
  z_snp <- ctwas::read_gwas(gwas,
                            id = 'rsid',
                            A1 = 'ALT',
                            A2 = 'REF',
                            beta = 'ES',
                            se = 'SE')
  return(z_snp)
}

################################################################################
################################################################################

prepare_gwas_znps <- function(GWAS_Path, snp_map) {

  z_snp <- extract_GWAS_z_snp(GWAS_Path)
  # Filter out missing, multiallelic, and strand ambiguous variants!!!
  # Convert variant IDs to match the refernce
  # Removes GWAS variants not present in LD reference!!!
  # Harmonises effect alleles with reference
  # Returns a read-to-use GWAS z-score table for cTWAS
  z_snp <- ctwas::preprocess_z_snp(z_snp, snp_map,
                                   drop_multiallelic = TRUE,
                                   drop_strand_ambig = TRUE,
                                   varID_converter_fun = ctwas::convert_to_ukb_varIDs)

  invisible(z_snp)

}

################################################################################
################################################################################

extract_name <- function(GWAS_Path) {

  parent_dir <- basename(dirname(GWAS_Path))
  study_name <- gsub("\\.(db|csv|tsv|vcf|gz)", "", basename(GWAS_Path))
  label <- paste(parent_dir, study_name, sep = "_")

  return(list(`parent_dir` = parent_dir,
              `study_name` = study_name,
              `label` = label))

}

################################################################################
################################################################################

save_file <- function(data, outpath, label) {

  # dir.create(glue("{outpath}/harmonised_gwas_zsnps/{dirname}/"),
  #            recursive = TRUE, showWarnings = FALSE)

  dir.create(glue("{outpath}/harmonised_gwas_zsnps/"),
             recursive = TRUE, showWarnings = FALSE)

  saveRDS(data,
          glue("{outpath}/harmonised_gwas_zsnps/{label}.rds"))

}

################################################################################
################################################################################
################################################################################
################################################################################

main <- function(GWAS_Path, SNP_map, genome_version, outpath) {

  z_snp <- prepare_gwas_znps(GWAS_Path, readRDS(SNP_map))

  label <- paste(extract_name(GWAS_Path)[['parent_dir']],
                 extract_name(GWAS_Path)[['study_name']],
                 sep = "/")

  save_file(z_snp, outpath, label)

}

if (sys.nframe() == 0) {
  suppressPackageStartupMessages({
    library("glue")
    library("ctwas")
    library("argparse")
    library("data.table")
  })

  ######################################################

  parser <- ArgumentParser(description = "Prepare reference files for cTWAS")

  parser$add_argument("--GWAS_Path", help = "region path", required = TRUE)

  parser$add_argument("--SNP_map", help = "SNP map for reference", required = TRUE)

  parser$add_argument("--genome_version", help = "genome_version",
                      required = FALSE,
                      default="b38")

  parser$add_argument("--chromosome", type = "integer",
                      help = "chromosme, for subsetting", required = FALSE,
                      default = NULL)

  parser$add_argument("--outpath", help = "outpath to .RDS files",
                      required = TRUE)

  ######################################################

  args <- parser$parse_args()

  ######################################################
  main(GWAS_Path = args$GWAS_Path,
       SNP_map = args$SNP_map,
       genome_version = args$genome_version,
       outpath = args$outpath)

}
