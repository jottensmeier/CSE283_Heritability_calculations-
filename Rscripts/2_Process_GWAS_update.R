################################################################################
################################################################################
check_col <- function(gwas) {
  id_col <- dplyr::case_when(
    "ID" %in% names(gwas) ~ "ID",
    "rsid" %in% names(gwas) ~ "rsid",
    TRUE ~ stop("No variant ID column found")
  )
  return(id_col)
}

read_VCF <- function(path) {
  gwas <- VariantAnnotation::readVcf(path)
  gwas <- as.data.frame(gwasvcf::vcf_to_tibble(gwas))

  id_col <- check_col(gwas)

  z_snp <- ctwas::read_gwas(gwas,
                            id = id_col,
                            A1 = "ALT",
                            A2 = "REF",
                            beta = "ES",
                            se = "SE")

  return(z_snp)
}

################################################################################
################################################################################

# Z scores were reconstructed from GWAS summary statistics using reported effect
# sizes and two-sided p-values under a standard normal approximation.
# Specifically, Z was computed as sign(BETA) × Φ⁻¹(1 − P/2).
# This approach is equivalent to the Wald statistic (BETA/SE) in the limit
# of infinite precision.
# Minor numerical differences relative to SE-based Z scores are expected due
# to rounding of p-values and finite-precision arithmetic and do not affect
# downstream analyses.

check_GWAS_cols <- function(df) {
  required_cols <- c("BETA", "P")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Cannot compute Z score due to missing required columns: ",
      paste(missing_cols, collapse = ", "))
  }
}

calculate_z_score <- function(df) {
  # ------------------------------------------------------------------
  # Z-score reconstruction from GWAS summary statistics
  #
  # Some processed GWAS summary statistics do not include per-variant
  # standard errors (SE) or test statistics, but do report effect sizes
  # (BETA) and p-values (P). CTWAS requires per-variant Z scores rather
  # than raw effect sizes.
  #
  # Under standard GWAS practice (Wald test with two-sided p-values),
  # Z scores can be reconstructed as:
  #
  #   Z = sign(BETA) * qnorm(1 - P / 2)
  #
  # Assumptions:
  #   - P is a two-sided p-value
  #   - BETA corresponds to the effect allele (EA)
  #   - Test statistic is approximately normal under the null
  #
  # Minor numerical differences relative to SE-based Z scores are
  # expected due to p-value rounding and finite precision.
  # ------------------------------------------------------------------
  if (!"Z" %in% colnames(df)) {
    # stop("Z-score not present - file not GWAS summary statistic")
    # MODIFY DF IN PLACE
    check_GWAS_cols(df)
    # Replace zero p-values to avoid infinite z scores
    df[!is.na(P) & P == 0, P := 1e-300]
    # Compute z scores
    df[, Z := sign(BETA) * qnorm(1-P/2)]
    # Sanity checks (optional but recommended)
    #   cor(Z, BETA) should be positive
    #   Z should be approximately N(0,1) excluding tails
  }
}

RSID_selector <- function(df) {
  # Trim spaces and convert empty strings to NA
  cols_to_clean <- c("RSID", "Original_RSID", "SNPID")
  # .SD = Subset of Data — all columns listed in .SDcols.
  # In this case, .SD will contain RSID, Original_RSID, and SNPID.
  # := is data.table syntax for by-reference assignment.
  # The parentheses around cols_to_clean indicate dynamic column names,
  # i.e., assign to the columns listed in that vector.
  # lapply(.SD, ...) applies the function to each column in
  df[, (cols_to_clean) := lapply(.SD, function(x) {
  # trims(x) -> removes leading and trailing whitespace (space, tabs)
  # from each string in the column.
  # fifelse(x == "", NA_character_, x) -> repleaces empty strings ("")
  # with NA (missing values).
  # fifelse is a fast, vectorised version of ifelse from data.table
  # NA_character_ ensures the column remains character type.
    x <- trimws(x)
    fifelse(x == "", NA_character_, x)
  # dplyr::coalesce(...) is a function that returns the first
  # non-missing value across its arguments.
  }), .SDcols = cols_to_clean]
  # Pick first non-missing ID
  df[, RSID_Sel := dplyr::coalesce(RSID, Original_RSID, SNPID)]
  return(df)
}

read_TSV <- function(path) {

  gwas <- fread(path, sep = "\t")
  # MODIFY DF IN PLACE
  calculate_z_score(gwas)
  # MODIFY DF IN PLACE
  RSID_selector(gwas)

  # cols = c("RSID", "Original_RSID", "NEA", "EA", "BETA", "SE")
  cols <- c("RSID_Sel", "NEA", "EA", "BETA", "Z")
  gwas <- gwas[, ..cols]

  setnames(
  gwas,
  old = cols,
  new = c("id", "A1", "A2", "ES", "z")
  )
  message("DF read successfully")
  # setDF(gwas)
  return(gwas)
}

################################################################################
################################################################################

# convert_to_ukb_varIDs <-function(varIDs, ref_format = "%s:%s_%s_%s") {
#     varID_list <- strsplit(varIDs, split = "_|:")
#     new_varIDs <- unlist(lapply(varID_list, function(x) {
#         if (length(x) >= 4) {
#             print(x)
#             sprintf(ref_format, x[1], x[2], x[3], x[4])
#         }
#         else {
#             x
#         }
#     }))
#     return(new_varIDs)
# }

prepare_gwas_znps <- function(GWAS_Path, SNP_map, key) {

  reader <- list(`VCF` = read_VCF, `TSV` = read_TSV)

  z_snp <- reader[[key]](GWAS_Path)
  message("Total SNPs in original file")
  message(capture.output(dim(z_snp)))

  # Filter out missing, multiallelic, and strand ambiguous variants!!!
  # Convert variant IDs to match the refernce
  # Removes GWAS variants not present in LD reference!!!
  # Harmonises effect alleles with reference
  # Returns a read-to-use GWAS z-score table for cTWAS
  z_snp <- ctwas::preprocess_z_snp(z_snp, SNP_map,
                                   drop_multiallelic = TRUE,
                                   drop_strand_ambig = TRUE,
                                   varID_converter_fun = convert_to_ukb_varIDs)

  message("Total SNPs in harmonised file")
  message(capture.output(dim(z_snp)))
  invisible(z_snp)
}

################################################################################
################################################################################

extract_name <- function(GWAS_Path) {

  parent_dir <- basename(dirname(GWAS_Path))
  print("parent_dir")
  print(parent_dir)
  study_name <- gsub("\\.(db|csv|tsv|vcf|gz)", "", basename(GWAS_Path))
  print("study_name")
  print(study_name)
  label <- paste(parent_dir, study_name, sep = "_")
  print("label")
  print(label)

  return(list(`parent_dir` = parent_dir,
              `study_name` = study_name,
              `label` = label))

}

################################################################################
################################################################################

save_file <- function(data, outpath, label, disease) {

  # dir.create(glue("{outpath}/harmonised_gwas_zsnps/{dirname}/"),
  #            recursive = TRUE, showWarnings = FALSE)

  out <- glue("{outpath}/harmonised_gwas_zsnps/{disease}")
  print(dir.exists(out))
  dir.create(out,
             recursive = TRUE, showWarnings = FALSE)

  saveRDS(data, glue("{out}/{label}.RDS"))

}

################################################################################
################################################################################
################################################################################
################################################################################

main <- function(GWAS_Path, SNP_map_path, genome_version, outpath, ftype="TSV", disease) {

  print("Executing Script")
  z_snp <- prepare_gwas_znps(GWAS_Path, readRDS(SNP_map_path), ftype)

  print("Extracting name")
  study_name <- extract_name(GWAS_Path)[["study_name"]]

  print("Saving file")
  save_file(z_snp, outpath, study_name, disease)

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

  parser$add_argument("--SNP_map_path", help = "SNP map for reference", required = TRUE)

  parser$add_argument("--genome_version", help = "genome_version",
                      required = TRUE)

  parser$add_argument("--chromosome", type = "integer",
                      help = "chromosme, for subsetting", required = FALSE,
                      default = NULL)

  parser$add_argument("--outpath", help = "outpath to .RDS files",
                      required = TRUE)

  parser$add_argument("--ftype", help = "file type to read",
                      required = FALSE,
                      default = "TSV")

  parser$add_argument("--disease", help = "file type to read",
                      required = FALSE,
                      default = NULL)

  ######################################################

  args <- parser$parse_args()

  ######################################################
  main(GWAS_Path = args$GWAS_Path,
       SNP_map_path = args$SNP_map,
       genome_version = args$genome_version,
       outpath = args$outpath,
       ftype = args$ftype,
       disease = args$disease)
}