################################################################################
################################################################################

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }

  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile)) {
    return(dirname(normalizePath(ofile)))
  }

  if (interactive()) {
    return(getwd())
  }

  stop("Cannot determine script directory")
}

################################################################################
################################################################################

read_gtex_eqtls <- function(path) {

  table <- read.csv(path, sep = "\t")

  extract_cols <- c("GENEID", "rsid", "SNPID", "ref", "alt", "BETA", "P")

  table <- table[, extract_cols]

  colnames(table) <-  c("gene",
                        "rsid",
                        "varID",
                        "ref_allele",
                        "eff_allele",
                        "weight",
                        "pval")

  return(table)
}

################################################################################
################################################################################

read_dice_eqtls <- function(path) {

  table <- read.delim(path, stringsAsFactors = FALSE)

  # TO REPLACE for DICE CELLTYPES!!
  extract_cols <- c("GENEID", "rsid", "SNPID", "ref", "alt", "BETA", "P")

  table <- table[, extract_cols]

  # TO REPLACE for DICE CELLTYPES!!
  colnames(table) <-  c("gene",
                        "rsid",
                        "varID",
                        "ref_allele",
                        "eff_allele",
                        "weight",
                        "pval")

  return(table)
}

read_selector <- list("GTEX" = read_gtex_eqtls, "DICE" = read_dice_eqtls)
################################################################################
################################################################################

prepare_gene_table <- function(gtf, gene_list) {

  df <- S4Vectors::mcols(gtf)[, c("gene_id", "gene_name", "gene_type")]
  df <- df[df$gene_id %in% unique(gene_list), ]
  df <- df[!duplicated(df), ]

  colnames(df) <- c("gene", "genename", "gene_type")

  return(data.frame(df))
}

################################################################################
################################################################################

extract_hblocks <- function(file_paths) {
  hblocks <- sapply(file_paths, function(path) {
  pieces <- as.list(strsplit(path, "\\.")[[1]])
  hblock <- pieces[length(pieces) - 1]
  return(hblock)
  })
}

Filter <- function(region, path_list) {
  tmp <- path_list[grepl(region, path_list)]
  matrix_path <- tmp[grepl(".RDS", tmp)]

  data.table(
    Chromosome = strsplit(matrix_path, "\\.|_")[[1]][7],
    region_id = region,
    LD_file =  matrix_path,
    SNP_file   = tmp[grepl(".Rvar", tmp)]
  )
}

extract_paths_regions_chrom <- function(high_level_chrom_path) {
  file_paths <- list.files(high_level_chrom_path)
  file_paths_long <- list.files(high_level_chrom_path, full.names = TRUE)
  hblocks <- extract_hblocks(file_paths)
  hblocks <- unique(unlist(unname(hblocks)))
  df <- rbindlist(lapply(hblocks, Filter, path_list = file_paths_long))
  return(df)
}

build_path_frame_chrom_all <- function(high_level_chrom_path) {
  print(high_level_chrom_path)
  file_paths <- list.files(high_level_chrom_path, full.names = TRUE)
  file_df <- rbindlist(lapply(file_paths, extract_paths_regions_chrom))
  return(file_df)
}

################################################################################
################################################################################

extract_id_from_db <- function(DB.path, SNPIDs, col = "RSID") {

  con <- DBI::dbConnect(RSQLite::SQLite(), DB.path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  if (!col %in% DBI::dbListFields(con, "Variants")) {
    stop("Error: column '", col, "' not found in Variants table")
  }

  if (length(SNPIDs) == 0) {
  stop("Error: no SNPIDs have been passed into 'extract_id_from_db' function.")
  }

  snps.use <- paste(DBI::dbQuoteString(con, SNPIDs), collapse = ",")

  query <- glue::glue(
    "SELECT * FROM Variants WHERE {col} IN ({snps.use})"
  )

  dbsnp <- DBI::dbGetQuery(con, query)

  if ("RSID" %in% colnames(dbsnp)) {
    data.table::setnames(dbsnp, "RSID", "rsnumber_156")
  }

  return(dbsnp)
}

extract_LD_info_from <- function(matrix_path, var_path, snps.keep) {
  matrix <- readRDS(matrix_path)
  var <- read.table(var_path, header = TRUE)

  # Find indices of SNPs to keep
  keep_idx <- which(var$id %in% snps.keep)
  # Subset rows and columns
  LDsub <- matrix[keep_idx, keep_idx, drop = FALSE]
  # Subset SNPIDs accordingly
  snp_subset <- var$id[keep_idx]
  # Get indices of the upper triangle (including diagonal if desired)
  ind <- which(upper.tri(LDsub, diag = TRUE), arr.ind = TRUE)
  # Create data.table
  LDlong <- data.table(
    RSID1 = snp_subset[ind[,1 ]],
    RSID2 = snp_subset[ind[,2 ]],
    LD   = LDsub[ind]
  )
  print("Done")
  print(matrix_path)
  print(dim(LDlong))
  return(LDlong)
}


high_level_chrom_path <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/data/LD_reference/LD_BioBank/"
chrom_paths <- build_path_frame_chrom_all(high_level_chrom_path)
rsids <- unique(qtl_table$rsid)

matrix_path<- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/data/LD_reference/LD_BioBank//1/ukb_b38_0.1_chr1.R_snp.100360849_101575460.RDS"
var_path <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/data/LD_reference/LD_BioBank//1/ukb_b38_0.1_chr1.R_snp.100360849_101575460.Rvar"

out <- lapply(seq_len(nrow(chrom_paths)), function(i) extract_LD_info_from(chrom_paths[i, ]$LD_file, chrom_paths[i, ]$SNP_file, rsids))

cov_table <- extract_LD_info_from(matrix_path, var_path, rsids)


order <- c("gene", "RSID1", "RSID2", "LD")
merged <- merge(cov_table, qtl_table[,c("rsid", "gene")], by.x ="RSID1", by.y ="rsid")
merged <- merged[,order]

rsid <- weights$rsid
varid <- weights$varID
snpd.DB.path <- "/home/jrocha/BioAdHoc/reference_panels/dbSNP/hg38/dbSNP156.db"


converted <- extract_id_from_db(snpd.DB.path, varid, col = 'SNPID')
weights_ss <- weights[weights$rsid %in% converted$rsnumber_156, ]

snps.keep <- weights_ss$rsid





################################################################################
################################################################################

parallele_processes <- function(qtl_paths, qtl_dbase, gtf, outdir, snpd.DB.path, ncores) {

  out <- mclapply(qtl_paths, function(qtl_path) {

  # Load QTL table
  qtl_table <- read_selector[[qtl_dbase]](qtl_path)
  outname <- tools::file_path_sans_ext(basename(qtl_path))

  # Prepare gene table using imported GTF
  gene_table <- prepare_gene_table(gtf, unique(qtl_table$gene))

  # Extract LD matrix
  # LD_infromation <- extractLD(snpd.DB.path, SNPIDs)

  # Build PredictDB
  ctwas::create_predictdb_from_QTLs(
      weight_table = qtl_table,
      gene_table   = gene_table,
      cov_table    = NULL,
      use_top_QTL  = FALSE,
      select_by    = "pval",
      outputdir    = outdir,
      outname      = outname
  )

  # Return path of created DB
  file.path(outdir, paste0(outname, ".db"))

  }, mc.cores = ncores)

  return(out)

}

################################################################################
################################################################################


################################################################################
################################################################################

main <- function(qtl_data_base, qtl_paths, gtfpath, outdir, ncores = 1) {

  # Import GTF once
  gtf <- rtracklayer::import(gtfpath, feature.type = "gene")

  qtl_paths <- path_parser(qtl_paths, qtl_data_base)
  print(qtl_paths)

  out <- parallele_processes(qtl_paths, qtl_data_base, gtf, outdir, ncores)

  return(out)
}

################################################################################
################################################################################

if (sys.nframe() == 0) {

  script_dir <- get_script_dir()
  source(file.path(script_dir, "utils.R"))

  suppressPackageStartupMessages({
    library("DBI")
    library("glue")
    library("ctwas")
    library("RSQLite")
    library("parallel")
    library("argparse")
    library("data.table")
    library("rtracklayer")
  })

  ######################################################

  parser <- ArgumentParser(description = "Prepare reference files for cTWAS")

  parser$add_argument("--qtl_data_base",
                      help = "qtl database origin - e.g. DICE or GTEX",
                      required = TRUE,
                      choices = c("DICE", "GTEX"))

  parser$add_argument("--qtl_paths",
                      nargs = "*",
                      help = "path to qtl yaml file or n number of paths",
                      required = TRUE)

  parser$add_argument("--gtf_path",
                      help = "path to GTF file.",
                      required = TRUE)

  parser$add_argument("--outpath", help = "outpath to .RDS files",
                      required = TRUE)

  parser$add_argument("--genome_version", help = "genome_version",
                      required = TRUE)

  parser$add_argument("--ncores", help = "number of processing cores",
                      default = 5,
                      type = "integer",
                      required = FALSE)

  ######################################################
  args <- parser$parse_args()
  ######################################################

  print(args)

  main(
    qtl_data_base = args$qtl_data_base,
    qtl_paths = args$qtl_paths,
    gtfpath = args$gtf_path,
    outdir = args$outpath,
    ncores = args$ncores
  )
}

# DB.path <- "/home/jrocha/BioAdHoc/reference_panels/dbSNP/hg38/dbSNP156.db"
################################################################################
################################################################################
# path <- "/home/jottensmeier/BioAdHoc/Projects/CTCF/Analysis_GRCh38_New/Intermediary/Input_Clean/filt_eQTL_catalogue/filtered/eQTL_Catalogue_1e_neg4_GTEx_ge_adipose_subcutaneous_cleaned.tsv"
# outdir <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/c_TWAS/results/eqtl_weights"
# tazle <- read_gtex_eqtls(path)
# gene_table <- prepare_gene_table(gtfp, unique(table$gene))
# gtfp <- paste0("/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_GALAXY/",
#                "reference/GRCh38-2020-A_build/gtf/",
#                "gencode.v32.primary_assembly.annotation.filtered.gtf")
# yaml::read_yaml()
################################################################################
################################################################################
# qtlPoc <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/c_TWAS/results/eqtl_weights/eQTL_Catalogue_1e_neg4_GTEx_ge_adipose_subcutaneous_cleaned.db"
# con <- dbConnect(RSQLite::SQLite(), dbname = qtlPoc)
# table <- dbReadTable(con, "weights")
# tab2 <-  dbReadTable(con, "extra")
# dbDisconnect(con)
################################################################################
################################################################################
# rename_cols <- function(table, mapper) {
#   cols <- colnames(table)
#   to_replace <- cols %in% names(mapper)
#   cols[to_replace] <- unname(mapper[cols[to_replace]])
#   colnames(table) <- cols
#   invisible(table)
# }
# ################################################################################
# ################################################################################
# read_gtex_eqtls <- function(path) {
#   table <- read.csv(path, sep = "\t")
#   mapper <- c("GENEID" = "gene",
#               "rsid" = "rsid",
#               "SNPID" = "varID",
#               "ref" = "ref_allele",
#               "alt" = "eff_allele",
#               "BETA" = "weight",
#               "P" = "pval")
#   table <- rename_cols(table, mapper)
#   table <- table[, unname(mapper), drop = FALSE]
#   invisible(table)
# }





################################################################################
################################################################################

#' @title Computes LD for weight variants using reference LD
#'
#' @param weights a list of preprocessed weights.
#'
#' @param region_info a data frame of region definitions.
#'
#' @param LD_map a data frame with filenames of LD matrices and SNP information for the regions.
#' Required when \code{load_predictdb_LD = FALSE}.
#'
#' @param snp_map a list of SNP-to-region map for the reference.
#' If NUll, it will reads SNP info from the "SNP_file" column of LD_map.
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param snpinfo_loader_fun a user defined function to load SNP information file,
#' if SNP information files are not in standard cTWAS reference format.
#'
#' @param ncore The number of cores used to parallelize computation.
#'
#' @importFrom parallel mclapply
#' @importFrom Matrix bdiag
#' @importFrom logging loginfo
#'
#' @return a list of processed weights, with LD of weight variants included.
#'
#' @export
# compute_weight_LD_from_ref <- function(weights,
#                                        region_info,
#                                        LD_map,
#                                        snp_map = NULL,
#                                        LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
#                                        LD_loader_fun = NULL,
#                                        snpinfo_loader_fun = NULL,
#                                        ncore = 1) {

#   LD_format <- match.arg(LD_format)

#   if (!inherits(weights,"list"))
#     stop("'weights' should be a list!")

#   if (!inherits(LD_map,"data.frame"))
#     stop("'LD_map' should be a data frame!")

#   weight_info <- lapply(names(weights), function(x){
#     as.data.frame(weights[[x]][c("chrom", "p0","p1", "molecular_id", "weight_name", "type","context")])})
#   weight_info <- do.call(rbind, weight_info)
#   weight_info$weight_id <- paste0(weight_info$molecular_id, "|", weight_info$weight_name)
#   # get the regions overlapping with each gene
#   for (k in 1:nrow(weight_info)) {
#     chrom <- weight_info[k, "chrom"]
#     p0 <- weight_info[k, "p0"]
#     p1 <- weight_info[k, "p1"]
#     idx <- which(region_info$chrom == chrom & region_info$start <= p1 & region_info$stop > p0)
#     weight_info[k, "region_id"] <- paste(sort(region_info[idx, "region_id"]), collapse = ",")
#   }

#   # compute LD for weight variants on each chromosome
#   chrs <- sort(unique(weight_info$chrom))
#   for (b in chrs) {
#     loginfo("Computing LD for variants in weights on chr%s", b)
#     weightinfo <- weight_info[weight_info$chrom == b, ]
#     if (nrow(weightinfo) > 0) {
#       weight_region_ids <- names(sort(-table(weightinfo$region_id)))
#       weight_LD_list <- mclapply_check(weight_region_ids, function(x){
#         # load the R_snp and SNP info for the region
#         # and extract LD for the weight variants
#         curr_region_LD_list <- list()
#         curr_region_ids <- unlist(strsplit(x, ","))
#         curr_region_idx <- match(curr_region_ids, LD_map$region_id)

#         LD_matrix_files <- unlist(strsplit(LD_map$LD_file[curr_region_idx], split = ","))
#         stopifnot(all(file.exists(LD_matrix_files)))

#         if (length(LD_matrix_files) > 1) {
#           R_snp <- lapply(LD_matrix_files, load_LD, format = LD_format, LD_loader_fun = LD_loader_fun)
#           R_snp <- suppressWarnings(as.matrix(bdiag(R_snp)))
#         } else {
#           R_snp <- load_LD(LD_matrix_files, format = LD_format, LD_loader_fun = LD_loader_fun)
#         }

#         # load SNP info of the region
#         # if snp_map is available, reads SNP info from snp_map;
#         # otherwise, reads SNP info from the "SNP_file" column of LD_map.
#         if (!is.null(snp_map)){
#           snpinfo <- as.data.frame(rbindlist(snp_map[curr_region_ids], idcol = "region_id"))
#         } else {
#           SNP_info_files <- LD_map$SNP_file[curr_region_idx]
#           stopifnot(all(file.exists(SNP_info_files)))
#           snpinfo <- read_snp_info_files(SNP_info_files, snpinfo_loader_fun = snpinfo_loader_fun)
#         }

#         rownames(R_snp) <- snpinfo$id
#         colnames(R_snp) <- snpinfo$id
#         weight_ids <- weightinfo[weightinfo$region_id == x, "weight_id"]

#         for (weight_id in weight_ids) {
#           wgt_snp_ids <- rownames(weights[[weight_id]]$wgt)
#           R_wgt <- R_snp[wgt_snp_ids, wgt_snp_ids, drop=FALSE]
#           curr_region_LD_list[[weight_id]] <- R_wgt
#         }
#         curr_region_LD_list
#       }, mc.cores = ncore, stop_if_missing = TRUE)

#       weight_LD_list <- unlist(weight_LD_list, recursive = FALSE)
#       for(weight_id in names(weight_LD_list)){
#         weights[[weight_id]][["R_wgt"]] <- weight_LD_list[[weight_id]]
#       }
#     }
#   }
#   return(weights)
# }


################################################################################
################################################################################
################################################################################
##################################### EXTRA ####################################
################################################################################
################################################################################
################################################################################

extract_hblocks <- function(path) {
  pieces <- as.list(strsplit(path, "\\.")[[1]])
  hblock <- pieces[length(pieces) - 1]
  return(hblock)
}

Filter <- function(region, path_list) {
  tmp <- path_list[grepl(region, path_list)]
  matrix_path <- tmp[grepl(".RDS", tmp)]

  data.table(
    Chromosome = strsplit(matrix_path, "\\.|_")[[1]][7],
    region_id = region,
    LD_file =  matrix_path,
    SNP_file   = tmp[grepl(".Rvar", tmp)]
  )
}

filter_region_paths <- function(region, path_list) {
  tmp <- path_list[grepl(region, path_list)]
  # Identify LD and SNP files
  ld_file <- tmp[grepl("\\.RDS$", tmp)]
  snp_file <- tmp[grepl("\\.Rvar$", tmp)]

  if (length(ld_file) != 1 || length(snp_file != 1)) {
    warning(
      sprintf(
        "Region %s: expected. 1 LD and 1 SNP file, found %d and %d"
        )
      )
  }
  chrom <- strsplit(basename(ld_file), "\\.|_")[[1]][7]

  data.table(
    Chromosome = chrom,
    region_id = region,
    LD_file = ld_file,
    SNP_file = snp_file
  )
}

extract_paths_regions_chrom <- function(low_level_chrom_path_hblock) {
  file_paths <- list.files(low_level_chrom_path_hblock, full.names = TRUE)
  hblocks <- sub(".*\\.(.*)\\..*$", "\\1", basename(file_paths))
  dt <- rbindlist(lapply(hblocks, filter_region_paths, path_list = file_paths))
  return(df)
}

