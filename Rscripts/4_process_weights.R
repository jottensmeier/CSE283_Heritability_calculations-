################################################################################
################################################################################

#########################
######## UTILITY ########
#########################

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

#########################
######## UTILITY ########
#########################

#########################
##### CONFIG PARSER #####
#########################

read_config <- function(path, key = NULL) {
  config <- yaml::read_yaml(path)
  if (!is.null(key) && key %in% names(config)) {
    config <- config[[key]]
  }
  return(config)
}

load_weight_spec <- function(weight_file,
                             model = NULL,
                             type = NULL,
                             context = NULL,
                             key = "prediction_models") {

  is_yaml <- is.character(weight_file) &&
    tools::file_ext(weight_file) %in% c("yaml", "yml")

  is_config_mode <- is_yaml && is.null(type) && is.null(context)
  is_inline_mode <- !is.null(type) && !is.null(context)

  if (is_config_mode == is_inline_mode) {
    stop("Specify either YAML config OR tyep/context (exclusively)")
  }

  if (is_config_mode) {
    print('Running Config Mode for weight harmonisation')
    config <-  read_config(weight_file, key)

  } else {
    print('Running Inline Mode for weight harmonisation')
    config <- list(list(
    model = model,
    type = type,
    context = context,
    path = weight_file
  ))

  return(config)
  }
}

#########################
######## WRITER  ########
#########################

save_weights_RDS <- function(data, outpath, label) {

  dir.create(glue("{outpath}/"),
             recursive = TRUE, showWarnings = FALSE)

  out <- glue("{outpath}/{label}.rds")

  saveRDS(data, out)
}

################################################################################
################################################################################

#########################
######## WRAPPER ########
#########################

process_weights_wrapper <- function(weight_file,
                                    region_info,
                                    snp_map,
                                    gwas_snp_ids,
                                    type,
                                    context,
                                    LD_map = NULL,
                                    load_predictdb_LD = TRUE) {

  weights <- ctwas::preprocess_weights(
    weight_file = weight_file,
    region_info = region_info,
    gwas_snp_ids = gwas_snp_ids,
    snp_map = snp_map,
    type = type,
    context = context,
    weight_name = paste0(type, "_", context),
    weight_format = "PredictDB",
    drop_strand_ambig = TRUE,
    filter_protein_coding_genes = FALSE,
    scale_predictdb_weights = FALSE,
    ##################################
    ## !!!Check these parameters!!! ##
    ##################################
    LD_map = LD_map,
    load_predictdb_LD = load_predictdb_LD,
    include_weight_LD = TRUE,
    ncore = 8,

    top_n_snps = 500
  )
  invisible(weights)
}

################################################################################
################################################################################

#########################
####### Processing ######
#########################

processing <- function(weight_file,
                       region_info,
                       snp_map,
                       gwas_snp_ids,
                       prediction_model,
                       type,
                       context,
                       outpath,
                       LD_map_path = NULL) {

  # Determine LD source
  LD_map <- if (!is.null(LD_map_path)) readRDS(LD_map_path) else NULL

  print(weight_file)
  print(glue('Null status of LD map path {is.null(LD_map_path)}.'))
  print(glue('Null status of LD map {is.null(LD_map)}.'))

  ## Calculate weights
  weight <- process_weights_wrapper(
    weight_file = weight_file,
    region_info = region_info,
    snp_map = snp_map,
    gwas_snp_ids = gwas_snp_ids,
    type = paste(prediction_model, type, sep = "_"),
    context = context,
    LD_map = LD_map,
    load_predictdb_LD = is.null(LD_map)
  )

  print("Processing finished")
  ## Save weights
  label <- paste(prediction_model, type, context, sep = "_")
  save_weights_RDS(weight, outpath, label)
  invisible(weight)
}

################################################################################
################################################################################

#########################
######## UTILITY ########
#########################

make_worker <- function(region_info,
                        snp_map,
                        gwas_snp_ids,
                        outpath,
                        LD_map_path = NULL) {

  LD_map <- if (!is.null(LD_map_path)) readRDS(LD_map_path) else NULL

  function(weight_file, context, type, prediction_model = NULL) {
    processing(
      weight_file = weight_file,
      region_info = region_info,
      snp_map = snp_map,
      gwas_snp_ids = gwas_snp_ids,
      prediction_model = prediction_model,
      type = type,
      context = context,
      outpath = outpath,
      LD_map_path = LD_map_path
    )
  }
}

#########################
#######  PARALLEL #######
#########################

single <- function(worker, weights, ncores) {
  print(weights)
  worker(
    weight_file = weights[[1]]$path,
    context = weights[[1]]$context,
    type = weights[[1]]$type,
    prediction_model = weights[[1]]$model
  )
}

parallel <- function(worker, weights, ncores) {
  # Multiple weight files.
  print("##################################")
  print(weights)
  print("##################################")
  multigroup_weights <- parallel::mclapply(
    weights, function(args) {
      print("##################################")
      print(args)
      print("##################################")
      worker(weight_file = args$path,
             context = args$context,
             type = args$type,
             prediction_model = args$model)
             },
    mc.cores = ncores
  )
  return(multigroup_weights)
}

dispatcher <- list(`single` = single,
                   `multiprocessing` = parallel)

################################################################################
################################################################################

main <- function(weight_files,
                 region_path,
                 chromosome,
                 snp_map_path,
                 gwas_snpid_path,
                 LD_map_path = NULL,
                 outpath,
                 ncores = 1)  {

  # Output directory.
  study_spec_outpath <- generate_outpath(
    outpath,
    file.path('harmonised_weights',
              extract_base_dir_label(gwas_snpid_path))
  )

  # Load GWAS SNP IDs.
  gwas_snp_ids <- readRDS(gwas_snpid_path)$id

  # Prepare regions.
  region_info <- prepare_region_info(region_path,
                                     chromosome)

  # Load SNP map.
  snp_map <- readRDS(snp_map_path)

  print(glue("Harmonising {length(weight_files)} weight files."))

  # Create worker (closure)
  worker <- make_worker(region_info,
                        snp_map,
                        gwas_snp_ids,
                        study_spec_outpath,
                        LD_map_path)

  print("Worker Loaded - beginning weight harmonisation")
  # Run in parallel
  proc_type <- if (length(weight_files) == 1) "single" else "multiprocessing"
  weights <- dispatcher[[proc_type]](worker, weight_files, ncores)

}

if (sys.nframe() == 0) {

  script_dir <- get_script_dir()
  source(file.path(script_dir, "utils.R"))

  suppressPackageStartupMessages({
    library("glue")
    library("ctwas")
    library("argparse")
    library("parallel")
    library("data.table")
  })

  ######################################################

  parser <- ArgumentParser(description = "Prepare reference files for cTWAS")

  parser$add_argument("--weight_files",
                      help = "weight file path per condition",
                      required = TRUE)

  parser$add_argument("--region_path",
                      help = "LD dir",
                      required = TRUE)

  parser$add_argument("--chromosome",
                      type = "integer",
                      help = "chromosme, for subsetting",
                      required = FALSE, default = NULL)

  parser$add_argument("--snp_map_path",
                      help = "snp_map",
                      required = TRUE)

  parser$add_argument("--gwas_snpid_path",
                      help = "GWAS SNPID references snp ids",
                      required = TRUE)

  parser$add_argument("--prediction_model",
                      help = "tissue type",
                      default = NULL,
                      required = FALSE)

  parser$add_argument("--type",
                      help = "qtl type",
                      required = FALSE,
                      default = NULL)

  parser$add_argument("--context",
                      help = "tissue type",
                      required = FALSE,
                      default = NULL)

  parser$add_argument("--LD_map_path",
                      help = paste0("A data frame with filenames of",
                                    "LD matrices and SNP information",
                                    "for the regions. Required when",
                                    "load_predictdb_LD = FALSE"),
                      required = FALSE,
                      default = NULL)

  parser$add_argument("--outpath",
                      help = "outpath to .RDS files",
                      required = TRUE)

  parser$add_argument("--ncores", help = "number of processing cores",
                      default = 1,
                      type = "integer",
                      required = FALSE)

  ######################################################
  args <- parser$parse_args()
  ######################################################
  print(args$weight_files)

  # Parse running mode - config or inline.
  weight_files <- load_weight_spec(
    weight_file = args$weight_files,
    model = args$prediction_model,
    type = args$type,
    context = args$context,
    key = "prediction_models"
  )

  main(
    weight_file = weight_files,
    region_path = args$region_path,
    chromosome = args$chromosome,
    snp_map_path = args$snp_map_path,
    gwas_snpid_path = args$gwas_snpid_path,
    LD_map_path = args$LD_map_path,
    outpath = args$outpath,
    ncores = args$ncores
  )
  ######################################################
}

    weight_files <- load_weight_spec(
    weight_file = weight_files,
    model = prediction_model,
    type = type,
    context = context,
    key = "prediction_models"
  )

#  main(
#     weight_file = weight_files,
#     region_path = region_path,
#     chromosome = NULL,
#     snp_map_path = snp_map_path,
#     gwas_snpid_path = gwas_snpid_path,
#     LD_map_path = LD_map_path,
#     outpath = outpath,
#     ncores = ncores
#   )

#   # for testing
#   zsnpid <- read.csv("/home/jottensmeier/BioAdHoc/Projects/DICE_LUNG/c_TWAS/results/harmonised_gwas_zsnps/LDL_Top_GTEX_eqtl_ukb-d-30780_irnt.vcf.gz.rds")$id
#   region_path <- "/mnt/biohome/jottensmeier/miniforge3/envs/cTWAS/lib/R/library/ctwas/extdata/ldetect/EUR.b38.ldetect.regions.RDS"
#   chromosome <- 16
#   region_info <- prepare_region_info(region_path, chromosome)
#   snp_map <- "/home/jottensmeier/BioAdHoc/Projects/DICE_LUNG/c_TWAS/results/snp_map/ukb_b38_0.1_chr16.R_snp.rds"
#   snp_map <- readRDS(snp_map)

#   weight_file <- system.file("extdata/sample_data", "expression_Liver.db", package = "ctwas")
# #   weight_adipose_file <- system.file("extdata/sample_data", "expression_Adipose_Subcutaneous.db", package = "ctwas")
#   context <- "Liver"
#   type <- "eqtl"

#   label = paste(context, type, sep = "_")

#   weights <- process_weights_wrapper(weight_file, region_info, zsnpid,
#                                      snp_map, context, type)


