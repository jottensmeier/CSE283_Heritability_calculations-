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
# Computation of z-score from molecular traits.
# This step performs the transcriptome-wide association stidy (TWAS) on each
# molecular trait, with a prediction model. The underlying calculation is based
# on S-PrediXcan.
# Notes -
# compute_gene_z
# convert SNP_level GWAS scores into gene (molecular trait) z-cores using the
# precomputed prediction models.
################################################################################
# For each gene g:
# Z_g = (w_gᵀ z_snp) / sqrt(w_gᵀ R_g w_g)
# wg = vector of SNP weights for gene g
# Z snp = GWAS SNP z-score
# Rg = LD matrix mong SNPs used by gene g
################################################################################
# This is the PrediXcan / S-PrediXcan test statistic, but computed explicitly.
# - Must ensure that variants are already harmonised. (preprocess_z_snp())
# - SNPs must have weight + LD information. (preprocess_weights())
# - Each weight entry must contain :
# - `wgt` : SNP
# filter_z_gene_by_group_size
################################################################################

extract_gene_score_wrapper <- function(z_snp,
                                       weights,
                                       ncore,
                                       min_group_size = 100) {

  z_gene <- ctwas::compute_gene_z(z_snp = z_snp,
                                  weights = weights,
                                  ncore = 6)

  z_gene <- ctwas::filter_z_gene_by_group_size(z_gene = z_gene,
                                               min_group_size = min_group_size)
  invisible(z_gene)
}

################################################################################
################################################################################
# TWAS is run region by region (haploblock by haploblock) - so data needs to be
#  preapred by haploblock.
# Assemble inputs needed for cTWAS fine-mapping by constructing `region_data`,
# a per-region container holding:
# SNPID and SNP Z=score,
# Gene IDs and gene Z-score,
# Group, type, and context annotation,
# Region boundaries and thinning metadata.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Construct a unified SNP reference table
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Assigns genes and SNPs to regions (chormosome-wide)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Supports SNP thinning for speed
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Explic handles boundary genes to prevent double countin
# (using prediction scores)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Enforces SNP count limit per region.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Injects Z-scores into region_data
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Structure ->
# region_data[[ region_id ]]
# gid, sid - gene and snp ids.
# z_gene, z_snp - zscore table,
# types, contexts, groups - model metadata (e.g. eqtl/sqtl, tissue, group)
# thin - thinning indicator.
# minpos, maxpos : effective region bounds.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Post data processing and computation of z-scores of molecular traits
# authors assemble the input data for all the regions using the function
# assemble_region_data. It assigns molecular traits, variants and their Z-scores
#  to each region.The assignment of a molecular trait to regions is based on the
# weights, i.e. variants in the prediction model, of the molecular trait.
# The function has a few arguments to control behaviour -
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# maxSNP sets a max on the number of varaints that can be in a single region
# to prevent memory issues during fine-mapping. If the number of SNPs in a
# region is more than maxSNP,
# it will, by default, choose the top maxSNPs with the highest Z-scores.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# min_group_size sets the minimum number of variables in a group. Groups with
# fewer, variables will be removed from the analysis, because cTWAS cannot
# reliably estimate parameters for groups with too few variables. ncore
# specifies the number of cores to use when parallelising over regions.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This function returns region_data, a list object containing input data
# (IDs of molecular traits and varaints, their Z-scores, genomic boundaries,
# etc.) for each of the regionss.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################
################################################################################

assemble_input_data_for_regions_wrapper <- function(region_info,
                                                    z_snp,
                                                    z_gene,
                                                    weights,
                                                    snp_map,
                                                    thin = 1,
                                                    maxSNP = 20000,
                                                    min_group_size = 100,
                                                    ncore = 6) {

  region_data <- ctwas::assemble_region_data(region_info = region_info,
                                             z_snp = z_snp,
                                             z_gene = z_gene,
                                             weights = weights,
                                             snp_map = snp_map,
                                             thin = thin,
                                             maxSNP = 20000,
                                             min_group_size = 100,
                                             adjust_boundary_genes = TRUE,
                                             ncore = 6)

  invisible(region_data)
}

################################################################################
############################ EXPRESSION MAXIMISATION ###########################
################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Use est_param() function to estimate two sets of parameters - the prior
# inclusion probabilies and the prior effect size variance for moelcualr traits
# and variants separately.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# cTWAS can analyse multiple groups of molecular traits. This step will estimate
# the prior parameters for each group. However, sharing of prior effect variance
# parameters among groups is allowed.
# This may improve the accuracy of parameter estimation.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# In cTWAS the latent variables are which SNPs or genes are truly causal
# in each region.
# Direcoty optimising the likelihood is infeasible because causuality is unkown
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# How EM works -
# 1) E-step (Expectation)
# Compute the expected probability that each SNP or gene is causal
# given the current parameter estimates.
# 2) M-step (Maximisation)
# Update parameters (e.g. group priors and effect-size variance) to maximise
# the expected likelihood from the E-step.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# These steps repeat until the log-likelihood converges or the maximum number
# of iterations is reached.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Why EM is used in cTWAS
# - Handles uncertainty in causal variables naturally.
# - Enables stables estimation of enrichment and variance parameters.
# - Scales well when many regions are analysed jointly.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# There are several optoins to control how parameters are shared among groups:
# This can be done by specifying the group_prior_var_structure:
# "shared_all" - (dflt) allows groups, including SNPs, to share same var param
# "shared_type" - allows all group in one molec QTL type to share the same
# variance param.
# "shared_context" - allows all non-snp groups (i.e. all molecular traits)
# to share the same variance parameter.
# "shared_nonSNP" - allows all non-snp groups (i.e. all molecular traits)
# to share the same varaince parameter.
# "independent" - allows all groups to have their own separate variance params.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This step will run the EM algorithm to estimate parameters and return the
# estimated parameters (group_prior, group_prior_var, etc).
# For technical reasons (see Methods of the cTWAS paper) it will run two rounds
# of EM algorithm.
# In complete data analysis, parameter estimation should be done using data
# from the entire genome.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################
################################################################################

estimate_parameters_wrapper <- function(region_data,
                                      group_prior_var_structure = "shared_all",
                                        niter_prefit = 3,
                                        niter = 50,
                                        min_group_size = 100,
                                        ncore = 6) {

  param <- ctwas::est_param(region_data = region_data,
                            group_prior_var_structure = group_prior_var_structure,
                            niter_prefit = niter_prefit,
                            niter = niter,
                            min_group_size = min_group_size,
                            ncore = ncore)

  invisible(param)

}

################################################################################
################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Screening regions
# assemble_region_data: key Concepts and workflow
# This code block defines the data assembly layer of cTWAS. Its pruprose is to
# construct a standardised `region_data` object that integrates SNPs, genes,
# z-scores, LD references, and grouping metadata for downstream EM estimation
# and fine-mapping.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# High-level Purpose -
# If thin <1 when preparing region_data, regions must be expanded to equal
# those of all snps for screening regions and fine-mapping. If the number
# fo SNPs in a region is more than masSNP, it will choose the top maxSNP
# SNPs with the highest Z-score. This is not needed when thin=1.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Screening regions
# Used screen regions function to perform a screening process to select regions
# with likely causal signals in molecular traits.
# cTWAS first selects regions with significant GWAS signals
# (min. p-value < 5e-8) . Then, it will perform an initial fine-mapping
# without LD (L=1) and compute the non-SNP PIP i.e. total PIP of all molecular
# traits, for each region (using estimated priors). It then
# selects regions by applying a min_nonSNP_PIP cutoff (0.5 by default). Only
# these selected regions would be subject to full fine-mapping analysis in the
# final step.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# A strong non-SNP signal in cTWAS means that a gene (or molecular trait)
# shows strong association with the phenotype after aggregating many weak GWAS
# SNP signals through an eQTL model, even though no single SNP is genome-wide
# significant.
# This occurs when genetic effects are polygenic and mediated through gene
# expresion: individual SNP effects are small, but their LD-aware weighted
# combination produces a strong gene-level signal. cTWAS is desigend to detect
# these cases, which standard SNP-wise GWAS often misses.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################
################################################################################

screen_regions_wrapper <- function(region_data,
                                   snp_map,
                                   z_snp,
                                   maxSNP,
                                   ncore,
                                   group_prior,
                                   group_prior_var,
                                   min_nonSNP_PIP = 0.5,
                                   thin = 1) {

  if (thin < 1){
    region_data <- ctwas::expand_region_data(region_data = region_data,
                                             snp_map = snp_map,
                                             z_snp = z_snp,
                                             maxSNP = 20000,
                                             ncore = 6)
  }

  screen_res <- ctwas::screen_regions(region_data,
                                      group_prior = group_prior,
                                      group_prior_var = group_prior_var,
                                      min_nonSNP_PIP = min_nonSNP_PIP,
                                      ncore = 6)

  invisible(screen_res)

  # invisible( list(`screened_region_data` = screen_res$screened_region_data,
  #                 `screen_summary` = screen_res$screen_summary ) )
}

################################################################################
############################  Fine-mapping regions  ############################
################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This code implements the core cTWAS fine-mapping workflow, integrating GWAS
# summary statistics, LD structure, and eQTL-derived gene models to identify
# likely causal genes and SNPs within genomic regions.
# The pipeline operates on predefined regions that contain both SNPs and genes
# each assigned to biologically meaningful groups (e.g. SNPs, tissue-specific
# genes). These groups allow cTWAS to apply strucutred priors reflecting
# different prior beliefs about SNP- vs gene-mediated effects.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# At a high level, the workflow proceeeds as follwos -
# 1) Input validation and filtering
# - Regions are checked for sufficient numbers of SNPs and genes. Regions
# that do not meet minimum thresholds are skipped to ensure model
# identifiability.
# 2) Prior construction
# Group-specific prior inclusion probabilities and effect-size variances are
# initialised. If not supplied, defaults or uniform priors are used.
# 3) Data extraction per region.
# For each region, SNP-gene z-scores, group labels, and annotations are
# assembled into a unified vector.
# 4) LD and correlation handling
# When LD is available, SNP-SNP, gene-gene, and SNP-gene correlations are
# combined into a single correlation matrix. A no-LD mode is also supported.
# 5) Model fitting (SuSiE/SER)
# The main model fits a multi-effect Bayesian sparese regression (SuSiE) to
# jointly fine-map SNPs and genes. A simiplified single-effect regression (SER)
# model is used for screening and EM steps.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################
################################################################################

Fine_map_regions_with_LD_wrapper <- function(screened_region_data,
                                             LD_map,
                                             weights,
                                             group_prior,
                                             group_prior_var,
                                             L = 5,
                                             ncore = 4) {

  res <- ctwas::finemap_regions(screened_region_data,
                                LD_map = LD_map,
                                weights = weights,
                                group_prior = group_prior,
                                group_prior_var = group_prior_var,
                                L = 5,
                                save_cor = TRUE,
                                cor_dir = "./cor_matrix",
                                ncore = 4,
                                force_compute_cor = TRUE)

  finemap_res <- res$finemap_res
  susie_alpha_res <- res$susie_alpha_res

  invisible(res)
}

################################################################################
################################################################################

ctwas_main_wrapper <- function(z_snp,
                               weights,
                               snp_map,
                               region_info,
                               LD_map,
                               maxSNP = 20000,
                               min_group_size = 100,
                               ncore = 6,
                               min_nonSNP_PIP = 0.5) {

  z_gene <- extract_gene_score_wrapper(z_snp = z_snp,
                                       weights = weights,
                                       ncore = ncore,
                                       min_group_size = min_group_size)

  region_data <- assemble_input_data_for_regions_wrapper(region_info = region_info,
                                                         z_snp = z_snp,
                                                         z_gene = z_gene,
                                                         weights = weights,
                                                         snp_map = snp_map,
                                                         thin = 1,
                                                         maxSNP = maxSNP,
                                                min_group_size = min_group_size,
                                                         ncore = ncore)

  params <- estimate_parameters_wrapper(region_data = region_data,
                                        group_prior_var_structure = "shared_all",
                                        niter_prefit = 3,
                                        niter = 50,
                                        min_group_size = min_group_size,
                                        ncore = ncore)

  group_prior <- params$group_prior
  group_prior_var <- params$group_prior_var

  screen_res <- screen_regions_wrapper(region_data = region_data,
                                       snp_map = snp_map,
                                       z_snp = z_snp,
                                       maxSNP = maxSNP,
                                       ncore = ncore,
                                       group_prior = group_prior,
                                       group_prior_var = group_prior_var,
                                       min_nonSNP_PIP = min_nonSNP_PIP)

  screened_region_data <- screen_res$screened_region_data
  # screen_summary <- screen_res$screen_summary

  res <- Fine_map_regions_with_LD_wrapper(screened_region_data = screened_region_data,
                                          LD_map = LD_map,
                                          weights = weights,
                                          group_prior = group_prior,
                                          group_prior_var = group_prior_var,
                                          L = 5,
                                          ncore = ncore)

  invisible(list(`z_gene` = z_gene,
                 `param` = params,
                 `finemap_res` = res$finemap_res,
                 `susie_alpha_res` = res$susie_alpha_res,
                 `region_data` = region_data,
                 `screen_res` = screen_res))

}

################################################################################
################################################################################

weights_reader <- function(weights_path) {
  #############
  # Decorator #
  #############
  multi_run <- function(weights_path) {
    weights <- list()
    for (weight_file in weights_path) {
        weights <- c(weights, readRDS(weight_file))
      }
    return(weights)
  }
  #############
  # Decorator #
  #############
  single_run <- function(path) { readRDS(path) }
  #############
  selector <- c(`multitissue` = multi_run,
                `single` = single_run)
  #############
  ctwas_setting <- if (length(weights_path) == 1) "single" else "multitissue"
  statement <- paste0("Running ",
                      ctwas_setting,
                      " for analysis of follwoing file/s ",
                      weights_path)

  print(statement)
  #############
  return(selector[[ctwas_setting]](weights_path))
}

################################################################################
################################################################################

save_cTWAS_objects <- function(ctwas_obj, fname, outpath) {
  sep <- .Platform$file.sep

  decorator <- function(sd, fname, ctwas_obj) {
    out <- glue("{outpath}{sep}{sd}{sep}")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)

    out_fname <- glue("{out}{fname}.RDS")
    saveRDS(ctwas_obj[[sd]], out_fname)
    }

  dir.create(outpath,
            recursive = TRUE, showWarnings = FALSE)

  subdirs <- list("z_gene",
                  "param",
                  "finemap_res",
                  "susie_alpha_res",
                  "region_data",
                  "screen_res")

  lapply(subdirs, decorator, fname=fname, ctwas_obj=ctwas_obj)
}

################################################################################
################################################################################

main <- function( z_snp_path,
                  weights,
                  snp_map_path,
                  LD_map_path,
                  region_data_path,
                  maxSNP = 20000,
                  min_group_size = 100,
                  ncore = 6) {

  output <- ctwas_main_wrapper(z_snp = readRDS(z_snp_path),
                               weights = weights,
                               snp_map = readRDS(snp_map_path),
                               LD_map = readRDS(LD_map_path),
                              region_info = readRDS(region_data_path),
                               maxSNP = 20000,
                               min_group_size = 100,
                               ncore = 6)

  # output <- ctwas_main_wrapper(z_snp = readRDS(z_snp_path),
  #                              weights = readRDS(),
  #                              snp_map = readRDS(snp_map_path),
  #                              LD_map = readRDS(LD_map_path),
  #                             region_info = readRDS(region_data_path),
  #                              maxSNP = 20000,
  #                              min_group_size = 100,
  #                              ncore = 6)
  invisible(output)
}

################################################################################
################################################################################

if (sys.nframe() == 0) {

  script_dir <- get_script_dir()
  source(file.path(script_dir, "utils.R"))

  suppressPackageStartupMessages({
    library("glue")
    library("ctwas")
    library("argparse")
    library("data.table")
  })

  ######################################################

  parser <- ArgumentParser(description = "Prepare reference files for cTWAS")

  parser$add_argument("--harmonised_gwas_score",
                      help = "LD dir",
                      required = TRUE)

  parser$add_argument("--harmonised_weights",
                      help = "weight file path per condition",
                      required = TRUE,
                      nargs = "+")

  parser$add_argument("--region_path",
                      help = "GWAS SNPID references snp ids",
                      required = TRUE)

  parser$add_argument("--snp_map",
                      help = "snp_map",
                      required = TRUE)

  parser$add_argument("--LD_map",
                      help = "LD_map",
                      required = TRUE)

  parser$add_argument("--maxSNP",
                      help = "",
                      required = FALSE,
                      type = "integer",
                      default = 20000)

  parser$add_argument("--min_group_size",
                      help = "",
                      required = FALSE,
                      type = "integer",
                      default = 100)

  parser$add_argument("--ncore",
                      help = "",
                      required = FALSE,
                      default = 6)

  parser$add_argument("--genome_version",
                      help = "genome_version",
                      required = FALSE,
                      default="b38")

  parser$add_argument("--chromosome",
                      type = "integer",
                      help = "chromosme, for subsetting",
                      required = FALSE,
                      default = NULL)

  parser$add_argument("--outpath",
                      help = "outpath to .RDS files",
                      required = TRUE)

  parser$add_argument("--fname",
                      help = "")
  ######################################################
  args <- parser$parse_args()
  ######################################################
  output <- main(z_snp_path = args$harmonised_gwas_score,
                 weights = weights_reader(args$harmonised_weights),
                 snp_map_path = args$snp_map,
                 LD_map_path = args$LD_map,
                 region_data_path = args$region_path,
                 maxSNP = args$maxSNP,
                 min_group_size = args$min_group_size,
                 ncore = args$ncore)
  ######################################################
  save_cTWAS_objects(output, args$fname, glue('{args$outpath}'))
  ######################################################
}


# out <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/data.preprocessed/cTWAS_output/"
# files <- sapply(dir(out), function(name) {glue("{out}{name}")})
# names(files) = gsub(".RDS", "", names(files))
# ctwas_obj <- lapply(files, readRDS)

# A_enet <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/results/harmonised_weights/ukb-d-30780_irnt/predixcan_enet_eqtl_Adipose.rds"
# L_enet <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/results/harmonised_weights/ukb-d-30780_irnt/predixcan_enet_eqtl_Liver.rds"
# A_mashr <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/results/harmonised_weights/ukb-d-30780_irnt/predixcan_mashr_eqtl_Adipose.rds"
# L_mashr <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/results/harmonised_weights/ukb-d-30780_irnt/predixcan_mashr_eqtl_Liver.rds"
# A_no_model <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/results/harmonised_weights/ukb-d-30780_irnt/no_model_eqtl_Adipose.rds"

# weights <- c(readRDS(A_enet), readRDS(L_enet), readRDS(A_mashr), readRDS(L_mashr), readRDS(A_no_model))

# harmonised_gwas_score <-  "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/results/harmonised_gwas_zsnps/ukb-d-30780_irnt.RDS"
# # harmpnised_weights <- "/home/jottensmeier/BioAdHoc/Projects/DICE_LUNG/c_TWAS/results/harmonised_weights/predixcan_elastic_net_Liver_eqtl.rds"
# snp_map_path <-   "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/results/snp_map/ukb_b38_0.1_chrom_all.rds"
# LD_map_path <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/results/LD_map/ukb_b38_0.1_chrom_all.rds"
# region_data_path <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/results/region_info/ukb_b38_0.1_chrom_all.rds"
# maxSNP <-  2000
# min_group_size <- 100
# ncore <- 6
# genome_version <- "b38"
# chromosome <- 16

# main(z_snp_path = harmonised_gwas_score,
#      weights_path = harmpnised_weights,
#      snp_map_path = snp_map_path,
#      LD_map_path = LD_map_path,
#      region_data_path = region_data_path,
#      maxSNP = maxSNP,
#      min_group_size = min_group_size,
#      ncore = ncore)

################################################################################
################################################################################

# ctwas_sumstats_wrapper <- function(z_snp,
#                                    weights,
#                                    region_info,
#                                    LD_map,
#                                    snp_map) {

#   ctwas_res <- ctwas_sumstats(z_snp,
#                               weights,
#                               region_info,
#                               LD_map,
#                               snp_map,
#                               thin = 1,
#                               maxSNP = 20000,
#                               min_group_size = 100,
#                               group_prior_var_structure = "shared_all",
#                               min_nonSNP_PIP = 0.5,
#                               min_abs_corr = 0.1,
#                               ncore = 6,
#                               ncore_LD = 4,
#                               save_cor = TRUE,
#                               cor_dir = "./cor_matrix",
#                               force_compute_cor = TRUE)

#   invisible(ctwas_res)
# }

