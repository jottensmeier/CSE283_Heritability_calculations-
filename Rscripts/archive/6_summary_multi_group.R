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

perpare_gene_annotation <- function(ens_db, finemap_gene_res) {
  finemap_gene_res <- subset(finemap_res, group != "SNP")
  gene_ids <- unique(finemap_gene_res$molecular_id)
  gene_annot <- ctwas::get_gene_annot_from_ens_db(ens_db, gene_ids)
  return(gene_annot)
}

################################################################################
################################################################################

add_gene_annotation <- function(finemap_res, snp_map, gene_annotation) {

  colnames(gene_annotation)[colnames(gene_annotation) == "gene_id"] <- "molecular_id"

  finemap_res <- anno_finemap_res(finemap_res,
                                  snp_map = snp_map,
                                  mapping_table = gene_annotation,
                                  add_gene_annot = TRUE,
                                  map_by = "molecular_id",
                                  drop_unmapped = FALSE,
                                  add_position = TRUE,
                                  use_gene_pos = "mid")

  return(finemap_res)
}

process_fine_map <- function(finemap_res_mt, snp_map, gene_annotation) {

  finemap_res_mt$pval <- z2p(finemap_res_mt$z)
  finemap_res_mt_ga <- add_gene_annotation(finemap_res_mt,
                                           snp_map,
                                           gene_annotation)

}

################################################################################
################################################################################

extract_heritability <- function(ctwas_parameters, group_label) {
  data <- data.table::data.table(
    category = names(ctwas_parameters$prop_heritability),
    percentage = ctwas_parameters$prop_heritability,
    group_label = group_label
  )
  data$percentage_label <- paste0(round(data$percentage * 100), "%")
  return(data)
}

herit_table <- function(param, gwas_n) {

  ctwas_parameters <- ctwas::summarize_param(param,
                                             gwas_n,
                                             enrichment_test = "fisher")
  print(names(ctwas_parameters))
  print(ctwas_parameters)
  data <- extract_heritability(ctwas_parameters, "multi_tissue")
  print(data)
  return(data)
}

# Generate Convergence plots
################################################################################
################################################################################

generate_frame <- function(name, param_sub_group_obj) {

  df <- data.table::data.table(name = names(param_sub_group_obj[[name]]),
                   score = unname(param_sub_group_obj[[name]]),
                   category = name)
  return(df)
}

summarise_parameters_single <- function(param, gwas_n) {

  ctwas_parameters <- ctwas::summarize_param(param,
                                             gwas_n,
                                             enrichment_test = "fisher")

  dfs <- lapply(names(ctwas_parameters),
                FUN = generate_frame,
                param_sub_group_obj = ctwas_parameters)

  # df <- dplyr::bind_rows(dfs)
  df <- rbindlist(dfs, use.names = TRUE, fill = TRUE)

  return(df)
}

################################################################################
################################################################################

extract_geneids <- function(ctwas_obj) {
  finemap_gene_res <- subset(ctwas_obj[["finemap_res"]], group != "SNP")
  return(finemap_gene_res$molecular_id)
}

extract_geneids_itter <- function(ctwas_obj) {
  all_gene_ids <- lapply(ctwas_obj, extract_geneids)
  return(unique(
         unlist(
                all_gene_ids)
                            )
                            )
}

################################################################################
################################################################################

summarise_ctwas_objects <- function(ctwas_obj_m,
                                    snp_map,
                                    gwas_n,
                                    gtf,
                                    label = "multitissue") {
  ##############################################################################
  ##############################################################################
  # What data went into the model?
  # Selected regions
  # Estimated number of causal signals (L)
  # Non-snp pips at region level
  # Primary use - 
  # Understanding why regions were included or excluded
  # Assessing model complexity (L per region)
  # QC for region filtering
  # Use for pipeline diagnostics, sensitivity analysis
  # and methodological validation.
  # region_data_mt <- ctwas_res$region_data
  ##############################################################################
  ##############################################################################
  ctwas_obj_m[["heritability_table"]] <- herit_table(
    ctwas_obj_m[["param"]],
    gwas_n
  )
  print(class(ctwas_obj_m[["heritability_table"]]))
  ctwas_obj_m[["heritability_table"]][, `group_level_label` := label]
  ##############################################################################
  ##############################################################################
  # Estimated global model parameters = how does genetics act?
  # Prior inclusion probabilities , effect size variance, enrichment,
  # PVE/ heritability
  # Primary use - architecture of genetic regulation
  # Comparing tissues / molecular types
  # Reporting mediation vs direct SNP effects
  ctwas_obj_m[["param"]] <- summarise_parameters_single(
    ctwas_obj_m[["param"]],
    gwas_n
  )
  ctwas_obj_m[["param"]][, `group_level_label` := label]
  ##############################################################################
  ##############################################################################
  # Z-score for molecular traits, (genes, spliicng etc...)
  # after integrating GWAS and QTL weights.
  # Can be used as a quick association at gene/trait level.
  # Diagnostic checks (direction, magnitude, consistency across tissues).
  # When to use -
  # Early exploration, comparing marginal signal vs fine-mapped signal,
  # sanity checks against TWAS-like results.
  # z_gene_mt <- ctwas_res[["z_gene"]]
  ##############################################################################
  ##############################################################################
  # Main causal inference output
  # Row-wise results for SNPs and molecular traits.
  # Z-score, PIPs, credible membership, posterior effect sizes.
  # Primary use -
  # Final result table, thresholding PIP > 0.8
  # Locus plots, Manuscript ready results
  gtf_ss <- prepare_gene_annotation(ctwas_obj_m[["finemap_res"]],
    gtf
  )
  ctwas_obj_m[["finemap_res"]] <- process_fine_map(
  ctwas_obj_m[["finemap_res"]],
  snp_map,
  gtf_ss
  )
  # print(class(ctwas_obj_m[["finemap_res"]]))
  # ctwas_obj_m[["finemap_res"]][, `group_level_label` := label]
  ctwas_obj_m[["finemap_res"]]$group_level_label <- label
  ##############################################################################
  ##############################################################################
  # Single effect posterior probabilities (a) from SuSie for each credible set.
  # Primary use - gene-level PIP aggregation
  # Multi-tissue / multi-trait inference.
  # susie_alpha_res_mt <- ctwas_res$susie_alpha_res
  ##############################################################################
  ##############################################################################
  # Result of region screening -
  # selected regions, estimated number of causal signals (L)
  # non-snp pips at region level
  # Prim use - understanding why regions were included or excluded
  # Assessing model complexity (L per region)
  # QC for region filtering
  # screen_res_mt <- ctwas_res$screen_res
  return(ctwas_obj_m)
}

#########################################
#    lapply inner function  decorator   #
#########################################
lapply_function_inner <- function(name,
                                  object,
                                  snp_map,
                                  gwas_n,
                                  gtf) {
print(name)
ctwas_obj <- object[[name]]
object[[name]] <- summarise_ctwas_objects(ctwas_obj,
                                          snp_map,
                                          gwas_n,
                                          gtf,
                                          label = name)

return(object[[name]])
}

summarise_ctwas_objects_iterator <- function(ctwas_obj,
                                             snp_map,
                                             gwas_n,
                                             gtf,
                                             place_holder = NULL) {

  #########################################
  # Return read and modified cTWAS object #
  #########################################
  # lapply(names(ctwas_obj),
  #        FUN = lapply_function_inner,
  #        object = ctwas_obj,
  #        snp_map = snp_map,
  #        gwas_n = gwas_n,
  #        gtf = gtf)

  result <- lapply(names(ctwas_obj),
                   FUN = lapply_function_inner,
                   object = ctwas_obj,
                   snp_map = snp_map,
                   gwas_n = gwas_n,
                   gtf = gtf)

  names(result) <- names(ctwas_obj)  # preserve names
  #########################################
  #########################################
  return(result)
}

################################################################################
################################################################################

#########################################
######     Function Dispensors     ######
#########################################

read_ctwas_disp <- list(
  `sngl_reader` = readRDS,
  `mlti_reader` = function(paths) {
    lapply(paths, readRDS)
  }
)

extract_geneid_disp <- list(`sngl_reader` = extract_geneids,
                            `mlti_reader` = extract_geneids_itter)

process_ctwas_disp <- list(`sngl_reader` = summarise_ctwas_objects,
                           `mlti_reader` = summarise_ctwas_objects_iterator)


################################################################################
################################################################################

main <- function(ctwas_obj_paths,
                 snp_map,
                 gwas_n,
                 gene_database,
                 label) {

  snp_map <- readRDS(snp_map_path)

  gwas_n <- 343621
  ##############################################################################
  ##############################################################################
  read_flag <- if (length(ctwas_obj_paths) > 1) "mlti_reader" else "sngl_reader"
  ##############################################################################
  ##############################################################################
  ctwas_obj <- read_ctwas_disp[[read_flag]](ctwas_obj_paths)

  gene_ids <- extract_geneid_disp[[read_flag]](ctwas_obj)

  gene_annot <- ctwas::get_gene_annot_from_ens_db(EnsDb.Hsapiens.v86, gene_ids)

  ctwas_obj <- process_ctwas_disp[[read_flag]](ctwas_obj,
                                               snp_map,
                                               gwas_n,
                                               gene_annot,
                                               NULL)

  extract_finemap_res <- rbindlist(lapply(ctwas_obj, function(obj) {obj[["finemap_res"]]}))

  herit_ext <- rbindlist(lapply(ctwas_obj, function(obj) {obj[["heritability_table"]]}))

  extract_parameters
}

################################################################################
################################################################################


if (sys.nframe() == 0) {

  script_dir <- get_script_dir()
  source(file.path(script_dir, "utils.R"))

  suppressPackageStartupMessages({
    library("DT")
    library("glue")
    library("dplyr")
    library("ctwas")
    library("enrichR")
    library("ggplot2")
    library("argparse")
    library("htmltools")
    library("data.table")
    library("EnsDb.Hsapiens.v86")
  })

  ######################################################

  parser <- ArgumentParser(description = "Prepare reference files for cTWAS")

  parser$add_argument("--inputdir",
                      help = "",
                      required = FALSE,
                      default = 6)

  parser$add_argument("--ncore",
                      help = "",
                      required = FALSE,
                      default = 6)

  parser$add_argument("--genome_version",
                      help = "genome_version",
                      required = TRUE)

  parser$add_argument("--outpath",
                      help = "outpath to .RDS files",
                      required = TRUE)

  ######################################################

  args <- parser$parse_args()

  output <- main()

  dir.create(args$outpath,
             recursive = TRUE, showWarnings = FALSE)

  saveRDS(output, glue('{args$outpath}{args$fname}.RDS'))
  ######################################################
}





#   ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res_param_allchrs.RDS", package = "ctwas"))
#   finemap_res_tmp <- ctwas_res[["finemap_res"]]
# finemap_res_tmp_liver <- finemap_res_tmp[
#   finemap_res_tmp$context %in% c("liver", "Liver"),
# ]
    
#   setorder(finemap_res_tmp, -"susie_pip")

#   z_gene_mt <- ctwas_res$z_gene
#   param_mt <- ctwas_res$param
#   finemap_res_mt <- ctwas_res$finemap_res
#   susie_alpha_res_mt <- ctwas_res$susie_alpha_res
#   region_data_mt <- ctwas_res$region_data
#   screen_res_mt <- ctwas_res$screen_res
#   ctwas_res_sum <- summarise_ctwas_objects(ctwas_res, snp_map, gwas_n, gene_annot)
#   finemap_res <- ctwas_res_sum[["finemap_res"]]

#   herit <- ctwas_res_sum[["heritability_table"]]
#   herit_cc <- rbindlist(list(herit_ext, herit))
#   write.csv(herit_cc, "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/summary/herit_table.csv")


#   finemap_res <- ctwas_res_sum[["finemap_res"]]
#   finemap_res_filt <- finemap_res[finemap_res$"susie_pip" > 0.8,]
#   ################################################################
#   ################################################################
#   extract_finemap_res_ss <- extract_finemap_res[chrom == 16,]


#   cc <- rbindlist(list(finemap_res, extract_finemap_res_ss))
#   cc_liver <- cc[context %in% c("liver", "Liver")]

#   cc_liver_pivot <- dcast(cc_liver,
#                           gene_name ~ group,
#                           value.var = "susie_pip",
#                           fun.aggregate = mean)

#   write.csv(cc_liver_pivot, "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/summary/susie_pip_comparison_liver.csv")



#   cc_adipose <- cc[context %in% c("adipose", "Adipose")]

#   cc_adipose_pivot <- dcast(cc_liver,
#                           gene_name ~ group,
#                           value.var = "susie_pip",
#                           fun.aggregate = mean)

#   write.csv(cc_liver_pivot, "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/summary/susie_pip_comparison_adipose.csv")
