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

paste_path <- function(name, path) { glue("{path}{name}") }

extract_path <- function(path, name) {
  base_path <- glue("{path}{name}/")
  files <- list.files(base_path)

  setNames(
    lapply(files, paste_path, path = base_path),
    files
  )
}

################################################################################
################################################################################

extract_gwas_n <- function(metatable_path, study) {
  table <- read.csv(metatable_path, sep=",")
  gwas_n <- table[table$Sample.name == study, ]$sample_size
  return(as.integer(gwas_n))
}

################################################################################
################################################################################

prepare_gene_annotation <- function(ens_db, finemap_res) {
  finemap_gene_res <- subset(finemap_res, group != "SNP")
  gene_ids <- unique(finemap_gene_res$molecular_id)
  gene_annot <- ctwas::get_gene_annot_from_ens_db(ens_db, gene_ids)
  return(gene_annot)
}

################################################################################
################################################################################

add_gene_annotation <- function(finemap_res, snp_map, gene_annotation) {

  colnames(gene_annotation)[colnames(gene_annotation) == "gene_id"] <- "molecular_id"

  finemap_res <- ctwas::anno_finemap_res(finemap_res,
                                         snp_map = snp_map,
                                         mapping_table = gene_annotation,
                                         add_gene_annot = TRUE,
                                         map_by = "molecular_id",
                                         drop_unmapped = FALSE,
                                         add_position = TRUE,
                                         use_gene_pos = "mid")

  return(finemap_res)
}

extract_finemap_df <- function(name, paths) {
  df <- readRDS(paths[[name]])
  df$pval <- ctwas::z2p(df$z)
  df$greater_group <- gsub(".RDS|.rds", "", name)
  return(df)
}


prepare_finemap_res <- function(input_path, snp_map) {

  paths <- extract_path(input_path, "/finemap_res")

  finemap_res_list <- lapply(names(paths),
                             extract_finemap_df,
                             paths = paths)

  # finemap_res <- dplyr::bind_rows(finemap_res_list)
  finemap_res <- data.table::rbindlist(finemap_res_list, use.names = TRUE, fill = TRUE)

  finemap_res <- data.frame(finemap_res)

  gtf <- prepare_gene_annotation(EnsDb.Hsapiens.v86, finemap_res)
  finemap_res <- add_gene_annotation(finemap_res, snp_map, gtf)

  return(finemap_res)
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

internal_lapply_herit_fun <- function(group_label, obj, gwas_n) {
  param <- readRDS(obj[[group_label]])

  ctwas_parameters <- ctwas::summarize_param(param,
                                             gwas_n,
                                             enrichment_test = "fisher")

  data <- extract_heritability(ctwas_parameters, gsub("\\.RDS", "", group_label))

  return(data)
}

herit_table <- function(input_path, gwas_n) {

  paths <- extract_path(input_path, "/param")

  dfs <- lapply(names(paths),
                FUN = internal_lapply_herit_fun,
                obj = paths,
                gwas_n = gwas_n)

  

  # df <- dplyr::bind_rows(df)
  df <- data.table::rbindlist(dfs, use.names = TRUE, fill = TRUE)

  return(df)
}

################################################################################
################################################################################

generate_frame <- function(name, param_sub_group_obj) {

  df <- data.table::data.table(name = names(param_sub_group_obj[[name]]),
                   score = unname(param_sub_group_obj[[name]]),
                   category = name)
  return(df)
}

summarise_parameters_single <- function(group_label, obj, gwas_n) {
  param <- readRDS(obj[[group_label]])

  ctwas_parameters <- ctwas::summarize_param(param,
                                             gwas_n,
                                             enrichment_test = "fisher")

  df <- data.table::rbindlist(

         lapply(names(ctwas_parameters),
                FUN = generate_frame,
                param_sub_group_obj = ctwas_parameters),

                use.names = TRUE, fill = TRUE
                )
  

  # gene names present in df (excluding SNP and NA)
  gene_names <- df[!is.na(name) & name != "SNP", unique(name)]

  if (length(gene_names) > 1L) {
    stop("Expected at most one tissue-context (non-SNP name) in df, found",
    length(gene_names), gene_names,
    ",")
  }

  tot_snp <- data.table::copy(df[category == "total_pve"])
  tot_snp[, name := "SNP"]
  tot_gene <- data.table::copy(df[category == "total_pve"])
  tot_gene[, name := gene_names]
  # total_pve rows expanded for each gene name

  df <- data.table::rbindlist(list(df[!is.na(name)], tot_snp, tot_gene))

  df[, `group_label` := gsub("\\.RDS$", "", group_label)]

  return(df)
}

summarise_parameters_iterative <- function(inputdir, gwas_n) {

  paths <- extract_path(inputdir, "/param")

  dfs <- lapply(names(paths),
                FUN = summarise_parameters_single,
                obj = paths,
                gwas_n = gwas_n)

  # df <- dplyr::bind_rows(dfs)
  dfl <- data.table::rbindlist(dfs, use.names = TRUE, fill = TRUE)

  dfl[, label :=
    fcase(
      name == "SNP", paste(name, group_label, sep = "_"),
      default = paste0("gene_", group_label)
    )
  ]

  dfw <- dcast(dfl, label ~ category, value.var = "score")
  return(list(`long` = df, `wide` = dfw))
}

################################################################################
################################################################################

# combined_pip_by_type <- combine_gene_pips(ctwas_obj[["no_model_eqtl_Adipose"]][["finemap_res"]],
#                                           # mapping_table = gtf,
#                                           group_by = "molecular_id",
#                                           by = "type",
#                                           method = "combine_cs",
#                                           filter_cs = FALSE,
#                                           include_cs_id = TRUE)

################################################################################
################################################################################

save_files <- function(List, study, out) {

  # Ensure trailing slash
  out <- if (endsWith(out, "/")) out else glue("{out}/")

  # Create study-specific output directory
  out_dir <- glue("{out}{study}")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Save each element of the list as TSV
  lapply(names(List), function(x) {
    file_path <- glue("{out_dir}/{x}.tsv")
    write.table(List[[x]], file=file_path, sep="\t", row.names=FALSE, quote=FALSE)
  })
}

################################################################################
################################################################################

find_gwas_n <- function(metatable_path, study, gwas_n) {

  if (!is.null(metatable_path) && !is.null(study)) {
    gwas_n <- extract_gwas_n(metatable_path, study)

  } else if (!is.null(gwas_n)) {
    gwas_n <- gwas_n

  } else {
    stop(
      "WARNING - Must provide EITHER metatable_path AND study name OR gwas_n"
    )
  }
  return(gwas_n)
}

################################################################################
################################################################################

main <- function(inputdir, snp_map_path, gwas_n) {

  snp_map <- readRDS(snp_map_path)

  finemap_res <- prepare_finemap_res(inputdir, snp_map)

  parameter_df_list <- summarise_parameters_iterative(inputdir, gwas_n)

  heritability_table <- herit_table(inputdir, gwas_n)

  # finemap_res_ss <- subset(finemap_res, group != "SNP" & susie_pip > 0.8 & !is.na(cs))
  message("Files processed")
  return( list(`finemap_res` = finemap_res,
               `parameters_long` = parameter_df_list[['long']],
               `parameters_wide` = parameter_df_list[['wide']],
               `heritability_table` = heritability_table)
  )
}

################################################################################
################################################################################

if (sys.nframe() == 0) {

  script_dir <- get_script_dir()
  source(file.path(script_dir, "utils.R"))

  suppressPackageStartupMessages({
    library("glue")
    library("dplyr")
    library("ctwas")
    library("enrichR")
    library("ggplot2")
    library("argparse")
    library("data.table")
    library("EnsDb.Hsapiens.v86")
  })

  ################################################################################
  ################################################################################

  parser <- ArgumentParser(description = "Prepare reference files for cTWAS")

  parser$add_argument("--inputdir",
                      help = "",
                      required = FALSE)

  parser$add_argument("--snp_map_path",
                      help = "",
                      required = FALSE)

  parser$add_argument("--metatable_path",
                      help = paste0("GWAS metatble containing ",
                                    "Sample.name and sample_size cols"))

  parser$add_argument("--study",
                      help = "",
                      required = FALSE,
                      default = NULL)

  parser$add_argument("--gwas_n",
                      help = "",
                      required = FALSE,
                      default = NULL,
                      type = "integer")

  parser$add_argument("--genome_version",
                      help = "genome_version",
                      required = TRUE)

  parser$add_argument("--outpath",
                      help = "outpath to .RDS files",
                      required = TRUE)

  ################################################################################
  ################################################################################

  args <- parser$parse_args()

  gwas_n <- find_gwas_n(metatable_path = args$metatable_path,
                        study = args$study,
                        gwas_n = args$gwas_n)

  output <- main(inputdir = args$inputdir,
                 snp_map_path = args$snp_map_path,
                 gwas_n)

  save_files(output, args$study, args$outpath)

  ################################################################################
  ################################################################################
}

# input_path <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/result.raw/ukb-d-30780_irnt"
# snp_map_path <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/data.preprocessed/snp_map/ukb_b38_0.1_chrom_all.rds"
# metatable_path <- "/home/acano/Bioadhoc/GWAS/data/meta_tables/hg38/45_diseases_metadata_original.csv"
# outpath <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/result.refined"

################################################################################
################################################################################
# a <- lapply(names(param), function(p) {
#     make_convergence_plots(p, gwas_n)})

# pdf("/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/summary/locusplot.pdf", width = 8, height = 10)
# lapply(param, function(p) {
#     make_convergence_plots(p, gwas_n)})
#  dev.off()