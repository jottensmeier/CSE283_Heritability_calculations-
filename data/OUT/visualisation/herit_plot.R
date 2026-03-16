library("DBI")
library("RSQLite")
library("ggplot2")
library("data.table")

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

##############################
##### PATH CONSTRUCTORS ######
##############################

get_D_S_C_name_from_ <- function(paths, r1 = 2, r2 = 0) {

  paste_helper <- function(el, r1, r2) {
        LEN <- length(el)
        LEN1 <- LEN - r1
        LEN2 <- LEN - r2
        el <- paste(el[LEN1:LEN2], collapse="|")
        return(el)}

  split <- strsplit(tools::file_path_sans_ext(paths), "/")

  names <- lapply(split,
                  paste_helper,
                  r1 = r1,
                  r2 = r2)

  return(unlist(names))
}

get_paths <- function(path_upper, ctwas_modality, file = "finemap_res") {

  if (!ctwas_modality %in% c("single", "multi")) {
    stop("Error - CTWAS modality does not exist - must be - single or multi")
  }

  paths <- list.files(file.path(path_upper, ctwas_modality),
                      recursive = TRUE,
                      full.names = TRUE)

  paths <- paths[ grepl(file, paths) &
                  grepl("summary",     paths) ]

  names(paths) <- get_D_S_C_name_from_(paths,
                                       r1 = 2,
                                       r2 = 1)

  return( paths )

}

# script_dir <- get_script_dir()
# source(file.path(script_dir, "utils.R"))

generate_plot <- function(df,
                          x = "context",
                          y = "N",
                          color = "context",
                          facet = "DISEASE_CAT") {

  context_colors <- c(
    Myeloid = "#B22222",
    NK = "#FFD700",
    B = "#32CD32",
    CD4 = "#00BFFF",
    CD8 = "#EE82EE",
    mashr_Adipose_Subcutaneous = "#8000FF",
    mashr_Liver = "#BFBFBF",
    mashr_Whole_Blood = "#164233",
    `mashr_Cells_EBV-transformed_lymphocytes` = "#FFD700",
    mashr_Lung = "#D3B2CB"
  )

  plt <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[color]])) +
    stat_summary(
      fun = mean,
      geom = "bar"
    ) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.2
    ) +
    geom_jitter(
      width = 0.2,
      height = 0,
      color = "black",
      size = 2
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = context_colors, drop = FALSE) +
    facet_wrap(vars(.data[[facet]]), scales = "free_y", ncol = 5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Tissue Type",
      y = "log enrichment - (πG/πV), (where πV is prior probability for a variant being causal)",
      title = "Enrichment of gene effects over non-gene-genetic effects."
    )

  return(plt)
}

prepare_plot_data <- function(df) {

  context_order <- c(
    "Myeloid",
    "NK",
    "B",
    "CD4",
    "CD8",
    "mashr_Adipose_Subcutaneous",
    "mashr_Liver",
    "mashr_Whole_Blood",
    "mashr_Cells_EBV-transformed_lymphocytes",
    "mashr_Lung"
  )

  df <- data.table::copy(df)

  df[, context := sub("^gene_", "", label)]
  df[, PROJECT := gsub("\\s+", "", PROJECT)]
  # df[, DISEASE_CAT := paste(disease, PROJECT, sep = "_")]
  df[, context := factor(context, levels = context_order)]

  return(df)
}

prepare_dice_param_table <- function(path, ctwas_modality, project, file) {
  
  params_dfs <- get_paths(path,
                             ctwas_modality = ctwas_modality,
                             file = file)

  dfs_D <- data.table::rbindlist(
    lapply(params_dfs, data.table::fread),
    idcol = "name" )
    dfs_D[, PROJECT := project ]

  return(dfs_D)

}


generate_bar_plot <- function(dfs, out) {

  dfs <- prepare_plot_data(dfs)

  plt <- generate_plot(dfs,
                          x = "context",
                          y = "score",
                          color = "context",
                          facet = "study_key")


  pdf(out, width = 14, height = 8)
  print(plt)
  dev.off()

}

prepare_table <- function(ctwas_path) {
      dfs_s <- prepare_dice_param_table(ctwas_path,
                                    "single",
                                    "single",
                                    "parameters_long")

    names(dfs_s) <- c("group", "name", "score", "category", "group_label", "label", "PROJECT")

    dfs_s <- dfs_s[name != "SNP"]
    dfs_s_h <- dfs_s[category == "prop_heritability"]


  dfs_m <- prepare_dice_param_table(ctwas_path,
                                    "multi",
                                    "multi",
                                    "parameters_long")
  names(dfs_m) <- c("group", "name", "score", "category", "group_label", "label", "PROJECT")

    dfs_m <- dfs_m[name!= "SNP"]
    dfs_m_h <- dfs_m[category == "prop_heritability"]

    return(
        list(
            `dfs_s` = dfs_s,
            `dfs_m` = dfs_m
        )
    )

}


plot_box <- function(df, colours) {

    ggplot(df, aes(x=PROJECT, y = score, fill = group)) +
    geom_boxplot(
      width = 0.2,
      alpha = 0.5,
      outlier.shape = NA,
      color = "black"
    ) +

    geom_jitter(
      aes(shape = PROJECT),
      width = 0.15,
      alpha = 0.6,
      size = 1.5
    ) +

    facet_wrap(~ group) +
    scale_fill_manual(values = colours, drop = FALSE) +
    scale_shape_manual(values = c(multi = 16, single = 17), drop = FALSE) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    labs(
      x = "cTWAS modality",
      y = "Percentage Heritability mediated by gene expression",
      title = "Distribution of GWAS SNP counts before and after harmonisation",
      shape = "Modality",
      fill = "Disease"
    )

}

plot_corr <- function(concat, colours) {

    df <- copy(concat)[category == "prop_heritability"]

    # reshape wide
    plot_df <- dcast(
    df,
    group + name ~ PROJECT,
    value.var = "score"
    )

    # optional: rename columns for clarity
    setnames(plot_df, old = c("single", "multi"),
            new = c("score_single", "score_multi"))

    # scatter plot
    ggplot(plot_df, aes(x = score_single, y = score_multi, color = group)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(~ group, ncol= 3) + 
    scale_color_manual(values = colours) + 
    labs(
        x = "Single-tissue heritability",
        y = "Combined-model heritability",
        title = "Single vs combined cTWAS heritability"
    ) +
    theme_classic()
}

plot_heat_map <- function(concat, mod) {

    heat_df <- concat[
  PROJECT == mod &
  category == "prop_heritability"
]

 plt <- ggplot(heat_df, aes(name, group, fill = score)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

  return(plt)
}


plot_top_10 <- function(concat,
                        disease_id = "LDL|ukb-d-30780_irnt",
                        top_n = 10,
                        order_by = "total") {
  library(data.table)
  library(ggplot2)

  df <- copy(concat)[
    group == disease_id &
    category == "prop_heritability" &
    PROJECT %in% c("single", "multi")
  ]

  df[, tissue := sub("^mashr_", "", name)]
  df[, tissue := sub("\\|.*$", "", tissue)]
  df[, tissue := gsub("_", " ", tissue)]

  wide <- dcast(
    df,
    tissue ~ PROJECT,
    value.var = "score",
    fun.aggregate = sum
  )

  if (!"single" %in% names(wide)) wide[, single := 0]
  if (!"multi" %in% names(wide))  wide[, multi := 0]

  wide[is.na(single), single := 0]
  wide[is.na(multi),  multi := 0]

  wide[, total := single + multi]

  if (order_by == "multi") {
    wide_top <- wide[order(-multi)][1:min(.N, top_n)]
    setorder(wide_top, multi)
  } else if (order_by == "single") {
    wide_top <- wide[order(-single)][1:min(.N, top_n)]
    setorder(wide_top, single)
  } else if (order_by == "total") {
    wide_top <- wide[order(-total)][1:min(.N, top_n)]
    setorder(wide_top, total)
  } else {
    stop("order_by must be 'single', 'multi', or 'total'")
  }

  wide_top[, tissue := factor(tissue, levels = tissue)]

  plot_df <- melt(
    wide_top,
    id.vars = "tissue",
    measure.vars = c("single", "multi"),
    variable.name = "PROJECT",
    value.name = "score"
  )

  plot_df[, PROJECT := factor(PROJECT, levels = c("single", "multi"))]

  ggplot(plot_df, aes(x = score, y = tissue, fill = PROJECT)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    labs(
      title = paste0("Top tissue contributions to heritability: ", disease_id),
      subtitle = paste0("Ordered by ", order_by, " contribution"),
      x = "Heritability contribution",
      y = NULL,
      fill = "Model"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      axis.text.y = element_text(size = 10)
    )
}

colours <- c(
  "IBD|ebi-a-GCST004131"  = "#1b9e77",
  "LDL|ukb-d-30780_irnt" = "#d95f02",
  "SBP|ukb-a-360"        = "#7570b3"
)


scatter <- function(ctwas_path) {

    dfs<- prepare_table(ctwas_path)
    dfs_s <- dfs[["dfs_s"]]
    dfs_m <- dfs[["dfs_m"]]


    concat <- data.table::rbindlist(list(dfs_m, dfs_s))
    concat$PROJECT <- factor(concat$PROJECT, levels = c("single", "multi"))
    plt <- plot_box(concat[category =="prop_heritability" ], colours)
    plt2 <- plot_corr(concat[category =="prop_heritability" ], colours)
    plt3 <- plot_heat_map(concat, "single")
    plt4 <- plot_heat_map(concat, "multi")
    plt5 <- plot_top_10(concat)

    pdf("/path/to/out/CSE284_presentation/box_plot.pdf",
    width = 12, height = 8)
    print(plt)
    print(plt2)
    print(plt3)
    print(plt4)
    print(plt5)
    dev.off()







}

# ctwas_path <- "/test_smk_single_49_gtex_49_subsetted_tissues_per_GWAS"




#   df_bind <- data.table::rbindlist(list(dfs_G, dfs_D))
#   setnames(df_bind, make.unique(names(df_bind)))
#   setnames(df_bind, c("name", "name.1"), c("study_key", "model_key"))

#   df_bind[, c("disease", "dataset") := tstrsplit(study_key, "|", fixed = TRUE)]
#   dfs <- df_bind[category == "log_enrichment"]
#   generate_bar_plot(dfs, paste(out, "log_enrichment_single_.pdf"))

#   dfs <- df_bind[category == "prop_heritability" & model_key != "SNP"]
#   generate_bar_plot(dfs, paste(out, "prop_heritability.pdf"))
