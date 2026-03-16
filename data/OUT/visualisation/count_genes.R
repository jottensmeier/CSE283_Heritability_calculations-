library("DBI")
library("ggplot2")
library("RSQLite")
library("data.table")

####################################################################
####################    COLOURS AND NAME MAPS   ####################
####################################################################


###################################################
####################    GTEX   ####################
###################################################

gtex_colours <- c(

# Adipose
Adipose_Subcutaneous = "#F4A582",
Adipose_Visceral_Omentum = "#D6604D",

# Adrenal
Adrenal_Gland = "#B2182B",

# Artery
Artery_Aorta = "#FDDBC7",
Artery_Coronary = "#F4A582",
Artery_Tibial = "#D6604D",

# Brain
Brain_Amygdala = "#6A3D9A",
Brain_Anterior_cingulate_cortex_BA24 = "#7B4AB2",
Brain_Caudate_basal_ganglia = "#8C5CC6",
Brain_Cerebellar_Hemisphere = "#9D70D8",
Brain_Cerebellum = "#AE84E8",
Brain_Cortex = "#6B9FE3",
Brain_Frontal_Cortex_BA9 = "#5C8ED6",
Brain_Hippocampus = "#4C7DC7",
Brain_Hypothalamus = "#3D6CB7",
Brain_Nucleus_accumbens_basal_ganglia = "#2E5CA8",
Brain_Putamen_basal_ganglia = "#1E4C97",
`Brain_Spinal_cord_cervical_c-1` = "#0F3C88",
Brain_Substantia_nigra = "#002F6C",

# Breast
Breast_Mammary_Tissue = "#F768A1",

# Cells
Cells_Cultured_fibroblasts = "#969696",
`Cells_EBV-transformed_lymphocytes` = "#636363",

# Colon
Colon_Sigmoid = "#1B9E77",
Colon_Transverse = "#66A61E",

# Esophagus
Esophagus_Gastroesophageal_Junction = "#D95F02",
Esophagus_Mucosa = "#E7298A",
Esophagus_Muscularis = "#7570B3",

# Heart
Heart_Atrial_Appendage = "#E41A1C",
Heart_Left_Ventricle = "#A50F15",

# Kidney
Kidney_Cortex = "#41AB5D",

# Liver
Liver = "#238B45",

# Lung
Lung = "#1F78B4",

# Salivary
Minor_Salivary_Gland = "#A6CEE3",

# Muscle
Muscle_Skeletal = "#FB9A99",

# Nerve
Nerve_Tibial = "#CAB2D6",

# Ovary
Ovary = "#FF7F00",

# Pancreas
Pancreas = "#FDBF6F",

# Pituitary
Pituitary = "#6A3D9A",

# Prostate
Prostate = "#B15928",

# Skin
Skin_Not_Sun_Exposed_Suprapubic = "#B2DF8A",
Skin_Sun_Exposed_Lower_leg = "#33A02C",

# Small intestine
Small_Intestine_Terminal_Ileum = "#1B9E77",

# Spleen
Spleen = "#762A83",

# Stomach
Stomach = "#E08214",

# Testis
Testis = "#1B7837",

# Thyroid
Thyroid = "#5E3C99",

# Uterus
Uterus = "#E7298A",

# Vagina
Vagina = "#C51B7D",

# Blood
Whole_Blood = "#2166AC"
)


group_colours <- c(
  `DICE-LUNG-combined-sc` = "#B22222",
  `DICE-LUNG-sc`          = "#FF8C00",
  `DICE-BLOOD-blk`        = "#FFD700",

  `GTEX-blk_enet`         = "#00BFFF",
  `GTEX-blk_mashr`        = "#00008B"
)

colours <- c( gtex_colours, group_colours)

####################################################################
# 1. Get file paths
####################################################################

get_paths <- function(path) {
    paths <- list.files(path, full.names=TRUE)
    paths <- paths[grepl("\\.db", paths)]

    Names <- lapply(strsplit(paths, "/"),
    FUN = function(x) {
        l <- length(x)
        return(
            tools::file_path_sans_ext(x[l])
            )
            }
            )

    names(paths) <- Names
return(paths)
}

####################################################################
# 2. Read db table
####################################################################

read_db <-  function(path) {
    con <- dbConnect(SQLite(), dbname = path)
    w <- dbReadTable(con, "weights")
    dbDisconnect(con)
    return(data.table(w))
}

####################################################################
# 3. manipulate dataframe - groupby counts
####################################################################

SNPs_per_gene <- function(df) {

  df[, length(unique(rsid)), by = gene]

}

####################################################################
# 4. Plot barplot 
####################################################################

plot_gene_c_bar <- function(df, colours) {

    ggplot(df, aes(x=celltype, y=V1, fill = celltype)) +
    geom_col() + 
    scale_fill_manual(values = colours) +
    geom_text(aes(label = V1), vjust = -0.3) +
    labs(
        x = "Cell type",
        y = "N unique genes",
        title = "N genes in prediction models"
    ) + 
    theme_bw() +
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold") 
    )
}

####################################################################
# 5. Plot histogram 
####################################################################

plot_hist <- function(df, colours) {

    ggplot(df, aes(x = V1, fill = celltype)) + 
      geom_histogram(bins = 30, color = "black") + 
      facet_wrap(~ celltype, scales = "free_y", ncol = 5) + 
      scale_fill_manual(values = colours) +
      labs(
        x = "Value",
        y = "Count",
        title = "Distribution of variants per gene - per celltype"
      ) +
      theme_bw() +
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold") 
    )

}

####################################################################
# 6. Plot Boxplot
####################################################################

plot_box <- function(df, out, colours) {

  ggplot(df, aes(x = group, y = V1)) +

    geom_violin(
      aes(fill = group),
      outlier.shape = NA,
      trim = FALSE,
      alpha = 0.5
    ) +

    geom_boxplot(
      width = 0.2,
      alpha = 0.5,
      outlier.shape = NA,
      color = "black"
    ) +

    geom_jitter(
      width = 0.15,
      alpha = 0.6,
      size = 1.5
    ) +

    scale_fill_manual(values = colours) +
    scale_color_manual(values = colours) +
    guides(color = "none") +
    theme_bw() +
    labs(
      x = "Model dataset",
      y = "Genes for testing",
      title = "Distribution of genes across prediction models"
    )
}
####################################################################
# 7. Generate counts table
####################################################################

generat_counts_table_iterative <- function(path) {
    
    paths <- get_paths(path)
    dfs <- lapply(paths, read_db)
    dfs <- lapply(dfs, SNPs_per_gene)
    snps_per_gene <- data.table::rbindlist(dfs, idcol = "celltype")
    snps_per_gene[, celltype := gsub("en_", "", celltype)]
    ngenes_per_cell <- snps_per_gene[ , length(unique(gene)), by = celltype]
    return(
        list(
            `snps_per_gene` = snps_per_gene,
            `ngenes_per_cell` = ngenes_per_cell
        )
    )

}

main <- function(path, out, colours) {

    paths <- get_paths(path)
    dfs <- lapply(paths, read_db)
    dfs <- lapply(dfs, SNPs_per_gene)
    snps_per_gene <- data.table::rbindlist(dfs, idcol = "celltype")
    snps_per_gene[, celltype := gsub("en_", "", celltype)]
    ngenes_per_cell <- snps_per_gene[ , length(unique(gene)), by = celltype]
    # return(
    #     list(
    #         `snps_per_gene` = snps_per_gene,
    #         `ngenes_per_cell` = ngenes_per_cell
    #     )
    # )

    plt1 <- plot_gene_c_bar(ngenes_per_cell, colours)
    plt2 <- plot_hist(snps_per_gene, colours)

    pdf(out, width = 15, height = 8)
    print(plt1)
    print(plt1 + coord_flip())
    print(plt2)
    print(plt2 + coord_flip())
    dev.off()

    return(ngenes_per_cell)
}


main_2_compare_gene_counts <- function(pathlist,  out, colours) {

    tables <- lapply(
        pathlist, function(path) 
    {generat_counts_table_iterative(path)[['ngenes_per_cell']]}
    )
    
    c_table <- data.table::rbindlist(tables, idcol = "group")

    npoints <- c_table[, .N, by = group]
    labels <- setNames(
    paste0(npoints$group, "\n(n=", npoints$N, ")"),
    npoints$group
    )

    plt <- plot_box(c_table, out, colours)
    plt <- plt + coord_cartesian(ylim = c(0, 20000)) + scale_x_discrete(labels = labels)
    pdf(out, width = 10, height = 5)
    print(plt)
    dev.off()
}


if (sys.nframe() == 0) {

library(argparse)

parser <- ArgumentParser()

parser$add_argument(
"--pathlist",
nargs = "+",
required = TRUE,
help = "List of input file paths"
)

parser$add_argument(
"--out",
required = TRUE,
help = "Output file path"
)

parser$add_argument(
"--colours",
nargs = "+",
required = TRUE,
help = "Vector of colours"
)

args <- parser$parse_args()


   main_2_compare_gene_counts(args$pathlist,  file.path(args$out, "PrediXcan_counts_genes_per_model_per_dataset.pdf"), args$colours)

}
