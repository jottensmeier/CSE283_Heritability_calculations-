
function (region_metatable, snpinfo_loader_fun = NULL)
{
    region_metatable <- as.data.frame(region_metatable)
    target_header <- c("chrom", "start", "stop")
    if (!all(target_header %in% colnames(region_metatable))) {
        stop("region_metatable needs to contain the following columns: ",
            paste(target_header, collapse = " "))
    }
    if (is.null(region_metatable$LD_file)) {
        stop("Please provide filenames of LD matrices in the 'LD_file' column of region_metatable")
    }
    if (is.null(region_metatable$SNP_file)) {
        stop("Please provide filenames of variant info in the 'SNP_file' column of region_metatable")
    }
    region_info <- region_metatable
    region_info$chrom <- readr::parse_number(as.character(region_info$chrom))
    region_info$start <- as.numeric(region_info$start)
    region_info$stop <- as.numeric(region_info$stop)
    if (is.null(region_info$region_id)) {
        region_info$region_id <- paste0(region_info$chrom, "_",
            region_info$start, "_", region_info$stop)
    }
    rownames(region_info) <- NULL
    region_info <- region_info[, c("chrom", "start", "stop",
        "region_id")]
    stopifnot(all(file.exists(region_metatable$SNP_file)))
    if (!is.null(snpinfo_loader_fun)) {
        snp_map <- lapply(region_metatable$SNP_file, snpinfo_loader_fun)
    }
    else {
        snp_map <- lapply(region_metatable$SNP_file, read_snp_info_file)
    }
    names(snp_map) <- region_metatable$region_id
    stopifnot(all(file.exists(region_metatable$LD_file)))
    LD_map <- data.frame(region_id = region_metatable$region_id,
        LD_file = region_metatable$LD_file, SNP_file = region_metatable$SNP_file)
    return(list(region_info = region_info, snp_map = snp_map,
        LD_map = LD_map))
}
