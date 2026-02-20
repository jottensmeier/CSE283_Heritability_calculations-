
load_weights <- function( weight_file,
                          weight_format = c("PredictDB", "FUSION"),
                          filter_protein_coding_genes = TRUE,
                          load_predictdb_LD = TRUE,
                          fusion_method = c("lasso", "enet", "top1", "blup", "bslmm", "best.cv"),
                          fusion_genome_version = NA,
                          ncore = 1) {

    weight_format <- match.arg(weight_format)
    fusion_method <- match.arg(fusion_method)

    if (weight_format == "PredictDB") {
        res <- load_predictdb_weights(weight_file,
                                      filter_protein_coding_genes = filter_protein_coding_genes, 
                                      load_predictdb_LD = load_predictdb_LD)
    }
    else if (weight_format == "FUSION") {
        res <- load_fusion_weights(weight_file, fusion_method = fusion_method,
            fusion_genome_version = fusion_genome_version, ncore = ncore)
    }
    return(res)
}


load_predictdb_weights <- function(weight_file,
                                   filter_protein_coding_genes = TRUE,
                                   load_predictdb_LD = TRUE) {

    stopifnot(file.exists(weight_file))
    loginfo("Load PredictDB weights")
    sqlite <- dbDriver("SQLite")
    db <- dbConnect(sqlite, weight_file)
    query <- function(...) dbGetQuery(db, ...)
    weight_table <- query("select * from weights")
    extra_table <- query("select * from extra")
    loginfo("Number of molecular traits in weights: %s", length(unique(weight_table$gene)))
    if (filter_protein_coding_genes) {
        if ("protein_coding" %in% extra_table$gene_type) {
            # loginfo("Limit to protein coding genes")
            extra_table <- extra_table[extra_table$gene_type == "protein_coding", , drop = FALSE]
            weight_table <- weight_table[weight_table$gene %in% extra_table$gene, ]
            # loginfo("Number of molecular traits in weights after filtering protein coding genes: %s",
                # length(unique(weight_table$gene)))
        }
        else {
            # loginfo("No 'protein_coding' in 'extra_table$gene_type'. Skipped filtering protein coding genes.")
        }
    }
    if (anyNA(weight_table)) {
        logwarn("weight_table contains NAs!")
    }
    if (load_predictdb_LD) {
        # loginfo("Load PredictDB LD")
        print(predictdb_LD_file)
        predictdb_LD_file <- paste0(file_path_sans_ext(weight_file),
            ".txt.gz")
        if (!file.exists(predictdb_LD_file)) {
            stop(paste("PredictDB LD file", predictdb_LD_file,
                "does not exist!"))
        }
        cov_table <- as.data.frame(fread(predictdb_LD_file, header = TRUE))
        if (anyNA(cov_table)) {
            logwarn("cov_table contains NAs!")
        }
    }
    else {
        cov_table <- NULL
    }
    dbDisconnect(db)
    return(list(weight_table = weight_table,
                extra_table = extra_table,
                cov_table = cov_table))
}

process_weights <- function(molecular_id,
                            type,
                            context,
                            weight_name,
                            weight_table,
                            cov_table,
                            snp_info,
                            weight_format = c("PredictDB", "FUSION"),
                            top_n_snps = NULL,
                            drop_strand_ambig = TRUE,
                            scale_predictdb_weights = TRUE,
                            verbose = FALSE) {

  if (weight_format == "FUSION") {
    scale_predictdb_weights <- FALSE
  }

  # Check LD reference SNP info
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info))){
    stop("snp_info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  if (verbose)
    loginfo("Processing weight for %s...", molecular_id)

  g.weight_table <- weight_table[weight_table$gene==molecular_id,]
  wgt.matrix <- as.matrix(g.weight_table[, "weight", drop = FALSE])
  rownames(wgt.matrix) <- g.weight_table$rsid

  chrom <- sapply(strsplit(g.weight_table$varID, "_"), "[[", 1)
  chrom <- unique(parse_number(chrom))
  if (length(chrom) > 1) {
    stop(sprintf("More than one 'chrom' in weight of %s!", molecular_id))
  }

  snp_idx <- match(g.weight_table$rsid, snp_info$id)
  snp_pos <- as.integer(snp_info$pos[snp_idx])
  snp_chrom <- snp_info$chrom[snp_idx]
  if (any(unique(snp_chrom) != chrom)) {
    stop(sprintf("'chrom' in weight of %s does not match with the LD reference!", molecular_id))
  }

  wgt.snpinfo <- data.frame(chrom = chrom,
                            id = g.weight_table$rsid,
                            cm = 0,
                            pos = snp_pos,
                            alt = g.weight_table$eff_allele,
                            ref = g.weight_table$ref_allele,
                            stringsAsFactors = FALSE)

  if (verbose)
    loginfo("Harmonize weight")
  harmonized_wgt_res <- harmonize_weights(wgt.matrix,
                                          wgt.snpinfo,
                                          snp_info,
                                          drop_strand_ambig = drop_strand_ambig)
  wgt.matrix <- harmonized_wgt_res$wgt.matrix
  wgt.snpinfo <- harmonized_wgt_res$wgt.snpinfo
  wgt.matrix <- wgt.matrix[wgt.matrix[, "weight"] != 0, , drop = FALSE]
  wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), , drop = FALSE]

  wgt.snp_ids <- intersect(rownames(wgt.matrix), snp_info$id)
  wgt <- wgt.matrix[match(wgt.snp_ids, rownames(wgt.matrix)), "weight", drop = FALSE]

  # use top n SNPs in weights
  if (!is.null(top_n_snps)) {
    wgt <- wgt[order(-abs(wgt[,"weight"])), ]
    wgt <- head(wgt,top_n_snps)
    wgt.snp_ids <- rownames(wgt)
  }

  # scale weights by variance from LD reference
  if (weight_format == "PredictDB") {
    if (scale_predictdb_weights){
      if (verbose) {
        loginfo("Scale weights by variance from LD reference")
      }
      wgt_snp_var <- snp_info$variance[match(wgt.snp_ids, snp_info$id)]
      wgt <- wgt*sqrt(wgt_snp_var)
    }
  }

  wgt.snpinfo <- wgt.snpinfo[match(wgt.snp_ids, wgt.snpinfo$id),]
  g.weight_table <- g.weight_table[match(wgt.snp_ids, g.weight_table$rsid),]

  n_wgt <- nrow(wgt)
  if (n_wgt > 0) {
    p0 <- min(wgt.snpinfo$pos, na.rm = TRUE)
    p1 <- max(wgt.snpinfo$pos, na.rm = TRUE)

    # Add LD matrix of weights (R_wgt)
    if (!is.null(cov_table)) {
      # compute LD of variants in weights using PredictedDB cov_table
      if (verbose) {
        loginfo("Compute LD of variants in weights (R_wgt) from PredictedDB LD")
      }
      g.cov_table <- cov_table[cov_table$GENE == molecular_id,]
      R_wgt <- get_weight_LD_from_predictdb(g.cov_table,
                                            g.weight_table,
                                            convert_cov_to_cor = TRUE)
      R_wgt <- R_wgt[wgt.snp_ids, wgt.snp_ids, drop=FALSE]
    } else{
      R_wgt <- NULL
    }
  } else {
    p0 <- p1 <- NA
    wgt <- R_wgt <- NULL
  }

  return(list("chrom" = chrom,
              "p0" = p0,
              "p1" = p1,
              "wgt" = wgt,
              "R_wgt" = R_wgt,
              "molecular_id" = molecular_id,
              "weight_name" = weight_name,
              "type" = type,
              "context" = context,
              "n_wgt" = n_wgt))
}

harmonize_weights <- function(wgt.matrix, wgt.snpinfo, snp_info, drop_strand_ambig = TRUE) {
    target_header <- c("chrom", "id", "pos", "alt", "ref")
    if (!all(target_header %in% colnames(snp_info))) {
        stop("snp_info needs to contain the following columns: ",
            paste(target_header, collapse = " "))
    }
    wgt.snpinfo <- wgt.snpinfo[match(rownames(wgt.matrix), wgt.snpinfo$id),
        ]
    snpnames <- intersect(wgt.snpinfo$id, snp_info$id)
    if (length(snpnames) != 0) {
        snps.idx <- match(snpnames, wgt.snpinfo$id)
        ld.idx <- match(snpnames, snp_info$id)
        qc <- allele.qc(wgt.snpinfo[snps.idx, ]$alt, wgt.snpinfo[snps.idx,
            ]$ref, snp_info[ld.idx, ]$alt, snp_info[ld.idx, ]$ref)
        ifflip <- qc[["flip"]]
        ifremove <- !qc[["keep"]]
        flip.idx <- snps.idx[ifflip]
        wgt.snpinfo[flip.idx, c("alt", "ref")] <- wgt.snpinfo[flip.idx,
            c("ref", "alt")]
        wgt.matrix[flip.idx, ] <- -wgt.matrix[flip.idx, ]
        if (drop_strand_ambig && any(ifremove)) {
            remove.idx <- snps.idx[ifremove]
            wgt.snpinfo <- wgt.snpinfo[-remove.idx, , drop = FALSE]
            wgt.matrix <- wgt.matrix[-remove.idx, , drop = FALSE]
        }
    }
    return(list(wgt.matrix = wgt.matrix, wgt.snpinfo = wgt.snpinfo))
}

allele.qc <- function (a1, a2, ref1, ref2) {

    a1 = toupper(a1)
    a2 = toupper(a2)
    ref1 = toupper(ref1)
    ref2 = toupper(ref2)
    ref = ref1
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip1 = flip
    ref = ref2
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip2 = flip
    snp = list()
    snp[["keep"]] = !((a1 == "A" & a2 == "T") | (a1 == "T" & 
        a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == 
        "C"))
    snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & 
        a2 == flip1)
    return(snp)
}

get_weight_LD_from_predictdb <- function (g.cov_table, g.weight_table, convert_cov_to_cor = TRUE) {
    if (convert_cov_to_cor) {
        g.cov_table <- convert_predictdb_cov_to_cor(g.cov_table)
    }
    if (any(grepl("rs", g.cov_table$RSID1))) {
        g.weight_table$varID <- g.weight_table$rsid
    }
    if (!all(g.weight_table$varID %in% unique(c(g.cov_table$RSID1, 
        g.cov_table$RSID2)))) {
        stop("Not all variants in weight_table are in cov_table!")
    }
    varIDs <- g.weight_table$varID
    n_wgt <- length(varIDs)
    if (n_wgt == 0) {
        return(NULL)
    }
    R_wgt <- diag(n_wgt)
    if (n_wgt > 1) {
        snp_pairs <- combn(length(varIDs), 2)
        R_snp_pairs <- apply(snp_pairs, 2, function(x) {
            match.idx <- which(g.cov_table$RSID1 == varIDs[x[1]] & 
                g.cov_table$RSID2 == varIDs[x[2]])
            if (length(match.idx) == 0) {
                match.idx <- which(g.cov_table$RSID1 == varIDs[x[2]] & 
                  g.cov_table$RSID2 == varIDs[x[1]])
            }
            if (length(match.idx) > 0) {
                g.cov_table[match.idx, "VALUE"]
            }
            else {
                NA
            }
        })
        R_wgt[t(snp_pairs)] <- R_snp_pairs
        R_wgt[t(snp_pairs[c(2, 1), ])] <- R_snp_pairs
    }
    rownames(R_wgt) <- g.weight_table$rsid
    colnames(R_wgt) <- g.weight_table$rsid
    return(R_wgt)
}

convert_predictdb_cov_to_cor <- function (cov_table) 
{
    stdev_table <- cov_table[cov_table$RSID1 == cov_table$RSID2, 
        ]
    stdev_table <- setNames(sqrt(stdev_table$VALUE), stdev_table$RSID1)
    cov_table$VALUE <- cov_table$VALUE/(stdev_table[cov_table$RSID1] * 
        stdev_table[cov_table$RSID2])
    return(cov_table)
}

library(logging)
library(readr)
weight_format <- "PredictDB"
verbose <- T
drop_strand_ambig <- F
top_n_snps <- NULL
scale_predictdb_weights = FALSE
weight_name <- "TMP"
type <- 'eqtl'
context <- 'something'

molecular_id <- weight_molecular_ids[[1]]

weights <- ctwas:::mclapply_check(weight_molecular_ids, function(molecular_id){
    process_weights(molecular_id,
                    type = type,
                    context = context,
                    weight_name = weight_name,
                    weight_table = weight_table,
                    cov_table = cov_table,
                    snp_info = snp_info,
                    weight_format = weight_format,
                    top_n_snps = top_n_snps,
                    drop_strand_ambig = drop_strand_ambig,
                    scale_predictdb_weights = scale_predictdb_weights,
                    verbose = verbose)
  }, mc.cores = ncore)

ncore=1



##################################################
##################################################
##################################################
compute_weight_LD_from_ref <- function (weights,
                                        region_info,
                                        LD_map,
                                        snp_map = NULL,
                                        LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                                        LD_loader_fun = NULL,
                                        snpinfo_loader_fun = NULL,
                                        ncore = 1) {
  # Define a function that will attach (region-specific) LD submatrices to each weight model.
  # Inputs:
  #   - weights: list of weight objects; each element must carry metadata (chrom/p0/p1/molecular_id/weight_name/type/context)
  #             and a SNP-by-weight matrix weights[[...]]$wgt whose rownames are SNP IDs.
  #   - region_info: data.frame defining LD regions (chrom, start, stop, region_id).
  #   - LD_map: data.frame mapping region_id -> LD_file (and optionally SNP_file) for the LD reference.
  #   - snp_map: optional pre-loaded SNP info by region (preferred when you already have SNP IDs / RSIDs without reading files).
  #   - LD_format: format of LD matrices on disk.
  #   - LD_loader_fun: optional custom loader for LD if LD_format == "custom".
  #   - snpinfo_loader_fun: optional custom loader for SNP info if SNP files are custom.
  #   - ncore: number of cores for parallel region loading.

    LD_format <- match.arg(LD_format)
    # Validate LD_format is one of the allowed strings; prevents typos and ensures downstream loaders know what to expect.

    if (!inherits(weights, "list"))
        stop("'weights' should be a list!")
    # Enforce that weights is a list (the function expects to iterate over named entries).

    if (!inherits(LD_map, "data.frame"))
        stop("'LD_map' should be a data frame!")
    # Enforce that LD_map is tabular; later code assumes column access like LD_map$region_id, LD_map$LD_file, LD_map$SNP_file.

    weight_info <- lapply(names(weights), function(x) {
        as.data.frame(weights[[x]][c("chrom", "p0", "p1", "molecular_id",
            "weight_name", "type", "context")])
    })
    # For each weight model (each list element in weights):
    #   - subset the model metadata columns that define genomic location and model identity
    #   - coerce to data.frame
    # Result: weight_info is a list of small data.frames (one per model).

    
    weight_info <- do.call(rbind, weight_info)
    # Row-bind (stack) the list of per-model metadata frames into a single data.frame:
    #   one row per model, with columns chrom/p0/p1/molecular_id/weight_name/type/context.

    weight_info$weight_id <- paste0(weight_info$molecular_id, "|", weight_info$weight_name)
    # Create a unique identifier for each model by concatenating molecular_id and weight_name.
    # This becomes the key used later to index into weights (weights[[weight_id]]).

    for (k in 1:nrow(weight_info)) {
        # Iterate over each model row in weight_info to assign which LD region(s) it overlaps.

        chrom <- weight_info[k, "chrom"]
        # Extract chromosome for this model.

        p0 <- weight_info[k, "p0"]
        # Extract start coordinate for this model’s window.

        p1 <- weight_info[k, "p1"]
        # Extract end coordinate for this model’s window.

        idx <- which(region_info$chrom == chrom & region_info$start <=
            p1 & region_info$stop > p0)
        # Identify rows of region_info whose region intervals overlap the model interval.
        # Overlap criteria (same chrom) AND (region start <= model end) AND (region stop > model start).

        weight_info[k, "region_id"] <- paste(sort(region_info[idx,
            "region_id"]), collapse = ",")
        # Assign the overlapping region_id(s) to this model.
        # If multiple regions overlap, sort them and collapse into a comma-separated string "r12,r13".
    }

    chrs <- sort(unique(weight_info$chrom))
    # Determine which chromosomes are represented among the models; process per chromosome.

    for (b in chrs) {
        # Loop over chromosomes to limit memory and I/O to one chromosome’s worth of regions at a time.

        loginfo("Computing LD for variants in weights on chr%s",
            b)
        # Emit a log line for progress monitoring.

        weightinfo <- weight_info[weight_info$chrom == b, ]
        # Subset the model metadata to just models on chromosome b.

        if (nrow(weightinfo) > 0) {
            # Safety check: only proceed if there are models on this chromosome.

            weight_region_ids <- names(sort(-table(weightinfo$region_id)))
            # Compute unique region-set strings used by models on this chromosome, ordered by frequency.
            # table(weightinfo$region_id) counts models per region-set (e.g., "r12" or "r12,r13").
            # -table(...) negates counts so sort() yields descending order.
            # names(...) extracts the region-set strings.

            weight_LD_list <- mclapply_check(weight_region_ids,
                function(x) {
                  # For one region-set string x (e.g., "r12,r13"):
                  #   - load the relevant LD matrix/matrices
                  #   - load/derive the SNP IDs for that LD
                  #   - subset LD to each model’s SNP list
                  #   - return a list mapping weight_id -> LD submatrix (R_wgt)

                  curr_region_LD_list <- list()
                  # Initialize an output container for this region-set: will store per-model LD submatrices.

                  curr_region_ids <- unlist(strsplit(x, ","))
                  # Split the region-set string into individual region IDs (vector).
                  # Example: "r12,r13" -> c("r12","r13").

                  curr_region_idx <- match(curr_region_ids, LD_map$region_id)
                  # Find the row indices in LD_map corresponding to each region ID.

                  LD_matrix_files <- unlist(strsplit(LD_map$LD_file[curr_region_idx], 
                    split = ","))
                  # Extract LD file path(s) for these regions.
                  # Each LD_map row may itself list multiple files (comma-separated), so split and flatten.

                  stopifnot(all(file.exists(LD_matrix_files)))
                  # Hard fail if any referenced LD file does not exist.

                  if (length(LD_matrix_files) > 1) {
                    # If multiple LD matrices are needed (e.g., multiple windows):
                    R_snp <- lapply(LD_matrix_files, load_LD,
                      format = LD_format, LD_loader_fun = LD_loader_fun)
                    # Load each LD matrix into a list (one per file).

                    R_snp <- suppressWarnings(as.matrix(bdiag(R_snp)))
                    # Combine matrices into one block-diagonal LD matrix:
                    #   assumes LD between SNP sets from different files is 0 (off-diagonal blocks are 0).
                    # Convert to a (dense) matrix; suppress warnings arising from coercion.
                  }
                  else {
                    # If only one LD file:
                    R_snp <- load_LD(LD_matrix_files, format = LD_format,
                      LD_loader_fun = LD_loader_fun)
                    # Load that LD matrix directly.
                  }

                  if (!is.null(snp_map)) {
                    # If SNP info is already supplied in memory (preferred when you’ve pre-harmonized IDs):
                    snpinfo <- as.data.frame(rbindlist(snp_map[curr_region_ids], 
                      idcol = "region_id"))
                    # Combine SNP info across the region IDs using rbindlist.
                    # idcol="region_id" adds a column identifying which region each SNP came from.
                    # snpinfo must include an 'id' column used below as the LD matrix dimnames.
                  }
                  else {
                    # Otherwise, load SNP info from files listed in LD_map:
                    SNP_info_files <- LD_map$SNP_file[curr_region_idx]
                    # Retrieve SNP-info file path(s) for the region(s).

                    stopifnot(all(file.exists(SNP_info_files)))
                    # Hard fail if SNP info files do not exist.

                    snpinfo <- read_snp_info_files(SNP_info_files,
                      snpinfo_loader_fun = snpinfo_loader_fun)
                    # Read SNP info; expected to return a data.frame with at least an 'id' column
                    # (the SNP IDs corresponding to rows/cols of the LD matrix).
                  }

                  rownames(R_snp) <- snpinfo$id
                  colnames(R_snp) <- snpinfo$id
                  # Assign SNP IDs as dimension names for the LD matrix.
                  # Critical: these IDs must match the SNP IDs used in weights[[...]]$wgt rownames
                  # for the subsetting below to work.

                  weight_ids <- weightinfo[weightinfo$region_id ==
                    x, "weight_id"]
                  # Identify which models (weight_id values) correspond to this exact region-set string x.

                  for (weight_id in weight_ids) {
                    # For each model in this region-set, subset the LD matrix to the SNPs used by that model:

                    wgt_snp_ids <- rownames(weights[[weight_id]]$wgt)
                    # Extract the SNP IDs used in the weight model (row names of the weight matrix).

                    R_wgt <- R_snp[wgt_snp_ids, wgt_snp_ids,
                      drop = FALSE]
                    # Subset the region LD matrix to the model’s SNP set (square submatrix).
                    # drop=FALSE preserves matrix shape even if only 1 SNP.

                    curr_region_LD_list[[weight_id]] <- R_wgt
                    # Store this model-specific LD submatrix in the output list keyed by weight_id.
                  }

                  curr_region_LD_list
                  # Return the list of {weight_id -> R_wgt} for this region-set.
                }, mc.cores = ncore, stop_if_missing = TRUE)
            # Run the above region-set processing in parallel over weight_region_ids.
            # mc.cores controls parallelism; stop_if_missing=TRUE likely enforces erroring on missing inputs.

            weight_LD_list <- unlist(weight_LD_list, recursive = FALSE)
            # Flatten the list-of-lists (one per region-set) into a single named list of matrices:
            #   names are weight_id, values are R_wgt matrices.

            for (weight_id in names(weight_LD_list)) {
                # Attach each computed LD submatrix back into the corresponding weights object:
                weights[[weight_id]][["R_wgt"]] <- weight_LD_list[[weight_id]]
                # Store model-specific LD as weights[[weight_id]]$R_wgt.
            }
        }
    }

    return(weights)
    # Return the original weights list with an added component R_wgt for each weight model.
}
