#' Pair metadata builder for siamese training
#'
#' Constructs a metadata table with sequence id pairs and binary labels that
#' indicate whether two sequences share at least one GO term. This table is
#' later consumed by \code{make_pair_generator()} and keeps storage demands low
#' because it omits the one-hot encoded tensors.
#'
#' @param cc3_file Path to the RDS file that maps sequence ids to GO terms.
#' @param output_file Destination for the resulting table. Defaults to
#'   \code{/vol/projects/BIFO/genomenet/twin_data/rds/CC/pairs_CC_YYYYMMDD.rds}.
#'   The file extension determines the format (\code{rds} or \code{csv}).
#' @param seq_id_filter Optional character vector of sequence ids to keep.
#'   When provided, only entries whose names match this vector are used.
#' @param max_pairs Maximum number of unique pairs to sample. Set to
#'   \code{NULL} to enumerate every possible pair (not recommended for large
#'   inputs). Defaults to 2e6.
#' @param seed Optional seed used for reproducible sampling.
#' @param verbose Whether to emit progress messages.
#'
#' @examplesIf FALSE
#' pairs_df <- build_pair_metadata(
#'   cc3_file = "/vol/projects/BIFO/genomenet/twin_data/preprocessed/GO_CC_transposed.rds",
#'   max_pairs = 1e5,
#'   seed = 123
#' )
#' head(pairs_df)
#'
#' @returns Invisibly returns the in-memory data frame (mainly for testing).
#' @export
build_pair_metadata <- function(cc3_file,
                                output_file = NULL,
                                seq_id_filter = NULL,
                                max_pairs = 2e6,
                                seed = 42,
                                verbose = TRUE) {
  if (!file.exists(cc3_file)) {
    stop("cc3_file does not exist: ", cc3_file)
  }
  
  cc3 <- readRDS(cc3_file)
  if (!is.list(cc3) || is.null(names(cc3))) {
    stop("cc3 must be a named list where each element contains GO terms.")
  }
  
  if (!is.null(seq_id_filter)) {
    seq_id_filter <- unique(seq_id_filter)
    present_ids <- seq_id_filter[seq_id_filter %in% names(cc3)]
    missing_ids <- setdiff(seq_id_filter, present_ids)
    if (length(missing_ids) > 0 && verbose) {
      warning("Ignoring ", length(missing_ids),
              " ids not found in cc3 (showing first 5): ",
              paste(head(missing_ids, 5), collapse = ", "))
    }
    if (length(present_ids) < 2) {
      stop("seq_id_filter retained fewer than two sequences.")
    }
    cc3 <- cc3[present_ids]
  }
  
  seq_ids <- names(cc3)
  n_seq <- length(seq_ids)
  if (n_seq < 2) stop("At least two sequences are required.")
  
  go_terms <- lapply(cc3, unique)
  total_pairs <- n_seq * (n_seq - 1) / 2
  
  if (is.null(output_file)) {
    today <- format(Sys.Date(), "%Y%m%d")
    output_file <- file.path(
      "/vol/projects/BIFO/genomenet/twin_data/rds/CC",
      paste0("pairs_CC_", today, ".rds")
    )
  }
  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  file_ext <- tolower(tools::file_ext(output_file))
  if (!file_ext %in% c("rds", "csv")) {
    stop("output_file must end with .rds or .csv")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  if (is.null(max_pairs) || max_pairs >= total_pairs) {
    target_pairs <- total_pairs
    sampling_mode <- "all"
    if (verbose) {
      message("Generating all ", format(total_pairs, big.mark = ","), " pairs.")
    }
  } else {
    target_pairs <- as.integer(max_pairs)
    sampling_mode <- "sample"
    if (verbose) {
      message("Sampling ", format(target_pairs, big.mark = ","), " pairs ",
              "out of ", format(total_pairs, big.mark = ","), ".")
    }
  }
  
  seq_id1 <- character(target_pairs)
  seq_id2 <- character(target_pairs)
  labels <- integer(target_pairs)
  
  shared_term <- function(i, j) {
    any(go_terms[[i]] %in% go_terms[[j]])
  }
  
  if (sampling_mode == "all") {
    idx <- 1L
    for (i in seq_len(n_seq - 1L)) {
      terms_i <- go_terms[[i]]
      for (j in (i + 1L):n_seq) {
        seq_id1[idx] <- seq_ids[i]
        seq_id2[idx] <- seq_ids[j]
        labels[idx] <- as.integer(any(terms_i %in% go_terms[[j]]))
        idx <- idx + 1L
      }
    }
  } else {
    seen <- new.env(parent = emptyenv())
    idx <- 0L
    while (idx < target_pairs) {
      sampled <- sample.int(n_seq, 2L, replace = FALSE)
      i <- min(sampled)
      j <- max(sampled)
      key <- paste0(i, "-", j)
      if (!exists(key, envir = seen, inherits = FALSE)) {
        assign(key, TRUE, envir = seen)
        idx <- idx + 1L
        seq_id1[idx] <- seq_ids[i]
        seq_id2[idx] <- seq_ids[j]
        labels[idx] <- as.integer(shared_term(i, j))
        if (verbose && idx %% 100000 == 0L) {
          message("Sampled ", format(idx, big.mark = ","), " pairs...")
        }
      }
    }
  }
  
  pairs_df <- data.frame(
    seq_id1 = seq_id1,
    seq_id2 = seq_id2,
    label = labels,
    stringsAsFactors = FALSE
  )
  
  if (file_ext == "csv") {
    utils::write.csv(pairs_df, output_file, row.names = FALSE)
  } else {
    saveRDS(pairs_df, output_file)
  }
  
  if (verbose) message("Saved pairs to: ", output_file)
  invisible(pairs_df)
}

#' Streaming pair generator for keras twin networks
#'
#' Creates a lightweight R generator that yields batches of sequence pairs by
#' reading one-hot encoded sequences from a shared lookup table on demand.
#' The generator matches the signature expected by keras training loops that
#' consume a list of inputs \code{list(X1, X2)} and labels.
#'
#' @param pairs_file Path to the metadata file produced by
#'   \code{build_pair_metadata()}. Supports RDS or CSV formats.
#' @param xmap_file Path to the RDS file that maps sequence ids to one-hot
#'   encoded matrices. Accepts either a named list of tensors or the bundled
#'   format used by \code{HMP2_proteinfamilies_CC_2200_onehot.rds}
#'   (\code{list(nucseq = ..., seqname = ...)}).
#' @param batch_size Number of pairs to emit per call.
#' @param shuffle Whether to reshuffle the order after every pass through the
#'   metadata table.
#' @param seed Optional seed for reproducible shuffling.
#'
#' @returns A closure with no arguments. Each invocation returns
#'   \code{list(list(X1, X2), Y)} where \code{X1} and \code{X2} both have shape
#'   \code{[batch_size, 3190, 4]} (or whatever dimensions are stored in
#'   \code{xmap_file}) and \code{Y} is a \code{[batch_size, 1]} matrix.
#'
#' @examplesIf file.exists("pairs_sample.rds") &&
#'   file.exists("HMP2_proteinfamilies_CC_2200_onehot.rds")
#' gen <- make_pair_generator(
#'   pairs_file = "pairs_sample.rds",
#'   xmap_file = "HMP2_proteinfamilies_CC_2200_onehot.rds",
#'   batch_size = 256
#' )
#' batch <- gen()
#' str(batch)
#'
#' @export
make_pair_generator <- function(pairs_file,
                                xmap_file,
                                batch_size = 512,
                                shuffle = TRUE,
                                seed = 42) {
  if (!file.exists(pairs_file)) stop("pairs_file does not exist: ", pairs_file)
  if (!file.exists(xmap_file)) stop("xmap_file does not exist: ", xmap_file)
  
  pairs <- load_pairs_table(pairs_file)
  if (!all(c("seq_id1", "seq_id2", "label") %in% names(pairs))) {
    stop("pairs_file must contain columns: seq_id1, seq_id2, label")
  }
  
  xmap <- load_xmap_table(xmap_file)
  if (is.null(names(xmap)) || any(names(xmap) == "")) {
    stop("xmap_file must provide a named list of sequence tensors.")
  }
  
  missing_ids <- unique(c(
    setdiff(pairs$seq_id1, names(xmap)),
    setdiff(pairs$seq_id2, names(xmap))
  ))
  if (length(missing_ids) > 0) {
    stop("Sequence ids missing in xmap: ", paste(head(missing_ids, 5), collapse = ", "))
  }
  
  n_pairs <- nrow(pairs)
  if (n_pairs == 0) stop("pairs_file contains no pairs.")
  
  sample_indices <- seq_len(n_pairs)
  if (!is.null(seed)) set.seed(seed)
  if (shuffle) sample_indices <- sample(sample_indices)
  cursor <- 1L
  
  example_tensor <- canonicalize_sequence(xmap[[pairs$seq_id1[1]]])
  template_dim <- dim(example_tensor)
  if (is.null(template_dim) || length(template_dim) != 2) {
    stop("Sequence tensors must collapse to 2D [length, channels].")
  }
  
  fetch_batch <- function() {
    if (cursor > n_pairs) {
      if (shuffle) {
        sample_indices <<- sample(sample_indices)
      }
      cursor <<- 1L
    }
    end <- min(cursor + batch_size - 1L, n_pairs)
    batch_idx <- sample_indices[cursor:end]
    cursor <<- end + 1L
    
    current_batch <- length(batch_idx)
    x_dim <- c(current_batch, template_dim)
    x1 <- array(0, dim = x_dim)
    x2 <- array(0, dim = x_dim)
    y <- matrix(0, nrow = current_batch, ncol = 1)
    
    for (k in seq_len(current_batch)) {
      idx <- batch_idx[k]
      id1 <- pairs$seq_id1[idx]
      id2 <- pairs$seq_id2[idx]
      x1[k, , ] <- canonicalize_sequence(xmap[[id1]])
      x2[k, , ] <- canonicalize_sequence(xmap[[id2]])
      y[k, 1] <- pairs$label[idx]
    }
    
    list(list(x1, x2), y)
  }
  
  fetch_batch
}

load_pairs_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    utils::read.csv(path, stringsAsFactors = FALSE)
  } else {
    readRDS(path)
  }
}

load_xmap_table <- function(path) {
  raw <- readRDS(path)
  
  if (is.list(raw) &&
      all(c("nucseq", "seqname") %in% names(raw))) {
    seqs <- raw[["nucseq"]]
    ids <- raw[["seqname"]]
    ids <- normalize_seq_ids(ids)
    
    if (length(seqs) != length(ids)) {
      stop("Mismatch between number of sequences and seqname entries in xmap.")
    }
    if (anyDuplicated(ids)) {
      dup_ids <- ids[duplicated(ids)]
      stop("Duplicated seq_ids detected in xmap: ",
           paste(head(dup_ids, 5), collapse = ", "), " ...")
    }
    names(seqs) <- ids
    return(seqs)
  }
  
  if (is.list(raw) && !is.null(names(raw))) {
    return(raw)
  }
  
  stop("Unsupported xmap object. Expected a named list or nucseq/seqname bundle.")
}

normalize_seq_ids <- function(ids) {
  if (is.list(ids)) {
    ids <- unlist(ids, use.names = FALSE, recursive = TRUE)
  }
  if (!is.character(ids)) {
    stop("seqname entries in xmap must be character values.")
  }
  ids
}

canonicalize_sequence <- function(tensor) {
  dims <- dim(tensor)
  if (is.null(dims)) {
    stop("Sequence tensor must have dimensions.")
  }
  canonical <- drop(tensor)
  canonical_dims <- dim(canonical)
  if (is.null(canonical_dims) || length(canonical_dims) != 2) {
    stop("Sequence tensors must reduce to 2D [length, channels]; ",
         "got dims: ", paste(dims, collapse = "x"))
  }
  canonical
}

# Example (replace file paths with actual locations):
# gen <- make_pair_generator(
#   pairs_file = "pairs_sample.rds",
#   xmap_file = "HMP2_proteinfamilies_CC_2200_onehot.rds",
#   batch_size = 256
# )
# batch <- gen()
# str(batch)
