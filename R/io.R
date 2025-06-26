# Helpers to simplify IO of common file types

#' writeXStringSet with identifiers
#'
#' Write a XStringSet to file (FASTA/FASTQ) with names as identifiers.
#'
#' Somehow the original `Biostrings::writeXStringSet()` does not write
#' identifiers. This function uses names provided by `name()` and write ID.
#'
#' @param x Object XString to write to file
#' @param filepath File to write to
#' @param append Must be False (does not support append)
#' @param compress Must be False (does not support compression)
#' @param ... Forwarded to `Biostrings::writeXStringSet()`
#'
#' @export
#'
#' @examples
#' #TODO
writeXStringSetNamed <- function(x, filepath, append = FALSE, compress = FALSE, ...){
  # Check no append nor compress
  if (append | compress) stop("append/compress not supported.")
  # Get names
  nm <- names(x)
  if (is.null(nm)) stop("x must be named.")
  # First call Biostrings::writeXStringSet()
  Biostrings::writeXStringSet(x, filepath, ...)
  # Next amend output file by adding identifiers to the columns
  tf <- tempfile()
  ic <- file(filepath, "r")
  oc <- file(tf, "w")
  idx <- 1
  while ( TRUE ) {
    line = readLines(ic, n = 1)
    if ( length(line) == 0 ) break
    if (line == ">" | line == "@"){
      writeLines(paste0(line, nm[idx]), con = oc)
      idx <- idx+1
    } else{
      writeLines(line, con = oc)
    }
  }
  close(ic)
  close(oc)
  file.copy(from = tf, to = filepath, overwrite = TRUE)
  file.remove(tf)
  return(invisible())
}

#' Read STAR Aligner Gene Counts
#'
#' @param named_paths Named list of count table paths.
#' @param clean Boolean, whether to remove entries starting with "N_".
#'
#' @returns `SummarizedExperiment` of counts. Metadata `sample` will be names.
#' @export
#'
#' @examples
#' # TODO
readStarGeneCounts <- function(named_paths, clean = FALSE){

  # Check that sample names are unique
  samples <- names(named_paths)
  if (length(unique(samples)) != length(samples))
    stop("Sample names (list names) must be unique.")

  # Read in counts
  data <- lapply(seq_along(named_paths), function(idx){
    fp <- named_paths[[idx]]
    df <- readr::read_delim(
      fp, delim = "\t", col_types = "ciii", skip = 0,
      col_names = c("symbol", samples[idx], "count.F", "count.R"))
    return(df[,1:2])
  })

  # Join counts
  data <- Reduce(function(x,y) dplyr::full_join(x,y,by = "symbol"), data)

  # Convert into matrix
  data.mat <- data[,-1]
  data.mat <- as.matrix(data.mat)
  rownames(data.mat) <- data$symbol
  if (!all(colnames(data.mat) == samples))
    stop("Sample name doesn't match. Submit issue.")
  colnames(data.mat) <- NULL

  # Remove uninformative symbols if asked
  if (clean){
    deselect <-
      row.names(data.mat) %in%
      c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
    data.mat <- data.mat[!deselect,]
  }

  # Construct SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = data.mat), colData = data.frame(sample = samples)
  )
  return(se)
}
