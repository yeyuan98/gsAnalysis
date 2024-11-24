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
