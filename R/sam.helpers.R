# Functions for SAM/BAM processing downstream of Rsamtools

#' Summarize CIGAR string by addition
#'
#' @param cigar CIGAR strings as character vector
#' @param which CIGAR operators to extract
#' @param FUNC Summarize function. Accept numeric vector and output length = 1.
#'
#' @return Numeric vector of sums of specified CIGAR operators.
#' @export
#'
#' @examples
#' cigar.test <- c(
#'   "148H47M1D113M1D34M1D39M460H",
#'   "50M"
#' )
#' bam_summary_cigar(
#'   # All Ops that consume reference
#'   # Therefore sum = reference lengths
#'   cigar.test, which = c("M", "D", "N", "=", "X")
#' )
bam_summary_cigar <- function(cigar, which, FUNC = sum){
  # Construct regex and lookup
  which <- paste(which, collapse = "|")
  which <- paste0("[0-9]+(?=", which, ")")
  cigar <- stringr::str_extract_all(cigar, which)
  # Convert and add
  for (i in 1:length(cigar)){
    curr <- as.numeric(cigar[[i]])
    cigar[[i]] <- FUNC(curr)
  }
  # Return
  return(unlist(cigar))
}

