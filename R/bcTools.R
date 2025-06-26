# Helper functions acting on Bioconductor standard objects

#' Subset SummarizedExperiment
#'
#' 1. Subset columns using `expr_colData`.
#' 2. Evaluate `expr_assay` on `assays()` of the subset.
#' 3. Apply function `fn_rowwise` to get subset vector.
#' 4. Use the subset vector to subset rows of the input `se`.
#'
#' @param se SummarizedExperiment object.
#' @param expr_assay Expression to evaluate on assays
#' @param expr_colData Expression to evaluate on colData
#' @param fn_rowwise Predicate to define rowwise filter
#'
#' @returns Subset SummarizedExperiment object.
#' @export
#'
#' @examples
#' # TODO
se_subset <- function(se, expr_assay, expr_colData, fn_rowwise){
  expr_colData <- substitute(expr_colData)
  expr_assay <- substitute(expr_assay)

  # colData subsetting
  colData <- as.data.frame(SummarizedExperiment::colData(se))
  sel_col <- eval(expr_colData, envir = colData)
  if (!any(sel_col)) stop("colData subset is empty. Check expression.")
  se_sub <- se[, sel_col]

  # assay evaluation
  assays <- as.list(SummarizedExperiment::assays(se_sub))
  assay_eval <- eval(expr_assay, envir = assays)
  if (!is.matrix(assay_eval) | any(dim(assay_eval) != dim(se_sub)))
    stop("Assay expression must evaluate to a matrix of identical shape.")

  # apply predicate rowwise function
  pd <- apply(assay_eval, 1, fn_rowwise)
  if (!is.logical(pd)) stop("fn_rowwise must be predicate.")

  # return the subset
  return(se[pd,])
}
