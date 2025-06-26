# TODO: TOO SLOW CURRENTLY.

# Somewhat customized differential analysis routines
#   beta-regression test - suitable for ratios

#' Differential analysis with beta regression
#'
#' This takes data.frame input.
#'
#' @param data Data to use.
#' @param design Design formula of the problem.
#' @param target Which target to test.
#' @param check Whether to perform sanity checks.
#'
#' @returns data.frame of test results.
#' @export
#'
#' @examples
#' # TODO.
diff_breg_single <- function(data, design, target, check = TRUE) {
  # Sanity checks
  if (check){
    if (!(target %in% names(data)))
      stop("Target must be present in data")
    if (length(unique(data[[target]])) != 2)
      stop("Target must have two levels.")
    if (!is.factor(data[[target]])){
      data[[target]] <- as.factor(data[[target]])
      message("Converted target to factor.")
    }
    if (!isa(design, "formula"))
      stop("Design must be formula.")
    if (!(target %in% all.vars(design)))
      stop("Target must be present in design.")
  }

  # Run beta reg
  tryCatch({

    # Fit model with custom formula
    model <- glmmTMB::glmmTMB(
      formula = design,
      family = glmmTMB::beta_family(),
      data = data,
      ziformula = ~1,  # Zero-inflation component
    )

    # Calculate contrasts
    emm <-
      emmeans::emmeans(model, specs = target, type = "response") |>
      emmeans::regrid()
    contrast <- pairs(emm) |> as.data.frame()

    data.frame(
      ratio_diff = contrast$estimate[1],
      std_error = contrast$SE[1],
      p_value = contrast$p.value[1]
    )
  }, error = function(e) {
    data.frame(
      ratio_diff = NA_real_,
      std_error = NA_real_,
      p_value = NA_real_
    )
  })
}


#' Differential analysis with beta regression
#'
#' @param se SummarizedExperiment containing data to test.
#' @param design Design formula.
#' @param target Which target to perform test.
#' @param n.cores Multicore by `parallel`. Either "auto", integer, or NA.
#' NA means no parallel processing.
#'
#' @returns data.frame of test results with FDR.
#' @export
#'
#' @examples
#' # TODO
diff_breg <- function(se, design, target, n.cores = "auto") {

  # Extract relevant data from se
  assayName <- all.vars(design)[1]
  if (!(assayName %in% names(SummarizedExperiment::assays(se))))
    stop("Data not present in assays(se).")
  assay <- SummarizedExperiment::assays(se)[[assayName]]
  colData <- SummarizedExperiment::colData(se) |> as.data.frame()

  # Sanity check
  if (!is.factor(colData[[target]]) | length(levels(colData[[target]])) != 2)
    stop("Target must be factor with two levels.")
  if (!isa(design, "formula"))
    stop("Design must be formula.")
  if (!(target %in% all.vars(design)))
    stop("Target must be present in design.")

  # Setup parallel computing
  if (n.cores == "auto"){
    cl <- parallel::makeCluster(parallel::detectCores() - 1)

  } else if (is.integer(n.cores)){

  } else if (is.na(n.cores)){

  } else stop("n.cores must be 'auto', integer, or NA.")

  # Rowwise run analysis
  res_list <- list()
  for (i in 1:nrow(assay)){
    # Prepare data
    data <- data.frame(a=assay[i,])
    names(data) <- assayName
    data <- cbind(data, colData)

    # Run single
    res_list[[i]] <- diff_breg_single(data, design, target, check = FALSE)
  }

  # Combine results and compute adjusted p-value
  res <- Reduce(rbind, res_list)
  res$p_adj <- stats::p.adjust(res$p_value, method = "BH")
  res <- cbind(data.frame(se_rowname = row.names(se)), res)

  return(res)
}
