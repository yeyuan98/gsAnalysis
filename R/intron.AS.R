# Summary functions of alternative splicing analysis routines
#   rMATS / rMATS-TURBO
#   IRFinder-S

#' Read rMATS output of different AS patterns
#'
#' @param outputs.dir Output directory of rMATS.
#' @param method Either "JC" or "JCEC", see details.
#'
#' @return Data frame of rMATS output
#' @export
#'
#' @details
#' This is a convenience function for reading rMATS output, merging different
#' splicing patterns into a single data frame for easier further analysis.
#' Refer to rMATS Github for details on results generated by rMATS:
#' https://github.com/Xinglab/rmats-turbo/blob/v4.3.0/README.md#output
#'
#' rMATS normalizes lengths of individual splicing variants. There are two
#' methods used for this normalization: JC and JCEC. Refer to the rMATS paper
#' for more details: https://doi.org/10.1073/pnas.1419161111
#'
#' rMATS coordinated are 0-based; the exact meaning of start and end varies by
#' splicing pattern type:
#'
#' - A3SS/A5SS: start/end = start/end of the long exon (inclusion form).
#' - SE: start/end = start/end of the skipped exon (inclusion form).
#' - MXE: start/end = start of the first exon / end of the second exon.
#' - RI: start/end = start/end of the retained intron (inclusion form).
#'
#' @examples
#' # TODO
rmats_read <- function(outputs.dir, method){

  # method is either JC or JCEC
  if (method != "JC" & method != "JCEC")
    stop("method must be either JC or JCEC")

  # AS patterns analyzed by rMATS
  patterns <- c("A3SS", "A5SS", "SE", "MXE", "RI")

  # Columns to get (shared by all patterns)
  columns <- c(
    "geneSymbol", "chr", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference",
    "IJC_SAMPLE_1", "SJC_SAMPLE_1",
    "IJC_SAMPLE_2", "SJC_SAMPLE_2"
    # start and end - definition depends on the specific pattern
    # ... not defined here
    # type - the specific AS pattern
    # ... not defined here
  )

  # Pattern-specific mapping of start and end
  mStart <- list(
    A3SS = "longExonStart_0base",
    A5SS = "longExonStart_0base",
    SE = "exonStart_0base",
    MXE = "1stExonStart_0base",
    RI = "upstreamEE"
  )
  mEnd <- list(
    A3SS = "longExonEnd",
    A5SS = "longExonEnd",
    SE = "exonEnd",
    MXE = "2ndExonEnd",
    RI = "downstreamES"
  )

  # Read and merge
  data <- lapply(patterns, function(p){

    # Read MATS file
    suppressMessages(
      df <- readr::read_delim(
        file.path(outputs.dir, paste0(p, ".MATS.", method, ".txt")),
        delim = "\t",
        show_col_types = FALSE
      )
    )

    # Extract columns of interest
    coi <- c(columns, mStart[[p]], mEnd[[p]])
    df <- df[,coi]
    df <- dplyr::rename(
      df,
      start_0base = mStart[[p]],
      end = mEnd[[p]]
    )
    df$Type <- p

    df
  })

  return(Reduce(rbind, data))
}


#' Assign filter labels of rMATS output
#'
#' @param df rMATS data read by `rmats_read()`.
#' @param readCov Minimum read coverage.
#' @param maxPSI Maximum percent spliced-in PSI.
#' @param minPSI Minimum PSI.
#' @param sigFDR FDR threshold for significance.
#' @param sigDeltaPSI PSI threshold for significance.
#' @param res_column Which column to store filter results
#'
#' @returns rMATS data table with filter label column `filter`.
#' @export
#'
#' @details
#' This function reproduces the data filter demonstrated in Xing Lab 2024
#' Nature Protocol tutorial. See Github repository
#' Xinglab/rmats-turbo-tutorial and refer to rmats_filtering.py.
#'
#' AS events that does not satisfy any of the following criteria will be given
#' the *remove* label:
#'
#' 1. Average read counts (IJC+SJC) in both Sample 1 and Sample 2 conditions
#' are at least `readCov`.
#' 2. Average PSI in both Sample 1 and Sample 2 conditions are between
#' `minPSI` and `maxPSI`.
#'
#' If all criteria above are satisfied, events will be labeled by statistical
#' significance:
#'
#' 1. `FDR < sigFDR` and `IncLevelDifference > sigDeltaPSI`: *up*.
#' 2. `FDR < sigFDR` and `IncLevelDifference < -sigDeltaPSI`: *down*.
#' 3. Otherwise: *ns*.
#'
#' @examples
#' # TODO
rmats_filter <- function(
    df, readCov = 10, maxPSI = .95, minPSI = .05,
    sigFDR = .1, sigDeltaPSI = .05, res_column = "filter"
){

  # Calculate parameters
  suppressWarnings({
    averageCountSample1 <-
      y3628::str_split_summary(df$IJC_SAMPLE_1, ",", \(x) sum(as.numeric(x)))+
      y3628::str_split_summary(df$SJC_SAMPLE_1, ",", \(x) sum(as.numeric(x)))
    averageCountSample1 <-
      averageCountSample1 /
      y3628::str_split_summary(df$IJC_SAMPLE_1, ",", length)

    averageCountSample2 <-
      y3628::str_split_summary(df$IJC_SAMPLE_2, ",", \(x) sum(as.numeric(x)))+
      y3628::str_split_summary(df$SJC_SAMPLE_2, ",", \(x) sum(as.numeric(x)))
    averageCountSample2 <-
      averageCountSample2 /
      y3628::str_split_summary(df$IJC_SAMPLE_2, ",", length)

    averagePsiSample1 <-
      y3628::str_split_summary(df$IncLevel1, ",", \(x) sum(as.numeric(x), na.rm = TRUE))
    averagePsiSample1 <-
      averagePsiSample1 /
      y3628::str_split_summary(df$IncLevel1, ",", length)

    averagePsiSample2 <-
      y3628::str_split_summary(df$IncLevel2, ",", \(x) sum(as.numeric(x), na.rm = TRUE))
    averagePsiSample2 <-
      averagePsiSample2 /
      y3628::str_split_summary(df$IncLevel2, ",", length)
  })

  # Designate group
  #   Here, following the XingLab standards.
  df[[res_column]] <- ifelse(
    averageCountSample1 >= readCov &
      averageCountSample2 >= readCov &
      averagePsiSample1 <= maxPSI & averagePsiSample2 <= maxPSI &
      averagePsiSample1 >= minPSI & averagePsiSample2 >= minPSI,
    ifelse(
      df$FDR <= sigFDR & df$IncLevelDifference >= sigDeltaPSI,
      "up", ifelse(
        df$FDR <= sigFDR & df$IncLevelDifference <= -sigDeltaPSI,
        "down", "ns"
      )
    ), "remove"
  )
  return(df)
}

#' Converts rMATS data frame to GenomicRange
#'
#' Convenience function for converting rMATS data to GRanges.
#'
#' @param df rMATS data read by `rmats_read()`.
#' @param ... <[`tidy-select`][dplyr::dplyr_tidy_select]> See details.
#'
#' Three columns are used to fill in GRanges information: `chr`, `start_0base`, `end`.
#'
#' Two columns are by default added to the metadata: `geneSymbol`, `Type`.
#'
#' Extra metadata columns are specified by the tidyselect ellipsis.
#'
#'
#' @return GRanges object
#' @export
#'
#' @examples
#' #TODO
rmats_toGRange <- function(df, ...){
  gr <- GenomicRanges::makeGRangesFromDataFrame(
    df, start.field = "start_0base", starts.in.df.are.0based = TRUE,
    seqnames.field = "chr", ignore.strand = TRUE
  )
  meta <- dplyr::select(df, c(dplyr::matches("geneSymbol|Type"), ...))
  GenomicRanges::mcols(gr) <- meta
  return(gr)
}


#' Read IRFinder-S output of a single sample
#'
#' @param result.dir IRFinder-S output directory
#' @param type.samples Which summary type to read in, see details
#'
#' @return tibble of IRFinder-S results
#' @export
#'
#' @details
#' Two sample types are output by IRFinder-S: `validated` and `full`.
#'
#' - `validated`: read the "IRFinder-IR-nondir-val.txt"
#' - `full`: read the "IRFinder-IR-nondir.txt"
#'
#' For column names please refer to the function body
#' @examples
#' #TODO
IRFinderS_read <- function(result.dir, type.samples = "validated"){

  # Report WARNINGS
  war <- readr::read_file(file.path(result.dir, "WARNINGS"))
  if ( stringr::str_length(war) == 0 ){
    message("No IRFinder warnings for sample = ", result.dir)
  } else{
    warning("IRFinder reported warnings for sample = ", result.dir, war, "\n")
  }

  # Read in the result table
  summary.colnames.basic <-
    c(
      "chr", "start", "end", "name", "not.used", "strand",
      "n.excludedBases", "coverage", "depth", "depth.25perc",
      "depth.50perc", "depth.75perc", "n.reads.exon2intronLeft",
      "n.reads.exon2intronRight", "depth.first50bp", "depth.last50bp",
      "n.reads.spliceLeft", "n.reads.spliceRight", "n.reads.spliceExact",
      "IR.ratio", "warning"
    )
  if (type.samples == "validated"){
    message("Reading entries from validated table")
    summary.table <- "IRFinder-IR-nondir-val.txt"
    summary.colnames <- summary.colnames.basic
    summary.colnames[5] <- "cnn.score"
  } else if (type.samples == "full"){
    message("Reading entries from all results table")
    summary.table <- "IRFinder-IR-nondir.txt"
    summary.colnames <- summary.colnames.basic
  } else{
    stop("Unsupported summary type.")
  }

  readr::read_table(
    file.path(result.dir, summary.table),
    col_names = summary.colnames, skip = 1,
    col_types = "ciicdcidddddiiddiiidc"
  ) |>
    dplyr::mutate(
      # Extract the gene symbol which is useful as identifier
      symbol = stringr::str_split(name, "/", simplify = TRUE)[,1]
    ) |>
    dplyr::mutate(
      start = start + 1 # IRFinder is ZERO-based for start!
    )
}


#' Read IRFinder-S output of multiple samples
#'
#' Report a "merged" SummarizedExperiment for differential analysis. Refer to
#' details section.
#'
#' The routine is divided into the following steps:
#'
#' First, read in retained introns from each sample by `IRFinderS_read`.
#'
#' `named.result.dirs` must be a named character vector whose:
#'
#' * names = group name
#' * values = full path to the sample result directories
#'
#' Introns that have warnings defined by `wl` is removed from each sample.
#' `wl` definition follows that of the `IRFinder Diff` routine.
#'
#' Unique introns are defined by columns specified in `join.columns`.
#'
#' Next, consensus introns are determined. If an intron is found in at least
#' `min.samples` of samples it is considered as a consensus intron.
#'
#' Finally, intron data from the 'full' data table are extracted,
#' which will be put as assays in the output SE object. Which data columns are
#' included is defined by `assay.columns`.
#'
#' ## `wl` filtering
#'
#' This filtering determines whether an intron is included in a sample.
#' It looks at the warning flags by IRFinder and possible values are:
#'
#' 0 -> no filter; 1 -> LowCover introns removed; 2 -> Lowcover & LowSplicing;
#' 3 -> LowCover & LowSplicing & MinorIsoform;
#' 4 -> LowCover & LowSplicing & MinorIsoform & NonUniformIntronCover
#'
#' @param named.result.dirs Named paths to result directories, see description.
#' @param assay.columns Which column(s) to read in as assays of the output.
#' @param join.columns Which columns define an intron.
#' @param type.samples Forwarded to `IRFinderS_read()`.
#' @param min.samples How to filter retained introns, see description.
#' @param wl How to filter retained introns, see description.
#'
#' @return `SummarizedExperiment` object of IRFinder results.
#' @export
#'
#' @examples
#' #TODO
IRFinderS_readSamples <- function(
    named.result.dirs,
    assay.columns = c("IR.ratio", "depth", "n.reads.spliceExact"),
    join.columns = c("chr", "start", "end", "strand", "symbol"),
    type.samples = "validated", min.samples = 3, wl = 0
){

  # Check for names
  if (is.null(names(named.result.dirs)))
    stop("Must provide sample group names as names.")
  message("Group names provided are: ",
          paste(unique(names(named.result.dirs)), collapse = ", "))

  # wl filter follows the IRFinder Diff definition
  if (wl %in% 0:4){
    wl <- switch(
      wl+1, "", "LowCover", c("LowCover", "LowSplicing"),
      c("LowCover", "LowSplicing", "MinorIsoform"),
      c("LowCover", "LowSplicing", "MinorIsoform", "NonUniformIntronCover"))
  } else stop("Invalid warning filter setting.")

  # Determine the introns that pass warning level filter
  if (type.samples %in% c("validated", "full")){
    # Read in samples
    dfs.val <- lapply(
      named.result.dirs, function(dir.path){
        suppressMessages(
          IRFinderS_read(dir.path, type.samples = type.samples)) |>
          dplyr::filter(!(warning %in% wl))
      })
    #   only keep the join columns
    dfs.val <- lapply(
      dfs.val, function(df){
        dplyr::select(df, dplyr::all_of(join.columns))
      }
    )
  } else stop("Unsupported sample type.")

  #   derive the 'consensus' join column table
  #     if a join combination exists >= min.samples samples, keep
  consensusFromVal <- lapply(
    dfs.val, function(df)
      with(df, paste(chr, start, end, strand, symbol, sep = "@@"))
  )
  consensusFromVal <- unlist(consensusFromVal)
  consensusFromVal <- table(consensusFromVal)
  consensusFromVal <-
    names(consensusFromVal)[consensusFromVal >= min.samples]
  consensusFromVal <-
    stringr::str_split(consensusFromVal, "@@", simplify = TRUE)
  consensusFromVal <- tibble::tibble(
    chr = consensusFromVal[,1], start = as.integer(consensusFromVal[,2]),
    end = as.integer(consensusFromVal[,3]), strand = consensusFromVal[,4],
    symbol = consensusFromVal[,5]
  )

  # Read in samples - full data tables
  dfs <- lapply(
    named.result.dirs, function(dir.path)
      suppressMessages(IRFinderS_read(dir.path, type.samples = "full")))
  #   Only grab the ones in consensus
  dfs <- lapply(
    dfs, function(df)
      dplyr::left_join(
        # preserve order of consensusFromVal
        consensusFromVal, df,
        by = c("chr", "start", "end", "strand", "symbol"))
  )

  # Construct SummarizedExperiment
  #   rowRanges
  gr <- GenomicRanges::makeGRangesFromDataFrame(consensusFromVal)
  names(gr) <- with(consensusFromVal, paste0(symbol,"_",start,"-",end))
  #   assay matrices
  assay.mats <- lapply(
    assay.columns, \(column){
      mat <- lapply(dfs, function(df) df[[column]])
      Reduce(cbind, mat)
    }
  )
  names(assay.mats) <- assay.columns
  #   colData
  coldat <- data.frame(group = names(named.result.dirs))
  #   construct object
  res <- SummarizedExperiment::SummarizedExperiment(
    assays = assay.mats, rowRanges = gr, colData = coldat
  )
  colnames(res) <- NULL
  return(res)
}
