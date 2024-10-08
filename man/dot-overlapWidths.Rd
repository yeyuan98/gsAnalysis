% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intron.properties.R
\name{.overlapWidths}
\alias{.overlapWidths}
\title{Count number of overlapping bases}
\usage{
.overlapWidths(query, subject)
}
\arguments{
\item{query}{GenomicRanges::Ganges of query.}

\item{subject}{GenomicRanges::Ganges of subject.}
}
\value{
Integer vector of numbers of overlapping bases. Guaranteed to be
the same order as the query ranges.
}
\description{
Count number of overlapping bases
}
\details{
'overlapping bases' is counted for each query against ALL subject ranges.
If a query overlaps with two subject ranges with 5 and 4 bases, the number
reported will be 5+4=9. Implementation is:
\enumerate{
\item \code{GenomicRanges::subtract(query, subject)}
\item Count up the width of the resulting GRangesList
\item (original width of query) - (width of subtracted query)
}
}
