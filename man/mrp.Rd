% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrp.R
\name{mrp}
\alias{mrp}
\title{Max-Ratio Projection function}
\usage{
mrp(data, cl, npc = 4, message = T, ...)
}
\arguments{
\item{data}{The dataset to apply MRP. Each row is an observation.}

\item{cl}{The class identification for each observation. The length of \code{cl} should be the same as the number of rows of \code{data}.}

\item{npc}{The number of max-ratio directions to be used in projecting the original data to the projected data.}

\item{message}{Logical. Wheather to show the accumulative variance explained by the projection directions or not.}
}
\value{
A list with the elements
\item{projected_df}{The projected data with selected number of max-ratio directions.}
\item{pccumvar}{The cummulative variance explained by the max-ratio principal components.}
}
\description{
function to project high-dimensional datasets to lower dimention with max-ratio projection.
}
