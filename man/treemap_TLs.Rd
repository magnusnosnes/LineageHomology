% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lineage_plots.R
\name{treemap_TLs}
\alias{treemap_TLs}
\title{treemap_TLs}
\usage{
treemap_TLs(Result_LineageHomology, date_gradient = NULL)
}
\arguments{
\item{Result_LineageHomology}{The output from LineageHomology.}

\item{date_gradient}{A user-specified range for the TMRCA values to be plotted. If not provided, the range will be calculated from the data.}
}
\value{
A treemap visualization of the results in the TL assignment from LineageHomology, showing the group number, group size, and TMRCA for each lineage.
}
\description{
treemap_TLs creates a treemap visualization of the Lineage Homology results, showing the group number, group size, and time-to-most-recent-common-ancestor (TMRCA) for each lineage.
The treemap is color-coded based on the TMRCA values, with a user-specified range or calculated range from the data.
}
