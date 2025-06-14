% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctsframework.R
\name{collate_csv_outputs}
\alias{collate_csv_outputs}
\title{Collate simulation results from multiple scenarios}
\usage{
collate_csv_outputs(
  module,
  name,
  ...,
  outputs = NULL,
  .path = file.path("cts_results", as.Date(Sys.time()))
)
}
\arguments{
\item{module}{A \code{character}. What is the module that generated the output?}

\item{name}{A \code{character}. The name given to the output in the \link{scenario}.}

\item{...}{Additional arguments (currently ignored).}

\item{outputs}{(optional) A \link{cts_output} or \link{cts_output_list}.}

\item{.path}{The path containing the simulation results. This directory
should contain the file \code{scenarios.csv} as well as subfolders for each
scenario, as generated by \code{\link[=CTS]{CTS()}}.}
}
\value{
A \link[data.table:data.table]{data.table::data.table}.
}
\description{
This function can be used after running multiple scenarios with \code{CTS} to
combine the results for further processing and output generation.
}
