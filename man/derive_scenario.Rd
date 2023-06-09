% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctsframework.R
\name{derive_scenario}
\alias{derive_scenario}
\title{Derive a new scenario from an existing one}
\usage{
derive_scenario(
  base_scenario,
  parameters = NULL,
  population = NULL,
  subject_model = NULL,
  treatment = NULL,
  pkpd_model = NULL,
  simulation = NULL,
  endpoints = NULL,
  summarization = NULL,
  ...,
  name,
  desc = NULL,
  seed = NULL
)
}
\arguments{
\item{base_scenario}{A \link{scenario} which serves as the base scenario from which to derive.}

\item{parameters}{\link{parameters} module (optional).}

\item{population}{\link{population} module (optional).}

\item{subject_model}{\link{subject_model} module (optional).}

\item{treatment}{\link{treatment} module (optional).}

\item{pkpd_model}{\link{pkpd_model} module (optional).}

\item{simulation}{\link{simulation} module (optional).}

\item{endpoints}{\link{endpoints} module (optional).}

\item{summarization}{\link{summarization} module (optional).}

\item{...}{Further argument passed to \code{\link[=scenario]{scenario()}}.}

\item{name}{Optional scenario name.}

\item{desc}{Optional scenario description.}

\item{seed}{RNG seed, for reproducibility.}
}
\description{
A new scenario can be derived from an existing one by overriding any of its
modules, or specific components of its modules.
}
