% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctsframework.R
\name{cts_module}
\alias{cts_module}
\alias{parameters}
\alias{population}
\alias{subject_model}
\alias{treatment}
\alias{pkpd_model}
\alias{simulation}
\alias{endpoints}
\alias{summarization}
\title{CTS modules}
\usage{
parameters(
  ...,
  outputs = NULL,
  name = NULL,
  desc = NULL,
  seed = NULL,
  collate = NULL
)

population(
  ...,
  outputs = NULL,
  name = NULL,
  desc = NULL,
  seed = NULL,
  collate = NULL
)

subject_model(
  ...,
  outputs = NULL,
  name = NULL,
  desc = NULL,
  seed = NULL,
  collate = NULL
)

treatment(
  ...,
  outputs = NULL,
  name = NULL,
  desc = NULL,
  seed = NULL,
  collate = NULL
)

pkpd_model(
  ...,
  outputs = NULL,
  name = NULL,
  desc = NULL,
  seed = NULL,
  collate = NULL
)

simulation(
  ...,
  outputs = NULL,
  name = NULL,
  desc = NULL,
  seed = NULL,
  collate = NULL
)

endpoints(
  ...,
  outputs = NULL,
  name = NULL,
  desc = NULL,
  seed = NULL,
  collate = NULL
)

summarization(
  ...,
  outputs = NULL,
  name = NULL,
  desc = NULL,
  seed = NULL,
  collate = NULL
)
}
\arguments{
\item{...}{A named list of expressions.}

\item{outputs}{(optional) A \link{cts_output} or \link{cts_output_list}.}

\item{name}{Modules can be given a name.}

\item{seed}{Modules can be given a RNG seed, for reproducibility.}
}
\value{
A module object.
}
\description{
Modules are the components that make up a CTS scenario. The power of the
framework comes from the ability to mix and match modules to easily create a
multitute of scenarios.
}
\details{
Each module has the same structure: it consists of a named list of
expressions that are run in sequence when the scenario is executed.
}
