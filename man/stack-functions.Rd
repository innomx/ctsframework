% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctsframework.R
\name{stack-functions}
\alias{stack-functions}
\alias{push}
\alias{pop}
\alias{top}
\title{Stack Functions}
\usage{
push(stack, value)

pop(stack)

top(stack)
}
\arguments{
\item{stack}{A list.}

\item{value}{Any value or object.}
}
\value{
\code{push} and \code{pop} return the stack in its state after the specified
operation. \code{top} return the value currently at the top of the stack, or null
if the stack is empty.
}
\description{
A simple stack implementation. Provides the basic functions \code{push}, \code{pop}
and \code{top}. The stack itself is implemented as a simple list.
}
\details{
In the event the stack is empty, both \code{pop} and \code{top} issue a warning.
}
\examples{
stack <- list()
stack <- stack |> push(3) |> push(5) |> push(7)
top(stack)
stack <- pop(stack)
top(stack)
stack <- pop(stack)
top(stack)
stack <- pop(stack)
top(stack)
stack <- pop(stack)
}
