\name{patlak.plot}
\alias{patlak.plot}
\title{Patlak plot implementation}
\usage{
  patlak.plot(input.function, tissue, time.start, time.end,
    plot = TRUE, ...)
}
\arguments{
  \item{input.function}{Input function TAC.}

  \item{tissue}{Tissue TAC.}

  \item{time.start}{Initial acquisition time for each
  frame.}

  \item{time.end}{Final acquisition time for each frame.}

  \item{plot}{Should the Patlak plot be displayed? Defaults
  to \code{TRUE}.}

  \item{...}{Additional parameters passed to the
  \code{plot} function.}
}
\description{
  Implements a Patlak plot linear analysis, another
  graphical method for analysis of tracers that can be
  modeled after an irreversible two-tissue compartment
  model.
}
\references{
  TODO
}

