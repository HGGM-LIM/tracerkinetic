\name{patlak.plot}
\alias{patlak.plot}
\title{Patlak plot implementation}
\usage{
  patlak.plot(input.function, tissue, time.start, time.end,
    plot = TRUE, maxerror = 0.1, ...)
}
\arguments{
  \item{input.function}{Input function TAC.}

  \item{tissue}{Tissue TAC.}

  \item{time.start}{Initial acquisition time for each
  frame.}

  \item{time.end}{Final acquisition time for each frame.}

  \item{plot}{Should the Patlak plot be displayed? Defaults
  to \code{TRUE}.}

  \item{maxerror}{\code{maxerror} parameter to be passed to
  \code{processlinear}.}

  \item{...}{Additional parameters passed to the
  \code{plot} function.}
}
\value{
  Returns a list with three fields: \code{kparms}, the
  computed kinetic parameters; \code{stderrors}, the
  standard errors for each parameter as a percentage;
  \code{fit}, the actual fitted object.
}
\description{
  Implements a Patlak plot linear analysis, another
  graphical method for analysis of tracers that can be
  modeled after an irreversible two-tissue compartment
  model.
}
\details{
  Please refer to the documentation of
  \code{\link{processlinear}} for implementation details.
}
\references{
  C. Patlak, R. Blasberg, and J. Fenstermacher, "Graphical
  evaluation of blood-to-brain transfer constants from
  multiple-time uptake data," J Cereb Blood Flow Metab,
  1983.
}
\seealso{
  \code{\link{logan.plot}}.
}

