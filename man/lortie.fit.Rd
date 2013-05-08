\name{lortie.fit}
\alias{lortie.fit}
\title{Lortie et al. (2007) model fit.}
\usage{
  lortie.fit(input.function, tissue, time.start, time.end,
    a = 0.77, b = 0.63, FLV.start = 0.1, FLV.lower = 0,
    FLV.upper = 1, flow.start = 0.5, flow.lower = 0,
    flow.upper = 8, k2.start = 0.1, k2.lower = 0,
    k2.upper = 8, left.ventricle = input.function,
    weight = time.end - time.start, plot = FALSE,
    interpolation.type = 1, ...)
}
\arguments{
  \item{input.function}{Input function TAC (typically, left
  ventricle).}

  \item{tissue}{Tissue (myocardium) TAC.}

  \item{time.start}{Initial time for each frame.}

  \item{time.end}{End time for each frame.}

  \item{a}{Parameter \code{a} of the Renkin-Crone equation
  implemented in \code{\link{flow2k1}}. Defaults to 0.77
  [unitless].}

  \item{b}{Parameter \code{b} of the Renkin-Crone equation
  implemented in \code{\link{flow2k1}}. Defaults to 0.63
  [ml/min/g].}

  \item{FLV.start,flow.start,k2.start}{Initial parameter
  values for the non-linear squares fitting.}

  \item{FLV.lower,flow.lower,k2.lower}{Lower bounds.}

  \item{FLV.upper,flow.upper,k2.upper}{Upper bounds.}

  \item{left.ventricle}{Left ventricle TAC, used for
  spill-over corretion.  Defaults to
  \code{input.function}.}

  \item{weight}{Weights for the non-linear least squares.
  Defaults to frame length.}

  \item{plot}{Should the results be plotted. Defaults to
  \code{FALSE}.}

  \item{interpolation.type}{Interpolation type selection.
  Passed to \code{\link{interpolate.tac}}. Defaults to 1.}

  \item{...}{If plotting is used, optional parameters are
  passed to the plot function.}
}
\value{
  Returns a list with three fields: \code{kparms}, the
  computed kinetic parameters; \code{stderrors}, the
  standard errors for each parameter as a percentage;
  \code{fit}, the actual fitted object.
}
\description{
  Fits the one-tissue compatment 82-rubidium model
  implemented in \code{\link{lortie.model}}.
}

