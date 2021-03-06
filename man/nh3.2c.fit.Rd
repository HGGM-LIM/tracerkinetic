\name{nh3.2c.fit}
\alias{nh3.2c.fit}
\title{Computes the two-tissue compartmental model parameters for the 13-ammonia
tracer.}
\usage{
  nh3.2c.fit(input.function, tissue, time.start, time.end,
    K1.start = 0.1, K1.lower = 0, K1.upper = 8,
    k2.start = 0.1, k2.lower = 0, k2.upper = 8,
    k3.start = 0.1, k3.lower = 0, k3.upper = 8,
    vB.start = 0.05, vB.lower = 0, vB.upper = 1,
    weight = NA, plot = FALSE, interpolation.type = 1, ...)
}
\arguments{
  \item{input.function}{Input function TAC.}

  \item{tissue}{Tissue TAC.}

  \item{time.start}{Initial acquisition time for each frame
  (in minutes).}

  \item{time.end}{Final acquisition time for each frame (in
  minutes).}

  \item{K1.start,k2.start,k3.start,vB.start}{Initial
  parameter values.}

  \item{K1.lower,k2.lower,k3.lower,vB.lower}{Parameter
  lower bounds.}

  \item{K1.upper,k2.upper,k3.upper,vB.upper}{Parameter
  upper bounds.}

  \item{weight}{Weights for the non-linear least squares.
  Defaults to no weighting. Use \code{"framelength"} to use
  the frame length as the weight value.}

  \item{plot}{\code{TRUE} if a plot is to be shown.
  Defaults to \code{FALSE}.}

  \item{interpolation.type}{Interpolation type selection.
  Passed to \code{\link{interpolate.tac}}. Defaults to 1.}

  \item{...}{Other parameters passed to the plot function,
  if used.}
}
\value{
  Returns a list with four fields: \code{kparms}, the
  computed kinetic parameters; \code{stderrors}, the
  standard errors for each parameter; \code{stderrorsp},
  the standard errors for each parameter as a percentage;
  \code{fit}, the actual fitted object.
}
\description{
  Automatically performs the metabolite correction
  implemented in \code{\link{nh3.vdhoff.correction}} to the
  \code{input.function} term. The \code{left.ventricle}
  term is not provided. All the other parameters are
  provided just like in the
  \code{\link{irreversible.2c.fit}} function and the return
  values are the same.
}

