\name{murase}
\alias{murase}
\title{Multilinear analysis for 1-compartment model.}
\usage{
  murase(input.function, tissue, time.start, time.end)
}
\arguments{
  \item{input.function}{Input function TAC.}

  \item{tissue}{Tissue TAC.}

  \item{time.start}{Initial acquisition time for each frame
  (in minutes).}

  \item{time.end}{Final acquisition time for each frame (in
  minutes).}
}
\value{
  A list with two fields: \code{kparms}, that contains the
  kinetic parameters, and \code{fit}, which contains the
  linear least squares fit.
}
\description{
  Multilinear analysis for 1-compartment model.
}
\references{
  K. Murase, "Efficient method for calculating kinetic
  parameters using T1-weighted dynamic contrast-enhanced
  magnetic resonance imaging.," Magn. Reson. Med., vol. 51,
  no. 4, pp. 858-62, Apr. 2004.
}

