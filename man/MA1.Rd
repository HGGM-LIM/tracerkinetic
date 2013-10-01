\name{MA1}
\alias{MA1}
\title{Multilinear analysis 1 (MA1) for reversible 2-tissue model.}
\usage{
  MA1(input.function, tissue, time.start, time.end)
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
  This method needs to find the optimal linear point. This
  is yet to be implemented and this method should not be
  used. The code is nevertheless included in this file for
  completeness.
}
\references{
  M. Ichise, H. Toyama, R. Innis, and R. Carson,
  "Strategies to improve neuroreceptor parameter estimation
  by linear regression analysis," J. Cereb. Blood Flow
  Metab., vol. 22, pp. 1271-1281, 2002.
}

