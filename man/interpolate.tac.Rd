\name{interpolate.tac}
\alias{interpolate.tac}
\title{Interpolates the given TAC to perform the correlation.}
\usage{
  interpolate.tac(tac, time.start, time.end,
    interpolation.type = 1)
}
\arguments{
  \item{tac}{The TAC to be interpolated.}

  \item{time.start}{Initial acquisition time for each
  frame.}

  \item{time.end}{Final acquisition time for each frame.}

  \item{interpolation.type}{Interpolation type. 1 for
  linear, 2 for "natural" splines. Defaults to 1.}
}
\value{
  A list containing \code{tsample}, the sampling time
  vector and the \code{y}, the interpolated TAC.
}
\description{
  This interpolation process (actually, an extrapolation)
  is needed to convert the given TAC into a signal with
  frames of equal length. This length is the minimum frame
  length for the given study.
}

