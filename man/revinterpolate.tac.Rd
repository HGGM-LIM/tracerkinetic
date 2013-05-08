\name{revinterpolate.tac}
\alias{revinterpolate.tac}
\title{Actual interpolation.}
\usage{
  revinterpolate.tac(tac, tsample, time.points)
}
\arguments{
  \item{tac}{The extrapolated TAC to interpolate}

  \item{tsample}{The sampling points used in the original
  extrapolation.}

  \item{time.points}{The new time points for which the
  signal is desired.}
}
\value{
  The new interpolated TAC.
}
\description{
  Typically used after \code{\link{interpolate.tac}}. This
  is just a simplified call to the \code{\link{approx}}
  function.
}

