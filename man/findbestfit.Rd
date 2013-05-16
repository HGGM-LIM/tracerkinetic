\name{findbestfit}
\alias{findbestfit}
\title{Finds optimal fitting point}
\usage{
  findbestfit(x, y, minpoints = 3, maxerror = 0.1)
}
\arguments{
  \item{x,y}{Regression data}

  \item{minpoints}{Minimum data points to be used for the
  regression.}

  \item{maxerror}{The maximum error rate allowed.}
}
\value{
  A list containing the time point chosen,
  \code{initial.time.point}, and the actual fitted object,
  \code{fit}.
}
\description{
  Finds the number of points that have to be used for an
  optimal fit between the \code{x} and \code{y} variables.
}
\details{
  This function finds the earliest point that allows to
  have a regression with less than 5% error between the
  chosen points (this parameter is controlled by the
  \code{maxerror} variable). If no point allows to have
  this error rate, the point that yields the minimum
  regression error is used.
}

