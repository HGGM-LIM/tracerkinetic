\name{processlinear}
\alias{processlinear}
\title{Helper function for linear models}
\usage{
  processlinear(termx, termy, plot = TRUE, ...)
}
\arguments{
  \item{termx,termy}{The x and y values of the graph.}

  \item{plot}{Should the result be shown? Defaults to
  \code{TRUE}.}

  \item{...}{Further parameters for the \code{plot}
  function.}
}
\value{
  Returns a list with three fields: \code{kparms}, the
  computed kinetic parameters; \code{stderrors}, the
  standard errors for each parameter as a percentage;
  \code{fit}, the actual fitted object.
}
\description{
  This function helps process linear models, such as
  \code{\link{patlak.plot}} or \code{\link{logan.plot}}. It
  finds the optimal fitting point by calling
  \code{\link{findbestfit}} and proceeeds to compute the
  kinetic parameters, that will mostly consist of the
  intercept and slope of the linear regression, along with
  their standard errors.
}
\details{
  It \code{plot} is set to \code{TRUE}, the regression will
  be shown on screen. The data points used to compute the
  regression are displayed with a filled style, while the
  ones that have not been used are displayed using open
  circles.
}

