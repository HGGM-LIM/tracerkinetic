\name{MLAIR2}
\alias{MLAIR2}
\title{Implements the Multiple linear analysis for irreversible radiotracer 2
(MLAIR2)}
\usage{
  MLAIR2(input.function, tissue, time.start, time.end)
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
  Implements the Multiple linear analysis for irreversible
  radiotracer 2 (MLAIR2)
}
\references{
  S. J. Kim, J. S. Lee, Y. K. Kim, J. Frost, G. Wand, M. E.
  McCaul, and D. S. Lee, “Multiple linear analysis
  methods for the quantification of irreversibly binding
  radiotracers.,” J. Cereb. Blood Flow Metab., vol. 28,
  no. 12, pp. 1965–77, Dec. 2008.
}

