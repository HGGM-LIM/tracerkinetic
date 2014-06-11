\name{lortie.model}
\alias{lortie.model}
\title{Reversible one-tissue compartment 82-rubidium model with flow correction}
\usage{
  lortie.model(input.function, FLV, flow, k2, time.start,
    time.end, a = 0.77, b = 0.63,
    left.ventricle = input.function,
    interpolation.type = 1)
}
\arguments{
  \item{input.function}{Input function TAC.}

  \item{FLV}{Fraction of blood (from left ventricle) in
  myocardium TAC.}

  \item{flow}{Myocardial blood flow.}

  \item{k2}{Rate constant.}

  \item{time.start}{Initial acquisition time for each frame
  (in minutes).}

  \item{time.end}{Final acquisition time for each frame (in
  minutes).}

  \item{a}{Parameter \code{a} of the Renkin-Crone equation
  implemented in \code{\link{flow2k1}}. Defaults to 0.77
  [unitless].}

  \item{b}{Parameter \code{b} of the Renkin-Crone equation
  implemented in \code{\link{flow2k1}}. Defaults to 0.63
  [ml/min/g].}

  \item{left.ventricle}{Left ventricle TAC, used for
  spill-over corretion.  Defaults to
  \code{input.function}.}

  \item{interpolation.type}{Interpolation type selection.
  Passed to \code{\link{interpolate.tac}}. Defaults to 1.}
}
\value{
  The TAC resulting of solving the model with the given
  parameters.
}
\description{
  This function implements the one-tissue compartment
  rubidium-82 model published in [1].
}
\references{
  [1] M. Lortie, R. S. B. Beanlands, K. Yoshinaga, R.
  Klein, J. N.  Dasilva, and R. a DeKemp, "Quantification
  of myocardial blood flow with 82Rb dynamic PET imaging.,"
  European journal of nuclear medicine and molecular
  imaging, vol. 34, no. 11, pp. 1765-74, Nov. 2007.
}

