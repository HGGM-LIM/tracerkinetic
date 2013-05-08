\name{reversible.1c.model}
\alias{reversible.1c.model}
\title{Implements a reversible one-tissue compartment model}
\usage{
  reversible.1c.model(input.function, K1, k2, vB,
    time.start, time.end, left.ventricle = input.function,
    interpolation.type = 1)
}
\arguments{
  \item{input.function}{Input function TAC.}

  \item{K1}{K1 parameter.}

  \item{k2}{k2 parameter.}

  \item{vB}{Fraction of blood in tissue.}

  \item{time.start}{Initial acquisition time for each
  frame.}

  \item{time.end}{Final acquisition time for each frame.}

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
  Implements a reversible one-tissue compartment model
}

