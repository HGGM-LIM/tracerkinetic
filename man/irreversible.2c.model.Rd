\name{irreversible.2c.model}
\alias{irreversible.2c.model}
\title{Implements an irreversible two-tissue compartment model}
\usage{
  irreversible.2c.model(input.function, K1, k2, k3, vB,
    time.start, time.end, left.ventricle = input.function,
    interpolation.type = 1)
}
\arguments{
  \item{input.function}{Input function TAC.}

  \item{K1}{K1 parameter.}

  \item{k2}{k2 parameter.}

  \item{k3}{k3 parameter.}

  \item{vB}{Fraction of blood in tissue.}

  \item{time.start}{Initial acquisition time for each frame
  (in minutes).}

  \item{time.end}{Final acquisition time for each frame (in
  minutes).}

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
  This model is basically the same as
  \code{\link{reversible.2c.model}} but the parameter
  \code{k4} is always 0. The simplified equation has been
  implemented in this function and can be fitted using
  \code{\link{irreversible.2c.fit}}.
}

