\name{flow2k1}
\alias{flow2k1}
\title{K1 to flow conversion}
\usage{
  flow2k1(flow, a = 0.77, b = 0.63)
}
\arguments{
  \item{flow}{The myocardial blood flow}

  \item{a}{Parameter \code{a} of the Renkin-Crone equation.
  Defaults to 0.77 [unitless].}

  \item{b}{Parameter \code{b} of the Renkin-Crone equation.
  Defaults to 0.63 [ml/min/g].}
}
\value{
  The K1 parameter value.
}
\description{
  Implements the Renkin-Crone equation that transforms flow
  values to K1 in the 82-rubidium tracer.
}

