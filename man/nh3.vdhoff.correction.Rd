\name{nh3.vdhoff.correction}
\alias{nh3.vdhoff.correction}
\title{NH3 metabolite correction}
\usage{
  nh3.vdhoff.correction(input.function, time.end,
    t0 = 0.48, T12 = 6.69)
}
\arguments{
  \item{input.function}{The input function TAC.}

  \item{time.end}{Final acquisition time for each frame.}

  \item{t0}{Empirical value for the t0 correction
  parameter. Defaults to 0.48 [min].}

  \item{T12}{Empirical value for the T1/2 correction
  parameter. Defaults to 6.69 [min].}
}
\description{
  Implements the 13-NH3 metabolite correction as published
  in [1]. This correction should be applied to the input
  function before using it as the input for any model.
}
\references{
  [1] J. van den Hoff et al., "[1-11C]Acetate as a
  Quantitative Perfusion Tracer in Myocardial PET," Journal
  of Nuclear Medicine, vol. 42, no. 8, pp. 1174-1182, Aug.
  2001.
}

