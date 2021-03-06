\name{MA2}
\alias{MA2}
\title{Multilinear analysis 2 (MA2) for reversible 2-tissue model.}
\usage{
  MA2(input.function, tissue, time.start, time.end)
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
  kinetic parameters (Vt, Vs, Vn), and \code{fit}, which
  contains the linear least squares fit.
}
\description{
  Multilinear analysis 2 (MA2) for reversible 2-tissue
  model.
}
\references{
  M. Ichise, H. Toyama, R. Innis, and R. Carson,
  "Strategies to improve neuroreceptor parameter estimation
  by linear regression analysis," J. Cereb. Blood Flow
  Metab., vol. 22, pp. 1271-1281, 2002.

  Y. Ikoma, H. Watabe, M. Shidahara, M. Naganawa, and Y.
  Kimura, "PET kinetic analysis: error consideration of
  quantitative analysis in dynamic studies.," Ann. Nucl.
  Med., vol. 22, no. 1, pp. 1-11, Jan. 2008.
}

