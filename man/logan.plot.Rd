\name{logan.plot}
\alias{logan.plot}
\title{Logan plot implementation}
\usage{
  logan.plot(input.function, tissue, time.start, time.end,
    plot = TRUE, ...)
}
\arguments{
  \item{input.function}{Input function TAC.}

  \item{tissue}{Tissue TAC.}

  \item{time.start}{Initial acquisition time for each
  frame.}

  \item{time.end}{Final acquisition time for each frame.}

  \item{plot}{Should the Logan plot be displayed? Defaults
  to \code{TRUE}.}

  \item{...}{Additional parameters passed to the
  \code{plot} function.}
}
\description{
  Implements a Logan plot, or Logan graphical analysis, a
  "graphical method of analysis applicable to ligands that
  bind reversibly to receptors or enzymes requiring the
  simultaneous measurement of plasma and tissue
  radioactivities for multiple times after the injection of
  a radiolabeled tracer is presented".
}
\references{
  J. Logan, J. S. Fowler, N. D. Volkow, A. P. Wolf, S. L.
  Dewey, D. J. Schlyer, R. R. MacGregor, R. Hitzemann, B.
  Bendriem, and S. J. Gatley, "Graphical analysis of
  reversible radioligand binding from time-activity
  measurements applied to [N-11C-methyl]-(-)-cocaine PET
  studies in human subjects.," Journal of cerebral blood
  flow and metabolism : official journal of the
  International Society of Cerebral Blood Flow and
  Metabolism, vol. 10, no. 5, pp. 740-7, Sep. 1990.
}

