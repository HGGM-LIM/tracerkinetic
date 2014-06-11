\name{plotfit}
\alias{plotfit}
\title{Kinetic model fit plot}
\usage{
  plotfit(tissue, fit, time, legend.pos = "bottomright",
    ...)
}
\arguments{
  \item{tissue}{The tissue TAC.}

  \item{fit}{A \code{fit} object from one of the
  \code{.fit} functions.}

  \item{time}{A time vector for plotting (x-axis).}

  \item{legend.pos}{Position for the legend in the plot.
  Defaults to \code{"bottomright"}. \code{NULL} to
  disable.}

  \item{...}{Other parameters passed to the plot function.}
}
\description{
  Plots the results of applying a given kinetic model to a
  tissue.
}

