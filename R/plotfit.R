#' Kinetic model fit plot
#' 
#' Plots the results of applying a given kinetic model to a tissue.
#' 
#' @param tissue The tissue TAC.
#' @param fit A \code{fit} object from one of the \code{.fit} functions.
#' @param time A time vector for plotting (x-axis).
#' @param legend.pos Position for the legend in the plot. Defaults to 
#'  \code{"bottomright"}. \code{NULL} to disable.
#' @param ... Other parameters passed to the plot function.

plotfit <- function(tissue, fit, time, legend.pos = "bottomright", ...) {    
    # Plots the results of fitting a given tissue TAC.
    #
    # Parameters:
    # - tissue: tissue TAC.
    # - fit: results of the NLS fit.
    # - time: time points for plotting (x axis).
    
    time.x <- c(time, time)
    activity <- c(tissue, predict(fit))
    plot(time.x, activity, type = "n", ...)
    grid()
    lines(time, predict(fit), type = "o", pch = 19)
    points(time, tissue)        
    if (!is.null(legend.pos)) {
        legend(legend.pos, c("Measured data", "Model"), 
               lty = c(NA, 1), pch = c(1, 19), inset = 0.02)
    }
}
