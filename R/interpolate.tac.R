#' Interpolates the given TAC to perform the correlation.
#' 
#' This interpolation process (actually, an extrapolation) is needed to convert
#' the given TAC into a signal with frames of equal length. This length is the
#' minimum frame length for the given study.
#' 
#' @param tac The TAC to be interpolated.
#' @param time.start Initial acquisition time for each frame.
#' @param time.end Final acquisition time for each frame.
#' @param interpolation.type Interpolation type. 1 for linear, 2 for "natural" 
#'   splines. Defaults to 1.
#'   
#' @return A list containing \code{tsample}, the sampling time vector and the 
#'   \code{y}, the interpolated TAC.

interpolate.tac <- function(tac, time.start, time.end, 
                            interpolation.type = 1) {
    
    # Initial condition is 0 for both the TAC and the time.end values
    tac0 <- c(0, tac)
    time.end0 <- c(0, time.end)    
    # Compute frame length for next step
    frame.length <- time.end - time.start    
    # The time.vector variable divides all the frames into the minimum 
    # time length unit, which is a second (1/60 minutes).             
    time.vector <- seq(0, max(time.end), 1/60)          
    # Interpolate according to user selection
    if (interpolation.type == 1)
        tac.inter <- approx(time.end0, tac0, xout = time.vector, rule = 2)
    else if (interpolation.type == 2)
        tac.inter <- spline(time.end0, tac0, xout = time.vector, 
                            method = "natural")                             
        
    return(list(tsample = time.vector, y = tac.inter$y))
}

#' Actual interpolation.
#' 
#' Typically used after \code{\link{interpolate.tac}}. This is just a simplified
#' call to the \code{\link{approx}} function.
#' 
#' @param tac The extrapolated TAC to interpolate
#' @param tsample The sampling points used in the original extrapolation.
#' @param time.points The new time points for which the signal is desired.
#'   
#' @return The new interpolated TAC.

revinterpolate.tac <- function(tac, tsample, time.points) {
        
    approx(tsample, tac, xout = time.points, rule = 2)$y
    
}
