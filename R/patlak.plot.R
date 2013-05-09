#' Patlak plot implementation
#' 
#' Implements a Patlak plot linear analysis, another graphical method for
#' analysis of tracers that can be modeled after an irreversible two-tissue
#' compartment model.
#' 
#' @param input.function Input function TAC.
#' @param tissue Tissue TAC.
#' @param time.start Initial acquisition time for each frame.
#' @param time.end Final acquisition time for each frame.
#' @param plot Should the Patlak plot be displayed? Defaults to \code{TRUE}.
#' @param ... Additional parameters passed to the \code{plot} function.
#' 
#' @references TODO
#' 
#' @export

patlak.plot <- function(input.function, tissue, time.start, time.end, 
                        plot = TRUE, ...) {
        
    # Compute data on both axis
    termy <- tissue / input.function
    dt <- time.end - time.start
    termx <- cumsum(input.function * dt) / input.function
    
    # Another option is using interpolation. The result is almost the same.
    #iif <- interpolate.tac(input.function, time.start, time.end)
    #dt <- iif$tsample
    #itermx <- cumsum(iif$y)
    #termx <- revinterpolate.tac(itermx, dt, time.end)
    #termx <- termx / input.function
    
    # Correct bad points (just for aesthetic purposes)
    bad <- termx < 0 | termy < 0
    termx[bad] <- 0
    termy[bad] <- 0
        
    # The fitting and plotting is done by this helper function
    res <- processlinear(termx, termy, plot, main = "Patlak", 
                         xlab = "Int(Cplasma) / Cplasma [min]",
                         ylab = "Ctissue / Cplasma [unitless]", ...)
    
    return(res)
}
