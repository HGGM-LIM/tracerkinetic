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
#' @details Please refer to the documentation of \code{\link{processlinear}}
#'   for implementation details. 
#'   
#' @return Returns a list with three fields: \code{kparms}, the computed kinetic
#'   parameters; \code{stderrors}, the standard errors for each parameter as a 
#'   percentage; \code{fit}, the actual fitted object.
#'   
#' @references C. Patlak, R. Blasberg, and J. Fenstermacher, "Graphical
#'   evaluation of blood-to-brain transfer constants from multiple-time uptake
#'   data," J Cereb Blood Flow Metab, 1983.
#' 
#' @seealso \code{\link{logan.plot}}.

patlak.plot <- function(input.function, tissue, time.start, time.end, 
                        plot = TRUE, ...) {
        
    # Compute data on both axis
    termy <- tissue / input.function
    dt <- time.end - time.start
    termx <- cumsum(input.function * dt) / input.function
        
    # The fitting and plotting is done by this helper function
    res <- processlinear(termx, termy, plot, main = "Patlak", 
                         xlab = "Int(Cplasma) / Cplasma [min]",
                         ylab = "Ctissue / Cplasma [unitless]", ...)
    
    return(res)
}
