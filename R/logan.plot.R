#' Logan plot implementation
#' 
#' Implements a Logan plot, or Logan graphical analysis, a "graphical method
#' of analysis applicable to ligands that bind reversibly  to receptors or
#' enzymes requiring the simultaneous measurement of plasma and tissue
#' radioactivities for multiple times after the injection of a radiolabeled 
#' tracer is presented".
#' 
#' @param input.function Input function TAC.
#' @param tissue Tissue TAC.
#' @param time.start Initial acquisition time for each frame.
#' @param time.end Final acquisition time for each frame.
#' @param plot Should the Logan plot be displayed? Defaults to \code{TRUE}.
#' @param ... Additional parameters passed to the \code{plot} function.
#' 
#' @details Please refer to the documentation of \code{\link{processlinear}}
#'   for implementation details. 
#' 
#' @return Returns a list with three fields: \code{kparms}, the computed kinetic
#'   parameters; \code{stderrors}, the standard errors for each parameter as a 
#'   percentage; \code{fit}, the actual fitted object.
#'   
#' @references J. Logan, J. S. Fowler, N. D. Volkow, A. P. Wolf, S. L.
#' Dewey, D. J. Schlyer, R. R. MacGregor, R. Hitzemann, B. Bendriem, and S. J.
#' Gatley, "Graphical analysis of reversible radioligand binding from
#' time-activity measurements applied to [N-11C-methyl]-(-)-cocaine PET studies
#' in human subjects.," Journal of cerebral blood flow and metabolism : official
#' journal of the International Society of Cerebral Blood Flow and Metabolism,
#' vol. 10, no. 5, pp. 740-7, Sep. 1990.
#' 
#' @seealso \code{\link{patlak.plot}}.

logan.plot <- function(input.function, tissue, time.start, time.end, 
                       plot = TRUE, ...) {
    
    # Compute data on both axis
    dt <- time.end - time.start
    termy <- cumsum(tissue * dt) / tissue
    termx <- cumsum(input.function * dt) / tissue
    
    # Correct bad points (just for aesthetic purposes)
    bad <- termx < 0 | termy < 0
    termx[bad] <- 0
    termy[bad] <- 0
    
    # The fitting and plotting is done by this helper function
    res <- processlinear(termx, termy, plot, main = "Logan", 
                         xlab = "Int(Cplasma) / Ctissue [min]",
                         ylab = "Int(Ctissue) / Ctissue [min]", ...)
        
    return(res)
}
