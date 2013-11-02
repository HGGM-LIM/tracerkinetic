#' Helper function for linear models
#' 
#' This function helps process linear models, such as \code{\link{patlak.plot}}
#' or \code{\link{logan.plot}}. It finds the optimal fitting point by calling
#' \code{\link{findbestfit}} and proceeeds to compute the kinetic parameters,
#' that will mostly consist of the intercept and slope of the linear regression,
#' along with their standard errors.
#' 
#' It \code{plot} is set to \code{TRUE}, the regression will be shown on screen.
#' The data points used to compute the regression are displayed with a filled
#' style, while the ones that have not been used are displayed using open
#' circles.
#' 
#' @param termx,termy The x and y values of the graph.
#' @param plot Should the result be shown? Defaults to \code{TRUE}.
#' @param ... Further parameters for the \code{plot} function.
#' 
#' @return Returns a list with four fields: \code{kparms}, the computed kinetic
#'   parameters; \code{stderrors}, the standard errors for each parameter;
#'   \code{stderrorsp}, the standard errors for each parameter as a percentage; 
#'   \code{fit}, the actual fitted object.

processlinear <- function(termx, termy, plot = TRUE, ...) {
    
    fit <- findbestfit(termx, termy)
    t0 <- fit$initial.time.point
    
    coefs <- as.numeric(coef(fit$fit))
    res <- data.frame(intercept = coefs[1], slope = coefs[2])
    stdabs <- coef(summary(fit$fit))[, 2]
    stderr <- data.frame(intercept = stdabs[1], slope = stdabs[2])
    stderrors <- (coef(summary(fit$fit))[, 2] / abs(res)) * 100
    stderrp <- data.frame(intercept = stderrors[1], slope = stderrors[2])
    
    if (plot == TRUE) {
        # Plotting style
        pchs <- rep(1, length(termx))
        pchs[t0:length(pchs)] <- 19
        plot(termx, termy, pch = pchs, ...)
        abline(fit$fit)
    }
    
    return(list(kparms = res, stderrors = stderr, stderrorsp = stderrp, 
                fit = fit$fit))
    
}
