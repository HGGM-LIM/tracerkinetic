#' Finds optimal fitting point
#' 
#' Finds the number of points that have to be used for an optimal fit between
#' the \code{x} and \code{y} variables.
#' 
#' @param x,y Regression data
#' @param minpoints Minimum data points to be used for the regression.
#' @param maxerror The maximum error rate allowed.
#' 
#' @details This function finds the earliest point that allows to have a 
#'   regression with less than 10% error between the chosen points (this 
#'   parameter is controlled by the \code{maxerror} variable). If no point
#'   allows to have this error rate, the point that yields the minimum 
#'   regression error is used.
#'   
#' @return A list containing the time point chosen, \code{initial.time.point}, 
#'   and the actual fitted object, \code{fit}.

findbestfit <- function(x, y, minpoints = 3, maxerror = 0.10) {
        
    fitdata <- data.frame(x, y)
    n <- nrow(fitdata)
    res <- rep(NA, n - minpoints + 1)
    
    for (i in 1:(n - minpoints + 1)) {
        fd <- fitdata[i:n, ]        
        lm1 <- lm(y ~ x, data = fd)    
        res[i] <- max(lm1$residuals / fd$y)
    }            
    
    suppressWarnings(respoint <- min(which(res < maxerror)))
    if (respoint == Inf) 
        respoint <- which(res == min(res))        
    lmres <- lm(y ~ x, data = fitdata[respoint:n, ])
    return(list(initial.time.point = respoint, fit = lmres))  
    
}
