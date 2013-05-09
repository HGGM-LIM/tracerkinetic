#' Finds optimal fitting 
#' 
#' Finds the number of points that have to be used for an optimal fit between
#' the \code{x} and \code{y} variables.

findbestfit <- function(x, y, minpoints = 3, maxerror = 0.05) {
    
    fitdata <- data.frame(x, y)
    n <- nrow(fitdata)
    res <- rep(NA, n - minpoints + 1)
    for (i in 1:(n - minpoints + 1)) {
        fd <- fitdata[i:n, ]        
        lm1 <- lm(y ~ x, data = fd)
        res[i] <- sum(lm1$residuals^2) / length(fd)
    }            
    
    suppressWarnings(respoint <- min(which(res < maxerror)))
    if (respoint == Inf) 
        respoint = which(res == min(res))
    lmres <- lm(y ~ x, data = fitdata[respoint:n, ])
    return(list(initial.time.point = respoint, fit = lmres))    
}
