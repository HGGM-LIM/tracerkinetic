
processlinear <- function(termx, termy, plot = TRUE, ...) {
    
    fit <- findbestfit(termx, termy)
    t0 <- fit$initial.time.point
    
    coefs <- as.numeric(coef(fit$fit))
    res <- data.frame(intercept = coefs[1], slope = coefs[2])
    stderrors <- (coef(summary(fit$fit))[, 2] / res) * 100
    stderr <- data.frame(intercept = stderrors[1], slope = stderrors[2])
    
    if (plot == TRUE) {
        # Plotting style
        pchs <- rep(1, length(termx))
        pchs[t0:length(pchs)] <- 19
        plot(termx, termy, pch = pchs, ...)
        abline(fit$fit)
    }
    
    return(list(kparms = res, stderrors = stderr, fit = fit$fit))
    
}