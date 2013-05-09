#' Implements a reversible one-tissue compartment model
#' 
#' @param input.function Input function TAC.
#' @param K1 K1 parameter.
#' @param k2 k2 parameter.
#' @param vB Fraction of blood in tissue.
#' @param time.start Initial acquisition time for each frame (in minutes).
#' @param time.end Final acquisition time for each frame (in minutes).
#' @param left.ventricle Left ventricle TAC, used for spill-over corretion.
#'   Defaults to \code{input.function}.
#' @param interpolation.type Interpolation type selection. Passed to
#'   \code{\link{interpolate.tac}}. Defaults to 1.
#' 
#' @return The TAC resulting of solving the model with the given parameters.

reversible.1c.model <- function(input.function, K1, k2, vB, 
                                time.start, time.end, 
                                left.ventricle = input.function,
                                interpolation.type = 1) {
    
    inter <- interpolate.tac(input.function, time.start, time.end, 
                             interpolation.type)    
    
    tsample <- inter$tsample
    input.function.inter <- inter$y
            
    # Convolve and return solution
    dt <- tsample[2] - tsample[1] # dt = frame length
    sol <- dt * convolve(K1*exp(-k2*tsample), rev(input.function.inter), 
                         type = "o")
    tsol <- seq(0, 2*max(tsample), dt)    
    (1 - vB) * revinterpolate.tac(sol, tsol, time.end) + vB * left.ventricle
}

#' Reversible one-tissue compartment model fit
#' 
#' Fits the model implemented in \code{\link{reversible.1c.model}}.
#' 
#' @param input.function Input function TAC.
#' @param tissue Tissue TAC.
#' @param time.start Initial acquisition time for each frame (in minutes).
#' @param time.end Final acquisition time for each frame (in minutes).
#' @param K1.start,k2.start,vB.start Initial parameter values.
#' @param K1.lower,k2.lower,vB.lower Parameter lower bounds.
#' @param K1.upper,k2.upper,vB.upper Parameter upper bounds.
#' @param left.ventricle Left ventricle TAC, used for spill-over corretion.
#'   Defaults to \code{input.function}.
#' @param weight Weight vector for the non-linear squares fit. Defaults to
#'  frame length.
#' @param plot \code{TRUE} if a plot is to be shown. Defaults to \code{FALSE}.
#' @param interpolation.type Interpolation type selection. Passed to
#'   \code{\link{interpolate.tac}}. Defaults to 1.
#' @param ... Other parameters passed to the plot function, if used.
#' 
#' @return Returns a list with three fields: \code{kparms}, the computed kinetic
#'   parameters; \code{stderrors}, the standard errors for each parameter as a 
#'   percentage; \code{fit}, the actual fitted object.

reversible.1c.fit <- function(input.function, tissue, time.start, time.end, 
                              K1.start = 0.5, K1.lower = 0, K1.upper = 1,
                              k2.start = 0.5, k2.lower = 0, k2.upper = 1,
                              vB.start = 0.05, vB.lower = 0, vB.upper = 1,
                              left.ventricle = input.function,
                              weight = time.end - time.start, plot = FALSE, 
                              interpolation.type = 1, ...) {   
        
    fit <- nlsLM(tissue ~ reversible.1c.model(input.function, K1, k2, 
                                              vB, time.start, time.end,
                                              left.ventricle,
                                              interpolation.type),
                 start = c(K1 = K1.start, k2 = k2.start, 
                           vB = vB.start),
                 weights = weight,                        
                 lower = c(K1 = K1.lower, k2 = k2.lower, 
                           vB = vB.lower),
                 upper = c(K1 = K1.upper, k2 = k2.upper, 
                           vB = vB.upper),
                 control = nls.lm.control(maxiter = 200))
    
    # Get parameters and standard errors.
    kparms <- coef(fit)    
    stderrors <- summary(fit)$coefficients[, 2]    
    stderrors <- (stderrors / kparms) * 100    
    
    if (plot == TRUE) {
        plotfit(tissue, fit, time.start, ...)
    }
    
    # Return fit parameters and fit object
    list(kparms = as.data.frame(t(kparms)), 
         stderrors = as.data.frame(t(stderrors)),
         fit = fit)
}
