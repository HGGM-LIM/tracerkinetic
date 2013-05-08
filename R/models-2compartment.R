#' Implements a reversible two-tissue compartment model
#' 
#' @param input.function Input function TAC.
#' @param K1 K1 parameter.
#' @param k2 k2 parameter.
#' @param k3 k3 parameter.
#' @param k4 k4 parameter.
#' @param vB Fraction of blood in tissue.
#' @param time.start Initial acquisition time for each frame.
#' @param time.end Final acquisition time for each frame.
#' @param left.ventricle Left ventricle TAC, used for spill-over corretion.
#'   Defaults to \code{input.function}.
#' @param interpolation.type Interpolation type selection. Passed to
#'   \code{\link{interpolate.tac}}. Defaults to 1.
#' 
#' @return The TAC resulting of solving the model with the given parameters.
#' 
#' @export
reversible.2c.model <- function(input.function, K1, k2, k3, k4, vB, 
                                time.start, time.end, 
                                left.ventricle = input.function, 
                                interpolation.type = 1) {
    
    inter <- interpolate.tac(input.function, time.start, time.end, 
                             interpolation.type)    
    tsample <- inter$tsample
    input.function.inter <- inter$y
        
    # Equations obtained from M. Bentourkia dn H. Zaidi, "Quantitative Analysis
    # in Nuclear Medicine Imaging", p. 398. Adapted to make the spill-over term
    # a true fraction.
    delta <- sqrt((k2 + k3 + k4)^2 - 4 * k2 * k4)
    a1 <- (k2 + k3  + k4 - delta) / 2
    a2 <- (k2 + k3  + k4 + delta) / 2
    
    term1 <- (K1 / (a2 - a1))*((k3 + k4 - a1)*exp(-a1*tsample) + 
                                   (a2 - k3 - k4)*exp(-a2*tsample))
    
    if(all(is.na(term1))) term1 <- rep(0, length(term1))
    
    # Equations obtained from Gunn et al., IEEE MIC 2011 - Kinetical Modeling,
    # slide 42.
    # Yield exactly the same results as the above implementation.
    # delta  <- same as above
    # theta1 <- a2
    # theta2 <- a1
    # phi1 <- (K1 * (theta1 - k3 - k4)) / delta
    # phi2 <- (K1 * (theta2 - k3 - k4)) / -delta    
    # term1 <- phi1 * exp(-theta1 * time.vector) + 
    #          phi2 * exp(-theta2 * time.vector)    
        
    term2 <- rev(input.function.inter)
    dt <- tsample[2] - tsample[1] # dt = frame length
    sol <- dt * convolve(term1, term2, type = "o")
    
    tsol <- seq(0, 2*max(tsample), dt) 
    (1 - vB) * revinterpolate.tac(sol, tsol, time.end) + 
        vB * left.ventricle
}

#' Reversible two-tissue compartment model fit
#' 
#' Fits the model implemented in \code{\link{reversible.2c.model}}.
#' 
#' @param input.function Input function TAC.
#' @param tissue Tissue TAC.
#' @param time.start Initial acquisition time for each frame.
#' @param time.end Final acquisition time for each frame.
#' @param K1.start,k2.start,k3.start,k4.start,vB.start Initial parameter values.
#' @param K1.lower,k2.lower,k3.lower,k4.lower,vB.lower Parameter lower bounds.
#' @param K1.upper,k2.upper,k3.upper,k4.upper,vB.upper Parameter upper bounds.
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
#'   
#' @export
reversible.2c.fit <- function(input.function, tissue, time.start, time.end, 
                              K1.start = 0.1, K1.lower = 0, K1.upper = 8,
                              k2.start = 0.1, k2.lower = 0, k2.upper = 8,
                              k3.start = 0.01, k3.lower = 0, k3.upper = 8,
                              k4.start = 0.01, k4.lower = 0, k4.upper = 8,
                              vB.start = 0.05, vB.lower = 0, vB.upper = 1,
                              left.ventricle = input.function, 
                              weight = time.end - time.start, plot = FALSE, 
                              interpolation.type = 1, ...) {
    
    fit <- nlsLM(tissue ~ reversible.2c.model(input.function, K1, k2, k3, k4, 
                                              vB, time.start, time.end, 
                                              left.ventricle,
                                              interpolation.type),
                 start = c(K1 = K1.start, k2 = k2.start, 
                           k3 = k3.start, k4 = k3.start, 
                           vB = vB.start),
                 weights = weight,                        
                 lower = c(K1 = K1.lower, k2 = k2.lower,
                           k3 = k3.lower, k4 = k4.lower, 
                          vB = vB.lower),
                 upper = c(K1 = K1.upper, k2 = k2.upper, 
                           k3 = k3.upper, k4 = k4.upper, 
                           vB = vB.upper),                        
                 control = nls.lm.control(maxiter = 200))
        
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

#' Implements an irreversible two-tissue compartment model
#' 
#' This model is basically the same as \code{\link{reversible.2c.model}} but
#' the parameter \code{k4} is always 0. The simplified equation has been
#' implemented in this function and can be fitterd using 
#' \code{\link{irreversible.2c.fit}}.
#' 
#' @param input.function Input function TAC.
#' @param K1 K1 parameter.
#' @param k2 k2 parameter.
#' @param k3 k3 parameter.
#' @param vB Fraction of blood in tissue.
#' @param time.start Initial acquisition time for each frame.
#' @param time.end Final acquisition time for each frame.
#' @param left.ventricle Left ventricle TAC, used for spill-over corretion.
#'   Defaults to \code{input.function}.
#' @param interpolation.type Interpolation type selection. Passed to
#'   \code{\link{interpolate.tac}}. Defaults to 1.
#' 
#' @return The TAC resulting of solving the model with the given parameters.
#' 
#' @export
irreversible.2c.model <- function(input.function, K1, k2, k3, vB,
                                  time.start, time.end,
                                  left.ventricle = input.function,
                                  interpolation.type = 1) {
    
    inter <- interpolate.tac(input.function, time.start, time.end, 
                             interpolation.type)    
    tsample <- inter$tsample
    input.function.inter <- inter$y    
    
    term1 <- ((K1 * k2) / (k2 + k3)) * exp(tsample * -(k2 + k3)) +
             (K1 * k3) / (k2 + k3)        
    term2 <- rev(input.function.inter)
        
    if(all(is.na(term1))) term1 <- rep(0, length(term1))    
    
    dt <- tsample[2] - tsample[1] # dt = frame length
    sol <- dt * convolve(term1, term2, type = "o")    
    tsol <- seq(0, 2*max(tsample), dt)     
    (1 - vB) * revinterpolate.tac(sol, tsol, time.end) + 
        vB * left.ventricle
}

#' Irreversible two-tissue compartment model fit
#' 
#' Fits the model implemented in \code{\link{irreversible.2c.model}}.
#' 
#' @param input.function Input function TAC.
#' @param tissue Tissue TAC.
#' @param time.start Initial acquisition time for each frame.
#' @param time.end Final acquisition time for each frame.
#' @param K1.start,k2.start,k3.start,vB.start Initial parameter values.
#' @param K1.lower,k2.lower,k3.lower,vB.lower Parameter lower bounds.
#' @param K1.upper,k2.upper,k3.upper,vB.upper Parameter upper bounds.
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
#'   
#' @export
irreversible.2c.fit <- function(input.function, tissue, time.start, time.end, 
                                K1.start = 0.1, K1.lower = 0, K1.upper = 8,
                                k2.start = 0.1, k2.lower = 0, k2.upper = 8,
                                k3.start = 0.1, k3.lower = 0, k3.upper = 8,
                                vB.start = 0.05, vB.lower = 0, vB.upper = 1,
                                left.ventricle = input.function,
                                weight = time.end - time.start, plot = FALSE, 
                                interpolation.type = 1, ...) {
            
    fit <- nlsLM(tissue ~ irreversible.2c.model(input.function, K1, k2, k3, 
                                                vB, time.start, time.end, 
                                                left.ventricle,
                                                interpolation.type),
                 start = c(K1 = K1.start, k2 = k2.start, 
                           k3 = k3.start, vB = vB.start),
                 weights = weight,                        
                 lower = c(K1 = K1.lower, k2 = k2.lower,
                           k3 = k3.lower, vB = vB.lower),
                 upper = c(K1 = K1.upper, k2 = k2.upper, 
                           k3 = k3.upper, vB = vB.upper),                        
                 control = nls.lm.control(maxiter = 200))
    
    kparms <- coef(fit)        
    stderrors <- summary(fit)$coefficients[, 2]    
    stderrors <- (stderrors / kparms) * 100       
    
    if (plot == TRUE) {
        plotfit(tissue, fit, time.start, ...)
    }
    
    # Return fit parameters and fit object
    list(kparms = as.data.frame(t(kparms)), 
         stderrors = stderrors,
         fit = fit)
}

#' Computes the one-tissue compartmental model parameters for the 13-ammonia
#' tracer. Automatically performs the metabolite correction implemented in 
#' \code{\link{nh3.vdhoff.correction}} to the \code{input.function} term. The
#' \code{left.ventricle} term is not provided. All the other parameters are
#' provided just like in the \code{\link{irreversible.2c.fit}} function and the
#' return values are the same.
#' 
#' @param input.function Input function TAC.
#' @param tissue Tissue TAC.
#' @param time.start Initial acquisition time for each frame.
#' @param time.end Final acquisition time for each frame.
#' @param K1.start,k2.start,k3.start,vB.start Initial parameter values.
#' @param K1.lower,k2.lower,k3.lower,vB.lower Parameter lower bounds.
#' @param K1.upper,k2.upper,k3.upper,vB.upper Parameter upper bounds.
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
#'   
#' @export
#' 
#' @export
nh3.2c.fit <- function(input.function, tissue, time.start, time.end, 
                       K1.start = 0.1, K1.lower = 0, K1.upper = 8,
                       k2.start = 0.1, k2.lower = 0, k2.upper = 8,
                       k3.start = 0.1, k3.lower = 0, k3.upper = 8,
                       vB.start = 0.05, vB.lower = 0, vB.upper = 1,                       
                       weight = time.end - time.start, plot = FALSE, 
                       interpolation.type = 1, ...) {
    
    corr.input.function <- nh3.vdhoff.correction(input.function, time.end)    
    irreversible.2c.fit(corr.input.function, tissue, time.start, time.end,
                        K1.start, K1.lower, K1.upper, 
                        k2.start, k2.lower, k2.upper,
                        k3.start, k3.lower, k3.upper,
                        vB.start, vB.lower, vB.upper,
                        left.ventricle = input.function, weight = weight,
                        plot = plot, interpolation.type = interpolation.type, 
                        ...)
}
