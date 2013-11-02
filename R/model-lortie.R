#' Reversible one-tissue compartment 82-rubidium model with flow correction
#' 
#' This function implements the one-tissue compartment rubidium-82 model 
#' published in [1].
#' 
#' @param input.function Input function TAC.
#' @param FLV Fraction of blood (from left ventricle) in myocardium TAC.
#' @param flow Myocardial blood flow.
#' @param k2 Rate constant.
#' @param time.start Initial acquisition time for each frame (in minutes).
#' @param time.end Final acquisition time for each frame (in minutes).
#' @param a Parameter \code{a} of the Renkin-Crone equation implemented in 
#'   \code{\link{flow2k1}}. Defaults to 0.77 [unitless].
#' @param b Parameter \code{b} of the Renkin-Crone equation implemented in 
#'   \code{\link{flow2k1}}. Defaults to 0.63 [ml/min/g].
#' @param left.ventricle Left ventricle TAC, used for spill-over corretion.
#'   Defaults to \code{input.function}.
#' @param interpolation.type Interpolation type selection. Passed to
#'   \code{\link{interpolate.tac}}. Defaults to 1.
#'   
#' @return The TAC resulting of solving the model with the given parameters.
#'   
#' @references [1] M. Lortie, R. S. B. Beanlands, K. Yoshinaga, R. Klein, J. N. 
#'   Dasilva, and R. a DeKemp, "Quantification of myocardial blood flow with 
#'   82Rb dynamic PET imaging.," European journal of nuclear medicine and 
#'   molecular imaging, vol. 34, no. 11, pp. 1765-74, Nov. 2007.

lortie.model <- function(input.function, FLV, flow, k2, time.start, time.end, 
                         a = 0.77, b = 0.63, left.ventricle = input.function, 
                         interpolation.type = 1) {        
       
    inter <- interpolate.tac(input.function, time.start, time.end,
                             interpolation.type)
    
    tsample <- inter$tsample    
    input.function.inter <- inter$y        
        
    # Convolve and return solution.
    # Uses flow for fitting instead of K1, as PMOD does. K1 value is 
    # computed from the flow parameter afterwards (see flow2k1 function).
    # Note that the flow2k1 function expects the flow in ml/min/g and this 
    # function computes the convolution in seconds, so unit conversion is needed
    # both for input and output (K1 is needed here in ml/min/g).
    K1 <- flow2k1(flow, a, b)
    dt <- tsample[2] - tsample[1] # dt = frame length
    sol <- convolve(K1*exp(-k2*tsample), dt * rev(input.function.inter), 
                         type = "o")        
    
    tsol <- seq(0, 2*max(tsample), dt)    
    FLV * left.ventricle + 
        (1 - FLV) * revinterpolate.tac(sol, tsol, time.end)
}    

#' K1 to flow conversion
#' 
#' Implements the Renkin-Crone equation that transforms flow values to K1
#' in the 82-rubidium tracer.
#' 
#' @param flow The myocardial blood flow
#' @param a Parameter \code{a} of the Renkin-Crone equation. Defaults to 0.77 
#'  [unitless].
#' @param b Parameter \code{b} of the Renkin-Crone equation. Defaults to 0.63
#'  [ml/min/g].
#'  
#' @return The K1 parameter value.

flow2k1 <- function(flow, a = 0.77, b = 0.63) {    
    (1 - a*exp(-b/flow))*flow     
}

#' Lortie et al. (2007) model fit.
#' 
#' Fits the one-tissue compatment 82-rubidium model implemented in
#' \code{\link{lortie.model}}.
#' 
#' @param input.function Input function TAC (typically, left ventricle).
#' @param tissue Tissue (myocardium) TAC.
#' @param time.start Initial time for each frame (in minutes).
#' @param time.end End time for each frame (in minutes).
#' @param a Parameter \code{a} of the Renkin-Crone equation implemented in 
#'   \code{\link{flow2k1}}. Defaults to 0.77 [unitless].
#' @param b Parameter \code{b} of the Renkin-Crone equation implemented in 
#'   \code{\link{flow2k1}}. Defaults to 0.63 [ml/min/g].
#' @param FLV.start,flow.start,k2.start Initial parameter values for the 
#'   non-linear squares fitting.
#' @param FLV.lower,flow.lower,k2.lower Lower bounds.
#' @param FLV.upper,flow.upper,k2.upper Upper bounds.
#' @param left.ventricle Left ventricle TAC, used for spill-over corretion.
#'   Defaults to \code{input.function}.
#' @param weight Weights for the non-linear least squares. Defaults to no
#'   weighting. Use \code{"framelength"} to use the frame length as the weight
#'   value.
#' @param plot Should the results be plotted. Defaults to \code{FALSE}.
#' @param interpolation.type Interpolation type selection. Passed to
#'   \code{\link{interpolate.tac}}. Defaults to 1.
#' @param ... If plotting is used, optional parameters are passed to the plot 
#'   function.
#'   
#' @return Returns a list with four fields: \code{kparms}, the computed kinetic
#'   parameters; \code{stderrors}, the standard errors for each parameter;
#'   \code{stderrorsp}, the standard errors for each parameter as a percentage; 
#'   \code{fit}, the actual fitted object.

lortie.fit <- function(input.function, tissue, time.start, time.end, 
                       a = 0.77, b = 0.63,
                       FLV.start = 0.1, FLV.lower = 0, FLV.upper = 1,
                       flow.start = 0.5, flow.lower = 0, flow.upper = 8,
                       k2.start = 0.1, k2.lower = 0, k2.upper = 8,     
                       left.ventricle = input.function,
                       weight = NA, plot = FALSE, 
                       interpolation.type = 1, ...) {
    
    if (length(weight) == 1) {
        # Set appropriate weight values
        if (is.na(weight) | is.null(weight))    # Constant weighting
            weight <- rep(1, length(tissue))
        else if (weight == "framelength")       # Frame length
            weight <- time.end - time.start
    }
        
    fit <- nlsLM(tissue ~ lortie.model(input.function, FLV, flow, k2, 
                                       time.start, time.end, 
                                       left.ventricle = left.ventricle, 
                                       interpolation.type = interpolation.type), 
                 start = list(FLV = FLV.start, flow = flow.start, 
                              k2 = k2.start),
                 weights = weight,                   
                 lower = c(FLV = FLV.lower, flow = flow.lower, 
                           k2 = k2.lower),
                 upper = c(FLV = FLV.upper, flow = flow.upper, 
                           k2 = k2.upper),
                 control = nls.lm.control(maxiter = 200))               
    
    # Get parameters and standard errors.
    kparms <- coef(fit)    
    stderrors <- summary(fit)$coefficients[, 2]    
    stderrorsp <- (stderrors / kparms) * 100         
    
    # Reorder
    stderrors <- stderrors[c(2, 3, 1)]
    stderrorsp <- stderrorsp[c(2, 3, 1)]
    kparms <- c(flow2k1(kparms[2]), kparms[c(2, 3, 1)])
    names(kparms) <- c("K1", "flow", "k2", "FLV")
        
    if (plot == TRUE) {        
        plotfit(tissue, fit, time.start, ...)
    }
    
    # Return a list with the result    
    list(kparms = as.data.frame(t(kparms)), 
         stderrors = as.data.frame(stderrors), 
         stderrorsp = as.data.frame(stderrorsp),
         fit = fit)
}
