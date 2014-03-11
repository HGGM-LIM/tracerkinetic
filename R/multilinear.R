#' Multilinear analysis 1 (MA1) for reversible 2-tissue model.
#'
#' @param input.function Input function TAC.
#' @param tissue Tissue TAC.
#' @param time.start Initial acquisition time for each frame (in minutes).
#' @param time.end Final acquisition time for each frame (in minutes).
#' 
#' @return A list with two fields: \code{kparms}, that contains the kinetic
#'   parameters, and \code{fit}, which contains the linear least squares fit.
#'   
#' @references M. Ichise, H. Toyama, R. Innis, and R. Carson, "Strategies to
#'   improve neuroreceptor parameter estimation by linear regression analysis,"
#'   J. Cereb. Blood Flow Metab., vol. 22, pp. 1271-1281, 2002.

MA1 <- function(input.function, tissue, time.start, time.end) {
    
    dt <- time.end - time.start
    
    x1 <- cumsum(input.function * dt)
    x2 <- cumsum(tissue * dt)
    
    # This method needs to find the linear part of the equation (T > t*). We
    # proceed using the same for loop that is used in the findbestfit
    # method
    maxerror <- 0.10
    minpoints <- 3
    fitdata <- data.frame(x1, x2, tissue)
    n <- nrow(fitdata)
    limit <- n - minpoints + 1    
    
    for (i in 1:limit) {
        fd <- fitdata[i:n, ]              
        ma1 <- lm(tissue ~ 0 + x1 + x2, data = fd)                  
        maxresid <- max(abs(ma1$residuals / fd$tissue))
        if (maxresid < maxerror) {            
            break
        }        
    }            
    
    coefs <- as.numeric(coef(ma1))    
    g1 <- coefs[1]
    g2 <- coefs[2]
    Vt <- -g1 / g2
    
    return(list(kparms = data.frame(Vt, g1, g2), fit = ma1))    
}

#' Multilinear analysis 2 (MA2) for reversible 2-tissue model.
#' 
#' @param input.function Input function TAC.
#' @param tissue Tissue TAC.
#' @param time.start Initial acquisition time for each frame (in minutes).
#' @param time.end Final acquisition time for each frame (in minutes).
#' 
#' @return A list with two fields: \code{kparms}, that contains the kinetic
#'   parameters (Vt, Vs, Vn), and \code{fit}, which contains the linear least 
#'   squares fit.
#'   
#' @references M. Ichise, H. Toyama, R. Innis, and R. Carson, "Strategies to
#'   improve neuroreceptor parameter estimation by linear regression analysis,"
#'   J. Cereb. Blood Flow Metab., vol. 22, pp. 1271-1281, 2002.
MA2 <- function(input.function, tissue, time.start, time.end) {
    
    dt <- time.end - time.start
    
    c1 <- cumsum(input.function * dt)
    c2 <- cumsum(tissue * dt) 
    
    x1 <- cumsum(c1 * dt)
    x2 <- cumsum(c2 * dt)
    x3 <- c2
    x4 <- c1
    
    ma2 <- lm(tissue ~ 0 + x1 + x2 + x3 + x4)    
    coefs <- as.numeric(coef(ma2))
    
    g1 <- coefs[1]
    g2 <- coefs[2]
    g3 <- coefs[3]
    g4 <- coefs[4]    
    Vt <- -g1 / g2
    Vs <- (-g1*(g1 + g3*g4) + g2*g4*g4)/(g2*(g1 + g3*g4))
    Vn <- Vt - Vs
    
    # This method might allow to extract the individual K rates
    # from the "g" parameters. This formulation is stille pending 
    # a validation and should NOT be used.
    K1 <- g4
    k3 <- -Vs*g2/g4
    k4 <- -1/((-(g3 + k3)/g2) + (Vt / K1) + k3/g2)
    k2 <- -g2/k4
    
    return(list(kparms = data.frame(Vt, Vs, Vn, K1, k2, k3, k4, g1, g2, g3, g4), 
                fit = ma2))
}

#' Implements the Multiple linear analysis for irreversible radiotracer 1
#' (MLAIR1)
#' 
#' @param input.function Input function TAC.
#' @param tissue Tissue TAC.
#' @param time.start Initial acquisition time for each frame (in minutes).
#' @param time.end Final acquisition time for each frame (in minutes).
#' 
#' @return A list with two fields: \code{kparms}, that contains the kinetic
#'   parameters, and \code{fit}, which contains the linear least squares fit.
#' 
#' @references S. J. Kim, J. S. Lee, Y. K. Kim, J. Frost, G. Wand, M. E. McCaul,
#'   and D. S. Lee, "Multiple linear analysis methods for the quantification of
#'   irreversibly binding radiotracers.," J. Cereb. Blood Flow Metab., vol. 28,
#'   no. 12, pp. 1965-77, Dec. 2008.
MLAIR1 <- function(input.function, tissue, time.start, time.end) {
    
    dt <- time.end - time.start
    
    c1 <- cumsum(input.function * dt)
    
    x1 <- input.function
    x2 <- c1
    x3 <- cumsum(tissue * dt)
    x4 <- cumsum(c1 * dt)
    
    mlair1 <- lm(tissue ~ 0 + x1 + x2 + x3 + x4)    
    coefs <- as.numeric(coef(mlair1))
        
    P1 <- coefs[1]
    P2 <- coefs[2]
    P3 <- coefs[3]
    P4 <- coefs[4]
    Ki <- - P4 / P3
    
    return(list(kparms = data.frame(Ki, P1, P2, P3, P4), fit = mlair1)) 
}

#' Implements the Multiple linear analysis for irreversible radiotracer 2
#' (MLAIR2)
#' 
#' @param input.function Input function TAC.
#' @param tissue Tissue TAC.
#' @param time.start Initial acquisition time for each frame (in minutes).
#' @param time.end Final acquisition time for each frame (in minutes).
#' 
#' @return A list with two fields: \code{kparms}, that contains the kinetic
#'   parameters, and \code{fit}, which contains the linear least squares fit.
#' 
#' @references S. J. Kim, J. S. Lee, Y. K. Kim, J. Frost, G. Wand, M. E. McCaul,
#'   and D. S. Lee, "Multiple linear analysis methods for the quantification of
#'   irreversibly binding radiotracers.," J. Cereb. Blood Flow Metab., vol. 28,
#'   no. 12, pp. 1965-77, Dec. 2008.
MLAIR2 <- function(input.function, tissue, time.start, time.end) {
    
    dt <- time.end - time.start
    
    c1 <- cumsum(input.function * dt)
    
    y <- cumsum(tissue * dt)
    x1 <- cumsum(c1 * dt)
    x2 <- c1
    x3 <- input.function
    x4 <- tissue
    
    mlair2 <- lm(y ~ 0 + x1 + x2 + x3 + x4)
    coefs <- as.numeric(coef(mlair2))
        
    P1 <- coefs[1]
    P2 <- coefs[2]
    P3 <- coefs[3]
    P4 <- coefs[4]        
    Ki <- P1
    
    return(list(kparms = data.frame(Ki, P1, P2, P3, P4), fit = mlair2))    
}

#' Multilinear analysis for 1-compartment model.
#' 
#' @param input.function Input function TAC.
#' @param tissue Tissue TAC.
#' @param time.start Initial acquisition time for each frame (in minutes).
#' @param time.end Final acquisition time for each frame (in minutes).
#' 
#' @return A list with two fields: \code{kparms}, that contains the kinetic
#'   parameters, and \code{fit}, which contains the linear least squares fit.
#' 
#' @references K. Murase, "Efficient method for calculating kinetic parameters 
#'   using T1-weighted dynamic contrast-enhanced magnetic resonance imaging.," 
#'   Magn. Reson. Med., vol. 51, no. 4, pp. 858-62, Apr. 2004.

murase <- function(input.function, tissue, time.start, time.end) {
    
    dt <- time.end - time.start
    x1 <- cumsum(input.function * dt)
    x2 <- cumsum(tissue * dt)
    x3 <- input.function
        
    lm1 <- nnls(cbind(x1, -x2, x3), tissue)
    coefs <- lm1$x
    
    K1 <- coefs[1] - coefs[2]*coefs[3]
    k2 <- coefs[2]
    vB <- coefs[3]
    
    return(list(kparms = data.frame(K1, k2, vB), fit = lm1))  
}
