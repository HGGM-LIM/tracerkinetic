#' Multilinear analysis 1 (MA1) for reversible 2-tissue model.
#'
#' This method needs to find the optimal linear point. This is yet to be 
#' implemented and this method should not be used. The code is nevertheless
#' included in this file for completeness.
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
    
    ma1 <- lm(tissue ~ x1 + x2)
    
    coefs <- as.numeric(coef(ma1))
    t1 <- coefs[2]
    t2 <- coefs[3]
    Vt <- -t1 / t2
    
    return(list(kparms = Vt, fit = ma1))    
}

#' Multilinear analysis 2 (MA2) for reversible 2-tissue model.
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
MA2 <- function(input.function, tissue, time.start, time.end) {
    
    dt <- time.end - time.start
    
    c1 <- cumsum(input.function * dt)
    c2 <- cumsum(tissue * dt)
    
    x1 <- cumsum(c1 * dt)
    x2 <- cumsum(c2 * dt)
    x3 <- c1
    x4 <- c2
    
    ma2 <- lm(tissue ~ x1 + x2 + x3 + x4)
    
    coefs <- as.numeric(coef(ma2))
    g1 <- coefs[2]
    g2 <- coefs[3]
    g3 <- coefs[4]
    g4 <- coefs[5]    
    Vt <- -g1 / g2
    
    return(list(kparms = Vt, fit = ma2))
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
    
    mlair1 <- lm(tissue ~ x1 + x2 + x3 + x4)
    
    coefs <- as.numeric(coef(mlair1))
    P1 <- coefs[2]
    P2 <- coefs[3]
    P3 <- coefs[4]
    P4 <- coefs[5]
    Ki <- - P4 / P3
    
    return(list(kparms = c(Ki, P1, P2, P3, P4), fit = mlair1)) 
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
    
    mlair2 <- lm(y ~ x1 + x2 + x3 + x4)
    coefs <- as.numeric(coef(mlair2))
    P1 <- coefs[2]
    P2 <- coefs[3]
    P3 <- coefs[4]
    P4 <- coefs[5]        
    Ki <- P1
    
    return(list(kparms = c(Ki, P1, P2, P3, P4), fit = mlair2))    
}
