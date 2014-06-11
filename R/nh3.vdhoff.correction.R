#' NH3 metabolite correction
#' 
#' Implements the 13-NH3 metabolite correction as published in [1]. This 
#' correction should be applied to the input function before using it as the 
#' input for any model.
#' 
#' @param input.function The input function TAC.
#' @param time.end Final acquisition time for each frame.
#' @param t0 Empirical value for the t0 correction parameter. Defaults to 0.48
#'   [min].
#' @param T12 Empirical value for the T1/2 correction parameter. Defaults to
#'  6.69 [min].
#'   
#' @references [1] J. van den Hoff et al.,
#'   "[1-11C]Acetate as a Quantitative Perfusion Tracer in Myocardial PET," 
#'   Journal of Nuclear Medicine, vol. 42, no. 8, pp. 1174-1182, Aug. 2001.

nh3.vdhoff.correction <- function(input.function, time.end, 
                                  t0 = 0.48, T12 = 6.69) {
    
    b <- -log(2) / T12
    
    # Create corrected array    
    input.function.corr <- input.function    
    post <- (time.end > t0)
    input.function.corr[post] <- exp(b * (time.end[post] - t0)) * 
                                 input.function[post]        
       
    return(input.function.corr)
}
