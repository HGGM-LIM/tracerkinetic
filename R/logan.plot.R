#' Logan plot implementation
#' 
#' Implements a Logan plot, or Logan graphical analysis, a "graphical method
#' of analysis applicable to ligands that bind reversibly  to receptors or
#' enzymes requiring the simultaneous measurement of plasma and tissue
#' radioactivities for multiple times after the injection of a radiolabeled 
#' tracer is presented".
#' 
#' @param input.function Input function TAC.
#' @param tissue Tissue TAC.
#' @param time.start Initial acquisition time for each frame.
#' @param time.end Final acquisition time for each frame.
#' @param plot Should the Logan plot be displayed? Defaults to \code{TRUE}.
#' @param ... Additional parameters passed to the \code{plot} function.
#' 
#' @references J. Logan, J. S. Fowler, N. D. Volkow, A. P. Wolf, S. L.
#' Dewey, D. J. Schlyer, R. R. MacGregor, R. Hitzemann, B. Bendriem, and S. J.
#' Gatley, "Graphical analysis of reversible radioligand binding from
#' time-activity measurements applied to [N-11C-methyl]-(-)-cocaine PET studies
#' in human subjects.," Journal of cerebral blood flow and metabolism : official
#' journal of the International Society of Cerebral Blood Flow and Metabolism,
#' vol. 10, no. 5, pp. 740-7, Sep. 1990.
#' 
#' @export
logan.plot <- function(input.function, tissue, time.start, time.end, 
                       plot = TRUE, ...) {
    
    tmid <- (time.start + time.end) / 2    
    input.function <- interpolate.tac(input.function, time.start, time.end)$y      
    tissue <- interpolate.tac(tissue, time.start, time.end)$y        
    input.function <- input.function[tmid + 1]
    tissue <- tissue[tmid + 1]                     
                                         
    frame.length <- time.end - time.start
    term1 <- cumsum(tissue * frame.length) / tissue 
    term2 <- cumsum(input.function * frame.length) / tissue             
    #term1 <- cumsum(tissue) / tissue 
    #term2 <- cumsum(input.function) / tissue             
    #term1 <- term1[tmid + 1]
    #term2 <- term2[tmid + 1]
    #l <- length(term1)    
    l <- length(term1)    
            
    # Find optimal time.regression point
    error <- rep(NA, l)
    time.regression <- l - 1
    
    for (i in 1:(l-2)) {
        lm1 <- lm(term1[(l - time.regression):l] ~ 
                      term2[(l - time.regression):l], 
                  weights = frame.length[(l - time.regression):l]) 
        error[i] <- summary(lm1)$sigma        
        time.regression <- time.regression - 1
    }        
    
    time.regression <- l - which(error == min(error, na.rm = TRUE))    
    
    # Do the regression from the chosen time point
    lm1 <- lm(term1[(l-time.regression):l] ~ term2[(l-time.regression):l], 
              weights = frame.length[(l-time.regression):l]) 
        
    # Plot, if requested
    if (plot) {         
        # Filled circles for taken points, empty circles for the others
        point.types <- rep(1, l)
        point.types[(l - time.regression):l] <- 19              
        # Adjust margins
        mar <- par("mar")
        mar[2] <- mar[2] + 3        
        par.def <- par(mar = mar)        
        # Set the x and y axis labels to the correct equations
        x.exp <- expression(frac(integral(C[plasma](u) * du, 0, t),
                                 C[tissue](t))
        )
        y.exp <- expression(frac(integral(C[tissue](u) * du, 0, t), 
                                 C[tissue](t))
                            )            
        # Plot with additional options (if provided)
        plot(term2, term1, pch = point.types, xlab = x.exp, ylab = y.exp, ...)  
        # Show regression line
        abline(lm1)
        # Return plot parameters to default values
        par(par.def)
    }        
    
    kparms = data.frame(Vt = as.numeric(lm1$coefficients[2]),
                        intercept = as.numeric(lm1$coefficients[1]))
    
    # Return slope as Vt, intercept and fit object
    list(kparms = kparms, fit = lm1)
}
