## (Generic) Function to calculate A values from microarray objects
##
## Parameters: raw    -> object of class maigesRaw or maiges
##             bkgSub -> string indicating the type of background
##                       subtraction. May be "none", "subtract",
##                       "half", "minimum", "movingmin", "edwards",
##                       "normexp" or "rma". Uses limma and defaults
##                       to "subtract".
##
## Gustavo H. Esteves
## 14/05/07
##
##


## Defining a default method
calcA.default <- function(object, ...)
    return(object@A)


## Defining a maigesRaw method
calcA.maigesRaw <- function(object, bkgSub="subtract", ...) {
    
    
    ## A basic test
    if(!is.element(bkgSub, c("none", "subtract", "half", "minimum", "movingmin",
    "edwards", "normexp")))
        stop("bkgSub must be one of 'none', 'subtract', 'half', 'minimum',
        'movingmin', 'edwards' or 'normexp'")
    
    
    ## Subtracting bkg and calculating A values
    tmp <- as(object, "RGList")
    tmp2 <- backgroundCorrect(tmp, method=bkgSub, verbose=FALSE)
    A <- (log2(tmp2$R)+log2(tmp2$G))/2
    
    return(A)
    
}
