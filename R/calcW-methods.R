## Function to calculate W values for microarray data
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
## Version: 1.0
##


## Defining a default method
calcW.default <- function(object, ...)
    return(object@W)


## Defining a maigesRaw method
calcW.maigesRaw <- function(object, bkgSub="subtract", ...) {
    
    ## A basic test
    if(!is.element(bkgSub, c("none", "subtract", "half", "minimum", "movingmin",
    "edwards", "normexp")))
        stop(" bkgSub must be one of 'none', 'subtract', 'half', 'minimum',
        'movingmin', 'edwards' or 'normexp'")
    
    ## Geting reference samples labeled with red channel
    idx <- (tolower(object@Rdye) == "red")
    
    ## Subtracting bkg
    tmp1 <- backgroundCorrect(as(object, "RGList"), method=bkgSub,
    verbose=FALSE)
    
    ## Calculating W values
    tmpM <- as.matrix(log2(tmp1$R))-as.matrix(log2(tmp1$G))
    tmpM[, idx] <- (-1)*tmpM[, idx]
    W <- tmpM
    
    return(W)
    
}
