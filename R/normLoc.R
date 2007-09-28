## Function to normalize the microarray data
##
## Parameters: obj -> object of type maigesRaw to be normalized
##             ... -> additional parameters for function normalizeWithinArrays
##
## Gustavo H. Esteves
## 15/05/07
##
## Version: 1.0
##
##


normLoc <- function(obj=NULL, ...) {
    
    
    ## converting the obj to an object of class limma
    if(class(obj) == "maigesRaw")
        toNorm <- as(obj, "RGList")
    else
        stop("The 'obj' object isn't of class 'maigesRaw'.")
    
    
    
    ##
    ## Defining a new object for normalized data
    ##
    norm <- new("maiges")
    norm@BadSpots <- obj@BadSpots
    norm@UseSpots <- obj@UseSpots
    norm@GeneGrps <- obj@GeneGrps
    norm@Paths <- obj@Paths
    norm@Layout <- obj@Layout
    norm@Glabels <- obj@Glabels
    norm@Slabels <- obj@Slabels
    norm@Date <- date()
    
    ## Picking R and packages version information
    tmp <- sessionInfo()
    vInfo <- list()
    vInfo$R.version <- tmp$R.version$version.string
    vInfo$BasePacks <- tmp$basePkgs
    tmp1 <- NULL
    for (i in 1:length(tmp$otherPkgs))
        tmp1 <- c(tmp1, paste(tmp$otherPkgs[[i]]$Package, "version",
        tmp$otherPkgs[[i]]$Version))
    
    vInfo$AddPacks <- tmp1
    
    norm@V.info <- vInfo
    
    
    
    
    ##
    ## Doing the normalization step
    ##
    
    ## Merging UseSpots and BadSpots
    for(i in 1:dim(obj@UseSpots)[2])
        obj@UseSpots[, i] <- obj@UseSpots[, i] & !obj@BadSpots
    
    
    ## Geting samples where Ref is labelled with ch2 to do the (-1) M multiplier
    idxMult <- tolower(getLabels(obj, "Ref")) == "red"
    
    
    ## Catching the method of normalization
    add <- list(...)
    if(is.null(add$method))
        type <- "printtiploess"
    else
        type <- add$method
    
    
    
    ## Normalization by the limma method
    ##
    ## Composite normalization (control spots) from limma. Here the
    ## parameter 'controlspots' must be specified to
    ## normalizeWithinArray function
    ##
    ## Control normalization (control spots) from limma. Here the
    ## parameter 'controlspots' must be specified to
    ## normalizeWithinArray function
    ##
    ## Robust spline normalization (control spots) from limma. Here the
    ## parameters 'df' and 'robust' must be specified to
    ## normalizeWithinArray function
    tmp <- normalizeWithinArrays(toNorm, toNorm$printer,
    weights=matrix(as.numeric(obj@UseSpots), dim(obj)[1], dim(obj)[2]), ...)
    
    tmp1 <- as.matrix(tmp$M)
    if(sum(idxMult) > 0)
        tmp1[, idxMult] <- (-1)*tmp1[, idxMult]
    
    norm@W <- tmp1
    norm@A <- as.matrix(tmp$A)
    norm@Notes <- paste(obj@Notes, ", normalized by ", type, " (from limma)", sep="")
    
    
    ## Returning the normalized object
    return(norm)
    
}
