## Function to scale the microarray data (by limma function)
##
## Parameters: obj -> object of type maigesRaw to be normalized
##             ... -> additional parameters for functions normalizeBetweenArrays
##
## Gustavo H. Esteves
## 15/05/07
##
##


normScaleLimma <- function(obj=NULL, ...) {
    
    
    ## converting the obj to an object of class limma
    if(class(obj) == "maigesRaw")
        toNorm <- as(obj, "RGList")
    else if(class(obj) == "maiges")
        toNorm <- as(obj, "MAList")
    else
        stop("The 'obj' object isn't of class 'maigesRaw' or 'maiges'.")
    
    
    
    ## Defining a new object for normalized data
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
    
    
    ## Geting samples where Ref is labelled with ch2 to do the (-1) M multiplier
    idxMult <- tolower(getLabels(obj, "Ref")) == "red"
    
    
    ## Catching the method of scale normalization and
    ## doing normalization by limma's method
    add <- list(...)
    if(is.null(add$method)) {
        type <- "Aquantile"
        tmp <- normalizeBetweenArrays(toNorm, ...)
    }
    else if(add$method == "vsn") {
        type <- add$method
        tmp <- normalizeVSN(toNorm)
    }
    else {
        type <- add$method
        tmp <- normalizeBetweenArrays(toNorm, ...)
    }
    
    
    
    ## Normalization by the limma method
    ##tmp <- normalizeBetweenArrays(toNorm, ...)
    tmp1 <- as.matrix(tmp$M)
    if (sum(idxMult) > 0)
        tmp1[, idxMult] <- (-1)*tmp1[, idxMult]
    norm@W <- tmp1
    norm@A <- as.matrix(tmp$A)
    norm@Notes <- paste(obj@Notes, ", scale adjust by ", type, " (from limma)",
    sep="")
    
    
    ## Returning the object
    return(norm)
    
}
