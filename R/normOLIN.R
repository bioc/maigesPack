## Function to normalize the microarray data by OLIN method
##
## Parameters: obj -> object of type maigesRaw to be normalized
##             ... -> additional parameters for function olin
##
## Gustavo H. Esteves
## 27/05/07
##
##


normOLIN <- function(obj=NULL, ...) {
    

    require(OLIN)
  
    ## converting the obj to an object of class limma
    if(class(obj) == "maigesRaw")
        toNorm <- as(obj, "marrayRaw")
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
    if(is.null(add$OSLIN))
        type <- "OLIN"
    else if(add$OSLIN)
        type <- "OSLIN"
    else 
        type <- "OLIN"
    
    
    
    ## Normalization by methods olin ou oslin from OLIN package
    tmp <- OLIN::olin(toNorm, weights=matrix(as.numeric(obj@UseSpots), dim(obj)[1],
    dim(obj)[2]), ...)
    
    tmp1 <- as.matrix(maM(tmp))
    if(sum(idxMult) > 0)
        tmp1[, idxMult] <- (-1)*tmp1[, idxMult]
    
    norm@W <- tmp1
    norm@A <- unname(maA(tmp))
    norm@Notes <- paste(obj@Notes, ", normalized by ", type, " (from OLIN)", sep="")
    
    
    ## Returning the object
    return(norm)
    
}
