## Function to construct the design matrices to fit linear models using
## DEgenesANOVA function
##
## Parameters: data      -> object of maiges class
##             factors   -> char strings specifying the sample labels
##                          to be used as factors
##             model     -> a formula specyfing the model to be fitted
##             contrasts -> string vector specifying the contrasts to
##                          be done. This is done by the limma's
##                          makeContrasts function. Pay attention that
##                          we use the treatment control parametrization.
##             ...       -> additional parameters for model.matrix function
##
## Gustavo Esteves
## 14/05/07
##
##


designANOVA <- function(data=NULL, factors=names(data@Slabels), model=NULL,
contrasts=NULL, ...) {
    
    
    ##
    ## Doing the group-means parametrization
    ##
    
    ## Getting the factors
    ##idx <- which(names(data@Slabels) %in% factors)
    ##tmpG <- which(apply(is.na(as.matrix(data@Slabels[,idx])), 1, sum) == 0)
    
    ## Defining the coefficients of the model
    ##if(length(factors) > 1) {
    ##  tmp <- paste(factors[1], data@Slabels[, idx[1]][tmpG], sep="")
    ##  for(i in 2:length(factors)) {
    ##    tmp <- paste(tmp, paste(factors[i], data@Slabels[, idx[i]][tmpG],
    ##    sep=""), sep=".")
    ##  }
    ##}
    ##else
    ##  tmp <- paste(factors, data@Slabels[, idx][tmpG], sep="")
    
    ## Constructing the model matrix
    ##desMatrix <- model.matrix(~0+factor(tmp, levels=unique(tmp)))
    ##colnames(desMatrix) <- unique(tmp)
    ##rownames(desMatrix) <- tmpG
    
    
    
    
    
    ##
    ## Doing the treatment-contrasts parametrization
    ##
    
    ## Getting old option for contrasts and setting for treatment
    oldContOpt <- options()$contrasts
    options(contrasts=c("contr.treatment", "contr.treatment"))  
    
    
    ## Getting the factors and removing any NA value
    idx <- names(data@Slabels) %in% factors
    tmpG <- which(apply(is.na(as.matrix(data@Slabels[, idx])), 1, sum) == 0)
    data <- data[, tmpG]
    
    toModelMatrix <- as.data.frame(data@Slabels[, idx])
    if(length(factors) > 1) 
        for(i in 1:length(factors))
            toModelMatrix[, i] <- as.factor(toModelMatrix[, i])
    else
        toModelMatrix[, 1] <- as.factor(toModelMatrix[, 1])
    
    names(toModelMatrix) <- names(data@Slabels)[idx]
    
    
    ## Defining a model if it is NULL
    if(is.null(model)) {
        model <- "~"
        for(i in factors)
            model <- paste(model, "+", i)
    }
    
    
    ## Constructing design matrix
    desMatrix <- model.matrix(eval(parse(text=model)), toModelMatrix, ...)
    tmp <- colnames(desMatrix)
    tmp <- gsub("[(]", "", tmp)
    tmp <- gsub("[)]", "", tmp)
    tmp <- gsub("[:]", ".", tmp)
    colnames(desMatrix) <- tmp
    
    
    ## setting contrast option to the old option again
    options(contrasts=oldContOpt)  
    
    
    
    
    ##
    ## Contructing the contrasts matrix
    ##
    if(is.null(contrasts))
        contrasts <- colnames(desMatrix)[-1]
    
    contMatrix <- makeContrasts(contrasts=contrasts, levels=desMatrix)
    
    
    ## converting toModelMatrix to character elements
    for(i in 1:length(toModelMatrix))
        toModelMatrix[[i]] <- as.character(toModelMatrix[[i]])
    
    
    ##
    ## Contructing the maigesANOVA object to return
    ##
    
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
    
    
    ## Generating an object of type maigesANOVA
    res <- new("maigesANOVA", W=data@W, A=data@A, SD=data@SD, IC1=data@IC1,
    IC2=data@IC2, BadSpots=data@BadSpots, Layout=data@Layout,
    GeneGrps=data@GeneGrps, Paths=data@Paths, Glabels=data@Glabels,
    Slabels=toModelMatrix, Notes=data@Notes, Dmatrix=desMatrix,
    Cmatrix=contMatrix, Date=date(), V.info=vInfo)
    
    return(res)
    
}
