## Function to do bootstrap lowess normalization...
##
## parameters: raw       -> maigesRaw object to be normalized
##             span      -> span for loess regression
##             propLoess -> proportion of spots (after filtering) to be used in
##                          each iteration of loess
##             nRep      -> number of repetitions for loess procedure
##             func      -> function to estimate the final W value. You must use
##                          'mean' or 'median' or 'none' (default)
##             bkgSub    -> bkg subtraction method, using backgroundCorrection
##                          from limma
##             ...       -> additional parameters for 'loessFit' function from
##                          limma package
##
## Gustavo H. Esteves
## 15/05/07
##
##


normRepLoess <- function(raw, span=0.4, propLoess=0.5, nRep=50, func="none",
bkgSub="none", ...) {
    
    
    ## Doing some basic tests
    if(class(raw) != "maigesRaw")
        stop("The 'raw' object isn't of class 'maigesRaw'.")
    
    if(!is.element(bkgSub, c("none","subtract","half","minimum","movingmin",
    "edwards","normexp")))
        stop("The bkgSub parameter must be 'none', 'subtract', 'half',
        'minimum', 'movingmin', 'edwards', 'normexp'.")
    
    if(!is.element(func, c("none","mean","median")))
        stop("The func parameter must be 'none','mean','median'.")
    
    if(nRep < 2)
        stop("nRep must be greater than 1.")
    
    
    
    ## Calculating W and A values
    W <- calcW(raw, bkgSub)
    A <- calcA(raw, bkgSub)
    
    
    ## Defining objects to normalized values and SDs and ICs
    wNorm <- NULL
    SDNorm <- NULL
    ICnorm1 <- NULL
    ICnorm2 <- NULL
    for (i in 1:ncol(W)) {
        wTmp <- W[, i]
        aTmp <- A[, i]
        weight <- as.matrix(as.numeric(raw@UseSpots[, i] & !raw@BadSpots))
        
        sig2noise <- apply(cbind(raw@Sf[, i]/raw@Sb[, i],
        raw@Rf[, i]/raw@Rb[, i]), 1, min)
        
        
        ## Doing a qualitative test
        if(unname(quantile(sig2noise, 0.5)) < 1)
            cat("\tAttention!!! More than 50% of spots with Sig < Bkg in two
            channels!\n")
        
        idxTmp <- which(weight > 0) ## index of good genes
        rmGenes <- round(length(idxTmp)*(1-propLoess)) ## n. of genes to sample
        wNormTmp <- NULL
        ## Doing loess repetition
        for (j in 1:nRep) { 
            weightTmp <- weight
            weightTmp[sample(idxTmp, rmGenes)] <- 0

            wNormTmp <- cbind(wNormTmp, unname(loessFit(wTmp, aTmp,
            weights=weightTmp, span=span, ...)$residuals))
            
        }
        
        ## Calculating the value to used as normalized
        if(func == "none") ## Global normalization
            wNorm <- cbind(wNorm, unname(loessFit(wTmp, aTmp, weights=weight,
            span=span, ...)$residuals))
        
        else ## Picking mean or median of the booted values
            wNorm <- cbind(wNorm, apply(wNormTmp, 1, eval(parse(text=func)),
            na.rm=TRUE))
        
        ## Calculating SDs and ICs for the normalized values
        SDNorm <- cbind(SDNorm, apply(wNormTmp, 1, sd, na.rm=TRUE))

        ICnorm1 <- cbind(ICnorm1, apply(wNormTmp, 1, quantile, probs=0.05,
        na.rm=TRUE))
        
        ICnorm2 <- cbind(ICnorm2, apply(wNormTmp, 1, quantile, probs=0.95,
        na.rm=TRUE))
        
    }
    
    
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
    
    
    ## Returning the object
    return(new("maiges", W=wNorm, A=A, SD=SDNorm, IC1=ICnorm1, IC2=ICnorm2,
    Layout=raw@Layout, GeneGrps=raw@GeneGrps, Paths=raw@Paths,
    Glabels=raw@Glabels, BadSpots=raw@BadSpots, UseSpots=raw@UseSpots,
    Slabels=raw@Slabels, Notes=raw@Notes, Date=date(), V.info=vInfo))
    
}
