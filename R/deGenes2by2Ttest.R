## Function to do statistical hipothesis test between 2 biological types by t
## test
##
## Parameters: data      -> object of class maiges
##             sLabelID  -> sample label to be used
##             sTypeComp -> list with char vectors specifiing
##                          the sample types to be compared
##             doClust   -> Logical. The object generated from this
##                          analysis will be used for cluster
##                          analysis? Defaults to TRUE
##             ...       -> additional parameters for function t.test
##
## Gustavo Esteves
## 14/05/07
##
## Version: 1.0
##


deGenes2by2Ttest <- function(data=NULL, sLabelID=names(data@Slabels)[1],
sTypeComp=NULL, doClust=TRUE, ...) {
    
    
    ## Making some basic tests
    if(is.null(data))
        stop("The data object must be specified (class maiges).")
    if(length(sTypeComp) > 2)
        stop("You must specify only 2 groups in sTypeComp.")
    if(length(unique(getLabels(data, sLabelID))) > 2 & is.null(sTypeComp))
        warning("Factor with more than 2 levels comparing the first 2
        by default.")
    
    
    ## Getting all sample types to construct the indexes
    tSamples <- getLabels(data, sLabelID)
    
    
    ## Defining the first 2 samples as default to be compared if sTypeComp
    ## is NULL
    if(is.null(sTypeComp)) {
        idx <- !is.na(unique(tSamples))
        grps <- c(paste(sLabelID, unique(tSamples)[idx][1], sep="."),
        paste(sLabelID, unique(tSamples)[idx][2], sep="."))
        
        sTypeComp <- list(unique(tSamples)[idx][1], unique(tSamples)[idx][2])
        names(sTypeComp) <- grps
    }
    
    
    ## Putting names in the sTypeComp if it's NULL
    if(is.null(names(sTypeComp)))
        names(sTypeComp) <- c("Grp.1", "Grp.2")
    
    
    ## Locate the indexes of the tissues
    indexes <- list()
    namesComp <- names(sTypeComp)
    for(i in 1:2) {
        tmp <- tSamples %in% sTypeComp[[i]]
        tSamples[tmp] <- namesComp[i]
        indexes[[namesComp[i]]] <- tmp
    }
    
    
    ## Defining the table to be used to do the tests
    wTable <- calcW(data)
    
    
    ## Defining a maigesDE class object
    if(doClust)
        dataComparisons <- new("maigesDEcluster")
    else
        dataComparisons <- new("maigesDE")
    
    resultFold <- NULL
    resultStat <- NULL
    resultP <- NULL
    
    ## Geting information about genes and samples...
    ## Pasting "(*)" at end of bad spots
    dataComparisons@GeneInfo <- data@Glabels
    for(i in 1:length(dataComparisons@GeneInfo))
        dataComparisons@GeneInfo[[i]][data@BadSpots] <-
        paste(dataComparisons@GeneInfo[[i]][data@BadSpots], "(*)")
    
    
    ## Getting information about the samples tested
    samp <- as.data.frame(getLabels(data, sLabelID))
    names(samp) <- sLabelID
    dataComparisons@SampleInfo <- samp
    if(doClust)
        dataComparisons@W <- wTable
    comparison <- NULL
    
    
    ## Defining the comparison
    comparison <- c(comparison, paste(namesComp[1], "-", namesComp[2], sep=""))
    sampleNumerator <- indexes[[namesComp[1]]]
    sampleDenumerator <- indexes[[namesComp[2]]]
    
    
    ## If only one gene, convert the vector to matrix
    if(is.vector(wTable))
        wTable <- t(as.matrix(wTable))
    
    
    ## Applying t-test
    tmpBasicTests <- t(apply(wTable, 1, function(x, obs1, obs2, ...) {
    tmpOut <- t.test(x[obs1], x[obs2], ...)
    meanDif <- mean(x[obs1], na.rm=TRUE)-mean(x[obs2], na.rm=TRUE)
    stat <- c(meanDif=meanDif, Statistic=unname(tmpOut$statistic),
    P.value=unname(tmpOut$p.value))
    return(stat) }, obs1=sampleNumerator, obs2=sampleDenumerator, ...))
    
    
    
    ## Constructing object to return
    resultFold <- cbind(resultFold, tmpBasicTests[, 1])
    resultStat <- cbind(resultStat, tmpBasicTests[, 2])
    resultP <- cbind(resultP, tmpBasicTests[, 3])
    
    colnames(resultFold) <- comparison
    colnames(resultStat) <- comparison
    colnames(resultP) <- comparison
    
    
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
    
    
    ##dataComparisons@testResults <- result
    dataComparisons@fold <- resultFold
    dataComparisons@stat <- resultStat
    dataComparisons@p.value <- resultP
    dataComparisons@test <- "T test"
    dataComparisons@factors <- sLabelID
    dataComparisons@Date <- date()
    dataComparisons@V.info <- vInfo
    
    return(dataComparisons)
    
}
