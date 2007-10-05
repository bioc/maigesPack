## Function to summarize measures of a maiges class object (both sample and
## genes)
##
## Parameters: object    -> an object of class maiges
##             gLabelID  -> Gene id to be used to summarize
##             sLabelID  -> Sample id to be used to summarize
##             funcS     -> function to be applied for sample replicates
##             funcG     -> function to be applied for gene replicates
##             rmBad     -> remove bad spots?
##             keepEmpty -> logical indicating if you want to maintain spots
##                           with empty information
##
##
## Gustavo Esteves
## 15/05/07
##
##


summarizeReplicates <- function(object=NULL, gLabelID="GeneName",
sLabelID="Sample", funcS="mean", funcG="mean", rmBad=TRUE, keepEmpty=TRUE) {
    
    
    ## Doing basic tests
    if (is.null(object))
        stop("You must specify an object of class maiges.")
    
    if (is.null(gLabelID))
        stop("You must specify a string in gLabelID.")
    
    if (is.null(sLabelID))
        stop("You must specify a string in sLabelID.")
    
    if(!is.null(funcS))
        if(!(funcS %in% c("mean", "median")))
            stop("funcS parameter must be 'mean' or 'median'.")
    
    if(!is.null(funcG))
        if(!(funcG %in% c("mean", "median")))
            stop("funcG parameter must be 'mean' or 'median'.")
    
    
    
    ## Removing bad spots
    if(rmBad)
        object <- object[!object@BadSpots, ]
    
    
    ## Geting Objects
    if(!keepEmpty) {
        idxCol <- getLabels(object, sLabelID) != ""
        idxRow <- getLabels(object, gLabelID, FALSE) != ""
        object <- object[idxRow, idxCol]
    }
    
    wObj <- calcW(object)
    aObj <- calcA(object)
    sd <- object@SD
    groups <- object@GeneGrps
    interval1 <- object@IC1
    interval2 <- object@IC2
    sampLab <- object@Slabels
    geneLab <- object@Glabels
    badSpots <- object@BadSpots
    useSpots <- object@UseSpots
    
    
    ## Resuming data replicates if this was specified (for samples)
    wTmp <- NULL
    aTmp <- NULL
    if(sum(dim(wObj) == dim(sd)) == 2)
        sdTmp <- NULL
    else
        sdTmp <- sd
    if(sum(dim(wObj) == dim(interval1)) == 2)
        ic1Tmp <- NULL
    else
        ic1Tmp <- interval1
    if(sum(dim(wObj) == dim(interval2)) == 2)
        ic2Tmp <- NULL
    else
        ic2Tmp <- interval2
    if(sum(dim(wObj) == dim(useSpots)) == 2)
        useSpotsTmp <- NULL
    else
        useSpotsTmp <- useSpots
    
    sampTmp <- data.frame()
    geneTmp <- data.frame()
    processed <- NULL
    if(!is.null(funcS)) {
        sampleLabels <- getLabels(object, sLabelID)
        for(i in 1:length(sampleLabels)) {
            if(sampleLabels[i] == "") {
                wTmp <- cbind(wTmp, wObj[, i])
                aTmp <- cbind(aTmp, aObj[, i])
                if(sum(dim(wObj) == dim(sd)) == 2)
                    sdTmp <- cbind(sdTmp, sd[, i])
                if(sum(dim(wObj) == dim(interval1)) == 2) {
                    ic1Tmp <- cbind(ic1Tmp, interval1[, i])
                    ic2Tmp <- cbind(ic2Tmp, interval2[, i])
                }
                if(sum(dim(wObj) == dim(useSpots)) == 2)
                    useSpotsTmp <- cbind(useSpotsTmp, useSpots[, i])
                sampTmp <- rbind(sampTmp, sampLab[i, ])
            }
            else if (!(sampleLabels[i] %in% processed)) {
                idx <- which(sampleLabels %in% sampleLabels[i])
                
                wTmp <- cbind(wTmp, apply(as.matrix(wObj[, idx]), 1,
                eval(parse(text=funcS)), na.rm=TRUE))
                
                aTmp <- cbind(aTmp, apply(as.matrix(aObj[, idx]), 1,
                eval(parse(text=funcS)), na.rm=TRUE))
                
                if(sum(dim(wObj) == dim(sd)) == 2)
                    sdTmp <- cbind(sdTmp, apply(as.matrix(sd[, idx]), 1, mean,
                    na.rm=TRUE))
                
                if(sum(dim(wObj) == dim(interval1)) == 2) {
                    ic1Tmp <- cbind(ic1Tmp,
                    apply(as.matrix(interval1[, idx]), 1, min, na.rm=TRUE))
                    
                    ic2Tmp <- cbind(ic2Tmp,
                    apply(as.matrix(interval2[, idx]), 1, max, na.rm=TRUE))
                    
                }
                if(sum(dim(wObj) == dim(useSpots)) == 2)
                    useSpotsTmp <- cbind(useSpotsTmp,
                    apply(as.matrix(useSpots[, idx]), 1, sum) == length(idx))
                
                sampTmp <- rbind(sampTmp, sampLab[idx[1], ])
                processed <- c(processed, sampleLabels[i])
            }
        }
        wObj <- wTmp
        aObj <- aTmp
        sd <- sdTmp
        interval1 <- ic1Tmp
        interval2 <- ic2Tmp
        useSpots <- useSpotsTmp
        sampLab <- sampTmp
        
        if(!is.null(funcG)) {
            ## Re-creating the NULL auxiliary objects
            wTmp <- NULL
            aTmp <- NULL
            if(sum(dim(wObj) == dim(sd)) == 2)
                sdTmp <- NULL
            else
                sdTmp <- sd
            if(sum(dim(wObj) == dim(interval1)) == 2)
                ic1Tmp <- NULL
            else
                ic1Tmp <- interval1
            if(sum(dim(wObj) == dim(interval2)) == 2)
                ic2Tmp <- NULL
            else
                ic2Tmp <- interval2
            if(sum(dim(wObj) == dim(useSpots)) == 2)
                useSpotsTmp <- NULL
            else
                useSpotsTmp <- useSpots
            
            geneTmp <- data.frame()
            
        }
    }
    
    ## Resuming data replicates if this was specified (for genes)
    groupsTmp <- NULL
    processed <- NULL
    bSpots <- logical()
    if(!is.null(funcG)) {
        geneLabels <- getLabels(object, gLabelID, FALSE)
        for(i in 1:length(geneLabels)) {
            if(geneLabels[i] == "") {
                wTmp <- rbind(wTmp, wObj[i, ])
                aTmp <- rbind(aTmp, aObj[i, ])
                if(sum(dim(object@GeneGrps)) != 0)
                    groupsTmp <- rbind(groupsTmp, groups[i, ])
                
                bSpots <- c(bSpots, badSpots[i])
                
                if(sum(dim(wObj) == dim(sd)) == 2)
                    sdTmp <- rbind(sdTmp, sd[i, ])
                if(sum(dim(wObj) == dim(interval1)) == 2) {
                    ic1Tmp <- rbind(ic1Tmp, interval1[i, ])
                    ic2Tmp <- rbind(ic2Tmp, interval2[i, ])
                }
                
                if(sum(dim(wObj) == dim(useSpots)) == 2)
                    useSpotsTmp <- rbind(useSpotsTmp, useSpots[i, ])

                geneTmp <- rbind(geneTmp, geneLab[i, ])
                
            }
            else if (!(geneLabels[i] %in% processed)) {
                
                idx <- which(geneLabels %in% geneLabels[i])
                
                if(length(idx) > 1) {
                    
                    cat(paste("\nGene", geneLabels[i], "with", length(idx),
                    "copies. Indexes:\n   "))
                    
                    cat(idx)
                    cat("\n")
                    
                    wTmp <- rbind(wTmp, apply(as.matrix(wObj[idx, ]), 2,
                    eval(parse(text=funcG)), na.rm=TRUE))
                    
                    aTmp <- rbind(aTmp, apply(as.matrix(aObj[idx, ]), 2,
                    eval(parse(text=funcG)), na.rm=TRUE))
                    
                    if(sum(dim(object@GeneGrps)) != 0)
                        groupsTmp <- rbind(groupsTmp,
                        as.logical(apply(groups[idx, ], 2, sum)))
                    
                    bSpots <- c(bSpots, as.logical(sum(badSpots[idx])))
                    
                    if(sum(dim(wObj) == dim(sd)) == 2)
                        sdTmp <- rbind(sdTmp, apply(as.matrix(sd[idx, ]), 2,
                        mean, na.rm=TRUE))
                    
                    if(sum(dim(wObj) == dim(interval1)) == 2) {

                        ic1Tmp <- rbind(ic1Tmp,
                        apply(as.matrix(interval1[idx, ]), 2, min, na.rm=TRUE))
                        
                        ic2Tmp <- rbind(ic2Tmp,
                        apply(as.matrix(interval2[idx, ]), 2, max, na.rm=TRUE))
                        
                    }
                    if(sum(dim(wObj) == dim(useSpots)) == 2)
                        useSpotsTmp <- rbind(useSpotsTmp,
                        apply(as.matrix(useSpots[idx, ]), 2, sum) ==
                        length(idx))
                    
                }
                else {
                    wTmp <- rbind(wTmp, wObj[idx, ])
                    aTmp <- rbind(aTmp, aObj[idx, ])
                    if(sum(dim(object@GeneGrps)) != 0)
                        groupsTmp <- rbind(groupsTmp, groups[idx, ])
                    
                    bSpots <- c(bSpots, badSpots[idx])

                    if(sum(dim(wObj) == dim(sd)) == 2)
                        sdTmp <- rbind(sdTmp, sd[idx, ])
                    
                    if(sum(dim(wObj) == dim(interval1)) == 2) {
                        ic1Tmp <- rbind(ic1Tmp, interval1[idx, ])
                        ic2Tmp <- rbind(ic2Tmp, interval2[idx, ])
                    }
                    if(sum(dim(wObj) == dim(useSpots)) == 2)
                        useSpotsTmp <- rbind(useSpotsTmp, useSpots[idx, ])
                }
                geneTmp <- rbind(geneTmp, geneLab[idx[1], ])
                processed <- c(processed, geneLabels[i])
            }
        }
        
        wObj <- wTmp
        aObj <- aTmp
        sd <- sdTmp
        interval1 <- ic1Tmp
        interval2 <- ic2Tmp
        badSpots <- bSpots
        useSpots <- useSpotsTmp
        
        if(sum(dim(object@GeneGrps)) != 0)
            groups <- groupsTmp
        geneLab <- geneTmp
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
    
    
    
    ## Creating a new object
    newObj <- new("maiges", W=wObj, A=aObj, SD=sd, IC1=interval1, IC2=interval2,
    Layout=object@Layout, BadSpots=badSpots, UseSpots=useSpots, GeneGrps=groups,
    Paths=object@Paths, Glabels=geneLab, Slabels=sampLab, Notes=object@Notes,
    Date=date(), V.info=vInfo)
    
    
    ## Returning another maiges object
    return(newObj)
    
}
