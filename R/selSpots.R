## Function to filter a marray object to apply to normalization step.
##
## Parameters: obj       -> object of type maigesRaw
##             sigNoise  -> real number indicating the cuttoff to remove spots
##                          below it
##             rmFlag    -> vector of simbols of Flags to be removed (for
##                          normalization)
##             gLabelsID -> character vector indicating the gene labels to be
##                          searched
##             remove    -> list of same length as GlabelsID with character
##                          vector indicating spots to be removed
##             badSpots  -> index of bad spots (numeric or logic), identifying
##                          bad spots. May be the names (from badLabel)
##             badLabel  -> char string specifying the label for badSpots
##
##
## Gustavo H. Esteves
## 27/05/07
##
##


selSpots <- function(obj=NULL, sigNoise=1, rmFlag=NULL, gLabelsID=c("Name"),
remove=list(c("BLANK","DAP","LYS","PHE","Q_GENE","THR","TRP")), badSpots=NULL,
badLabel=NULL) {
    
    
    ## atualizing BadSpots slot
    if(!is.null(badSpots)) {
        if(is.character(badSpots)) {
            idxTmp <- getLabels(obj, badLabel, FALSE) %in% badSpots
            obj@BadSpots <- obj@BadSpots | idxTmp
        }
        else if (is.numeric(badSpots)) {
            tmp <- rep(FALSE, dim(obj)[1])
            tmp[badSpots] <- TRUE
            obj@BadSpots <- obj@BadSpots | badSpots
        }
        else if (is.logical(badSpots))
            obj@BadSpots <- obj@BadSpots | badSpots
        else
            stop("The argument 'badSpots' must be logical, numerical
            or character.")
        
    }
    
    
    ## removing spots according to gLabelsID and remove parameters
    tmp <- NULL
    if(!is.null(gLabelsID)) {
        for (i in 1:length(gLabelsID)) {
            labTmp <- getLabels(obj, gLabelsID[i], FALSE)
            tmp <- c(tmp, which(labTmp %in% remove[[i]]))
        }
        obj@UseSpots[tmp, ] <- FALSE
    }
    
    
    ## removing spots with signal/back <= sigNoise
    if(!is.null(sigNoise)) 
        for (i in 1:dim(obj)[2]) {
            sig2noise <- apply(cbind(obj@Sf[, i]/obj@Sb[, i],
            obj@Rf[, i]/obj@Rb[, i]), 1, min)
            
            tmp <- sig2noise <= sigNoise
            obj@UseSpots[tmp, i] <- FALSE
        }
    
    
    ## removing flags
    if(!is.null(rmFlag)) 
        for (i in 1:dim(obj)[2]) {
            tmp <- obj@Flag[, i] %in% rmFlag
            obj@UseSpots[tmp, i] <- FALSE
        }
    
    obj@Date <- date()
    
    
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
    
    obj@V.info <- vInfo
    
    
    ## return maigesRaw object with UseSpots slot atualized
    return(obj)
    
}
