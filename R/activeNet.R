## Function to search by gene networks with altered profiles (functional
## classification of gene networks).
##
## Parameters: data        -> object of maiges class
##             sLabelID    -> Label identification to be used
##             samples     -> a list with character vectors specifying
##                            the groups that must be compared
##             type        -> type of correlation to be calculated. May
##                            be 'Rpearson' (default), 'pearson',
##                            'kendall', 'spearman'.
##             bRep        -> number of bootstraps to be done in the
##                            correlation test
##             alternative -> type of test. May be 'less' or 'greater'.
##                            less test test the activity of the
##                            original network, greater test the
##                            activity of the 'reversed' network
##             adjP        -> type of p-value adjustment. May be
##                            "Bonferroni", "Holm", "Hochberg",
##                            "SidakSS", "SidakSD", "BH",
##                            "BY" or "none". Defaults to "none"
##
##
##
##
## Gustavo H. Esteves
## 27/05/07
##
## Version: 1.1
##


activeNet <- function(data=NULL, samples=NULL, sLabelID="Classification",
type="Rpearson", bRep=1000, alternative="greater", adjP="none") {
    
    
    ## Making some basic initial tests...
    if(is.null(data))
        stop("You MUST specify a valid 'data' object.")
    
    if(!(sLabelID %in% names(data@Slabels)))
        stop(paste(sLabelID, "is't a valid sample label."))
    
    if(!(type %in% c("Rpearson", "pearson", "kendall", "spearman", "MI")))
        stop("Parameter 'type' must be 'Rpearson', 'pearson', 'kendall',
        'spearman' or 'MI'.")
    
    if(!(alternative %in% c("less", "greater")))
        stop("Parameter 'alternative' must be 'less' or 'greater'.")  
    
    if(!is.element(adjP, c("Bonferroni", "Holm", "Hochberg", "SidakSS",
    "SidakSD", "BH", "BY", "none")))
        stop("Incorrect value for parameter adjP see help pages for or
        'mt.rawp2adjp'")
    
    
    ## Getting gene and sample labels
    if(is.null(samples)) {
        samples <- getLabels(data, sLabelID)
        idx <- !is.na(samples)
        samples <- samples[idx]
        data <- data[, idx]
    }
    else {
        tmp <- getLabels(data, sLabelID)
        idx <- !is.na(tmp)
        tmp <- tmp[idx]
        data <- data[, idx]
        for (i in 1:length(samples)) {
            idxTmp <- is.element(tmp, samples[[i]])
            tmp[idxTmp] <- names(samples)[i]
        }
        samples <- tmp
    }
    
    
    ## Removing samples with less than 2 observations
    conditions <- unique(samples)
    for(i in conditions) {
        idx <- which(samples == i)
        if (length(idx) < 2) {
            warning(paste("Sample", i,
            "was removed because it presents less than 2 elements."))
            
            samples <- samples[-idx]
        }
    }
    conditions <- unique(samples)
    
    
    ## Geting colnames of GeneGrps slot as modules
    paths <- names(data@Paths)[-1]
    
    
    ## Defining the table to be used
    table <- calcW(data)
    genes <- getLabels(data, data@Paths[[1]], FALSE)
    genes[data@BadSpots] <- paste(genes[data@BadSpots], "(*)")
    
    nGenes <- length(table[, 1]) ## Number of Genes
    nCond <- length(conditions) ## Number of Arrays
    nPaths <- length(paths) ## Number of modules (gene groups)
    
    
    ## Defining a matrix with the scores...
    pathScores <- matrix(NA, nrow=nCond, ncol=nPaths)
    pathPvalues <- matrix(NA, nrow=nCond, ncol=nPaths)
    for (j in 1:nPaths) {
        for (i in 1:nCond) {
            scoreTmp <- NULL
            sampIdx <- samples %in% conditions[i]
            for(k in nodes(data@Paths[[j+1]])) {
                tmp1 <- edges(data@Paths[[j+1]], k)[[1]]
                tmp2 <- edgeWeights(data@Paths[[j+1]], k)[[1]]
                if(length(tmp1) > 0) {
                    for(node in 1:length(tmp1)) {
                        if(tmp2[node] > 0)
                            typeTest <- "greater"
                        if(tmp2[node] < 0)
                            typeTest <- "less"
                        if(tmp2[node] == 0)
                            stop("type of iteration equal zero.")
                        idx1 <- getLabels(data, data@Paths[[1]], FALSE) %in% k
                        idx2 <- getLabels(data, data@Paths[[1]], FALSE) %in%
                        tmp1[node]
                        
                        if(type == "MI") {
                            pValue <- bootstrapMI(table[idx1, sampIdx],
                            table[idx2, sampIdx], bRep=bRep)
                        }
                        else {
                            pValue <- bootstrapCor(table[idx1, sampIdx],
                            table[idx2, sampIdx], type=type, bRep=bRep,
                            alternative=typeTest)
                        }
                        
                        if(pValue == 0)
                            pValue <- 1/(bRep+1)
                        
                        scoreTmp <- c(scoreTmp, pValue)
                        
                    }
                }
            }
            
            ## Adjusting p-values if specified by the user.
            if(adjP != "none") {
                tmp <- multtest::mt.rawp2adjp(scoreTmp, proc=adjP)
                scoreTmp <- tmp$adjp[order(tmp$index), 2]
            }
            score <- sum(-log(scoreTmp))/(length(scoreTmp))
            pathScores[i, j] <- score
            if(alternative == "greater")
                pathPvalues[i, j] <- pgamma(score, shape=length(scoreTmp),
                scale=(1/length(scoreTmp)), lower.tail=FALSE)
            
            else
                pathPvalues[i, j] <- pgamma(score, shape=length(scoreTmp),
                scale=(1/length(scoreTmp)), lower.tail=TRUE)
            
        }
    }
    rownames(pathScores) <- conditions
    colnames(pathScores) <- paths
    rownames(pathPvalues) <- conditions
    colnames(pathPvalues) <- paths
    
    
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
    return(new("maigesActNet", scores=pathScores, Pvalues=pathPvalues,
    Date=date(), V.info=vInfo))
    
}
