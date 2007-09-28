## Function to construct relevance networks for 2 biological types (it is the
## old RelNetwork function)
##
##
## Parameters: data     -> Object of class maiges
##             gLabelID -> Identification of gene name label ID.
##             sLabelID -> sample label identification to be used
##             geneGrp  -> character string (or numeric index) specifying the
##                         gene group where to calculate the correlation
##                         values between them. If NULL (together with path)
##                         all genes are used.
##             path     -> character string (or numeric index) specifying the
##                         pathway where to calculate the correlation values
##                         between them. If NULL (together with geneGrp) all
##                         genes are used.
##             samples  -> a named list with two character vectors specifying
##                         the two groups that must be compared
##             type     -> type of correlation to be calculated. May be
##                         'Rpearson' (default), 'pearson', 'kendall',
##                         'spearman' or 'MI'
##             ...      -> additional parameters for functions robustCorr or cor
##
##
##
## Gustavo H. Esteves
## 15/05/07
##
## Version: 1.0
##


relNetworkM <- function(data=NULL, gLabelID="GeneName",
sLabelID="Classification", geneGrp=NULL, path=NULL, samples=NULL,
type="Rpearson", ...) {
    
    
    ## Doing a simple test...
    if(is.null(data))
        stop("The data MUST be specified!!!")
    
    ##if(!(type %in% c("Rpearson", "pearson", "kendall", "spearman", "MI")))
    ##    stop("Parameter 'type' must be 'Rpearson', 'pearson', 'kendall',
    ##    'spearman' or 'MI'.")
    if(!(type %in% c("Rpearson", "pearson", "kendall", "spearman")))
        stop("Parameter 'type' must be 'Rpearson', 'pearson', 'kendall' or
        'spearman'.")
    
    if(!is.null(geneGrp) & !is.null(path))
        stop("You must specify only one of geneGrp and path, or leave
        both NULL.")
    
    
    ## Getting gene and sample labels
    allGenes <- getLabels(data, gLabelID, FALSE)
    allSamples <- getLabels(data, sLabelID)
    
    
    ## removing genes and samples not named or NA from the data object
    idxGtmp <- which(allGenes != "" & !is.na(allGenes))
    idxStmp <- which(allSamples != "" & !is.na(allSamples))
    data <- data[idxGtmp, idxStmp]
    
    
    ## atualizing all gene and sample labels
    allGenes <- getLabels(data, gLabelID, FALSE)
    allGenes[data@BadSpots] <- paste(allGenes[data@BadSpots], "(*)")
    allSamples <- getLabels(data, sLabelID)
    
    
    ## Constructing the table to be used
    table <- calcW(data)
    if(is.null(geneGrp) & is.null(path))
        genes <- allGenes
    else if(!is.null(geneGrp)) {
        if(!is.numeric(geneGrp)) {
            geneGrp <- which(colnames(data@GeneGrps) == geneGrp)
            genes <- allGenes[data@GeneGrps[, geneGrp]]
        }
        else
            genes <- allGenes[data@GeneGrps[, geneGrp]]
    }
    else if(!is.null(path)) {
        genePath <- which(nodes(data@Paths[[path]]) %in% getLabels(data,
        data@Paths[[1]], FALSE))
        
        genes <- allGenes[genePath]
    }
    
    
    ## Testing the length of genes
    if (length(genes) < 2)
        stop("There is less than 2 elements in this group of genes.")
    
    
    ## Getting the first 2 sample types to use in case of NULL samples
    ## parameter
    if(is.null(samples)) {
        samples <- list(unique(allSamples)[1], unique(allSamples)[2])
        names(samples) <- unique(allSamples)[1:2]
    }
    
    
    ## Finding indexes of the samples specified and respective lengths
    idxTmp <- allSamples %in% c(samples[[1]], samples[[2]])
    table <- table[, idxTmp]
    sTypes <- allSamples[idxTmp]
    
    
    ## Constructing the table to be used, averaging eventual replicates...
    tableTmp <- NULL
    for(i in unique(genes)) {
        idx <- allGenes %in% i
        if(sum(idx) > 1)
            tableTmp <- rbind(tableTmp, apply(table[idx, ], 2, mean,
            na.rm=TRUE))
        
        else
            tableTmp <- rbind(tableTmp, table[idx, ])
    }
    
    
    ## Naming tableTmp
    rownames(tableTmp) <- unique(genes)
    colnames(tableTmp) <- sTypes
    colnames(tableTmp)[sTypes %in% samples[[1]]] <- names(samples)[1]
    colnames(tableTmp)[sTypes %in% samples[[2]]] <- names(samples)[2]
    
    
    ## Finding indexes of the samples specified and respective lengths
    sample1 <- which(sTypes %in% samples[[1]])
    sample2 <- which(sTypes %in% samples[[2]])
    n1 <- length(sample1)
    n2 <- length(sample2)
    
    
    ## Calculating the correlation coefficients
    if(type == "Rpearson") {
        corObj1 <- robustCorr(tableTmp[, sample1], ...)[[1]]
        corObj2 <- robustCorr(tableTmp[, sample2], ...)[[1]]
        ## Calculating the p-values of the differences
        difPvalue <- compCorr((n1-1), corObj1, (n2-1), corObj2)[[2]]
    }
    ##else if(type == "MI") {
    ##    corObj1 <- MI(tableTmp[, sample1], ...)
    ##    corObj2 <- MI(tableTmp[, sample2], ...)
    ## Calculating the p-values of the differences
    ##    difPvalue <- compCorr(n1, corObj1, n2, corObj2)[[2]]
    ##}
    else {
        corObj1 <- cor(t(tableTmp[, sample1]), method=type, ...)
        corObj2 <- cor(t(tableTmp[, sample2]), method=type, ...)
        ## Calculating the p-values of the differences
        difPvalue <- compCorr(n1, corObj1, n2, corObj2)[[2]]
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
    
    
    ## Defining the object to return
    result <- new("maigesRelNetM", W=tableTmp, Corr1=corObj1, Corr2=corObj2,
    DifP=difPvalue, Date=date(), types=names(samples), Slabel=sLabelID,
    V.info=vInfo)
    
    return(result)
    
}
