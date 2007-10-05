## Generic function to do som cluster analysis
##
## Parameters: data     -> maiges (or maigesRaw) class object
##             group    -> Type of grouping: by rows ('R') or columns ('C')
##             distance -> char string giving the type of distance to use from
##                         function
##             method   -> char string specifying the linkage method for the
##                         hierarquical cluster
##             sampleT  -> list with 2 vectors. The 1st one specify sample types
##                         to be colored and the 2nd one specify the respective
##                         colors. If NULL (default) black is used to all
##             doHier   -> logical indicating if you want to do the hierarquical
##                         branch in the opositte dimension of
##                         clustering. Defaults to FALSE and is only applicable
##                         for SOM or KM
##             sLabelID -> Sample label id
##             gLabelID -> Gene label id
##             idxTest  -> index of the test to be sorted
##             adjP     -> method of p-value adjustment
##             nDEgenes -> number of DE genes to be selected. If in (0,1) all
##                         genes with p.value <= n.dif.genes will be used
##             ...      -> additional parameters for som function
##
## Gustavo H. Esteves
## 27/05/07
##
##


somMde <- function(data, group=c("C", "R")[1], distance="correlation",
method="complete", sampleT=NULL, doHier=FALSE, sLabelID="SAMPLE",
gLabelID="GeneName", idxTest=1, adjP="none", nDEgenes=0.05, ...) {
    
    
    ## Some some basic tests
    if(!(group %in% c("C", "R")))
        stop("Parameter 'group' must be 'C' or 'R'.")
    
    if(!is.element(adjP, c("Bonferroni", "Holm", "Hochberg", "SidakSS",
    "SidakSD", "BH", "BY", "none")))
        stop("Incorrect value for parameter adjP see help pages for
        or 'mt.rawp2adjp'")
    
    if(!(method %in% c("ward", "single", "complete", "average", "mcquitty",
    "median", "centroid")))
        stop("Parameter method must be 'ward', 'single', 'complete', 'average',
        'mcquitty', 'median' or 'centroid'.")
    
    if(!(distance %in% c("euclidean", "correlation")))
        stop("Parameter 'distance' must be 'euclidean' or 'correlation'.")
    
    
    ## Adjusting p-values if specified by the user.
    if(adjP != "none") {
        tmp1 <- multtest::mt.rawp2adjp(data@p.value[, idxTest], proc=adjP)
        data@p.value[, idxTest] <- tmp1$adjp[order(tmp1$index), 2]
    }
    
    
    ## Selecting the genes according to the p-values
    if(!is.null(nDEgenes)) {
        if(nDEgenes > 0 & nDEgenes < 1)
            idxTmp <- which(data@p.value[, idxTest] <= nDEgenes)
        else
            idxTmp <- sort(data@p.value[, idxTest],
            index.return=TRUE)$ix[1:nDEgenes]
        
        if(length(idxTmp) == 0)
            stop("There are not genes satisfying your criteria.")
        
        tmpMatrix <- data@W[idxTmp, ]
        colnames(tmpMatrix) <- getLabels(data, sLabelID)
        rownames(tmpMatrix) <- getLabels(data, gLabelID, FALSE)[idxTmp]
    }
    else {
        tmpMatrix <- data@W
        colnames(tmpMatrix) <- getLabels(data, sLabelID)
        rownames(tmpMatrix) <- getLabels(data, gLabelID, FALSE)
    }
    
    idx <- apply(!is.finite(tmpMatrix), 1, sum) == 0
    tmpMatrix <- tmpMatrix[idx, ]
    
    
    ## Defining limits for ploting
    ##zLim <- range(tmpMatrix)
    ##zLim <- c(-max(abs(zLim)), max(abs(zLim)))
    
    add <- list(...)
    if(sum(c("xdim", "ydim", "topol") %in% names(add)) < 3)
        stop("You must specify the parametes 'xdim', 'ydim' and 'topol' for
        som function!")
    
    if(group == "C") {
        if(distance == "correlation")
            SOM <- som::som(som::normalize(t(tmpMatrix)), ...)
        else
            SOM <- som::som(t(tmpMatrix), ...)
        
        groups <- paste(SOM$visual$x, SOM$visual$y, sep="")
        
        heatmapsM(tmpMatrix, distfun=function(c) amap::Dist(c, method=distance),
        hclustfun=function(d) hclust(d, method=method), groups, sampleT=sampleT,
        doHier=doHier, col=maigesPack:::greenRed())
        
    }
    else {
        if(distance == "correlation")
            SOM <- som::som(som::normalize(tmpMatrix), ...)
        else
            SOM <- som::som(tmpMatrix, ...)
        
        groups <- paste(SOM$visual$x, SOM$visual$y, sep="")
        
        heatmapsM(tmpMatrix, distfun=function(c) amap::Dist(c, method=distance),
        hclustfun=function(d) hclust(d, method=method), groups, sampleT=sampleT,
        doHier=doHier, col=maigesPack:::greenRed())
        
    }
    
    invisible(SOM)
    
}
