## Function to do hierarquical cluster analysis
##
## Parameters: data      -> object of class maigesDEcluster
##             group     -> Type of grouping: by rows ('R'), columns ('C') or
##                          both ('B')
##             distance  -> char string giving the type of distance to use
##                          from function Dist (lib amap)
##             method    -> char string specifying the linkage method for the
##                          hierarquical cluster
##             doHeat    -> logical indicanting to do or not the heatmap
##             sLabelID  -> Sample label id
##             gLabelID  -> Gene label id
##             idxTest   -> index of the test to be sorted
##             adjP      -> method of p-value adjustment
##             nDEgenes  -> number of DE genes to be selected. If in (0,1) all
##                          genes with p.value <= n.dif.genes will be used
##             ...       -> additional parameters for heatmap
##
## Gustavo H. Esteves
## 27/05/07
##
##


hierMde <- function(data, group=c("C", "R", "B")[1],
distance="correlation", method="complete", doHeat=TRUE, sLabelID="SAMPLE",
gLabelID="GeneName", idxTest=1, adjP="BH", nDEgenes=0.05, ...) {
    
    
    ## Some some basic tests
    if(!(group %in% c("C", "R", "B")))
        stop("Parameter 'group' must be 'C', 'R' or 'B'.")
    
    if(!(method %in% c("ward", "single", "complete", "average", "mcquitty",
    "median", "centroid")))
        stop("Parameter method must be 'ward', 'single', 'complete', 'average',
        'mcquitty', 'median' or 'centroid'.")
    
    if(!(distance %in% c("euclidean", "maximum", "manhattan", "canberra",
    "binary", "pearson", "correlation", "spearman")))
        stop("Parameter 'distance' must be 'euclidean', 'maximum', 'manhattan',
        'canberra', 'binary', 'pearson', 'correlation', 'spearman'.")
    
    if(!is.element(adjP, c("Bonferroni", "Holm", "Hochberg", "SidakSS",
    "SidakSD", "BH", "BY", "none")))
        stop("Incorrect value for parameter adjP see help pages for or
        'mt.rawp2adjp'")
    
    
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
    zLim <- range(tmpMatrix)
    zLim <- c(-max(abs(zLim)), max(abs(zLim)))
    
    
    if(group == "R") {
        if(doHeat)
            heatmap(tmpMatrix, distfun=function(c) amap::Dist(c,
            method=distance), hclustfun=function(d) hclust(d, method=method),
            Colv=NA, zlim=zLim, col=maigesPack:::greenRed(), scale="none", ...)
        
        else
            plot(hclust(amap::Dist(tmpMatrix, method=distance), method=method))
        
    }
    else if(group == "C") {
        if(doHeat)
            heatmap(tmpMatrix, distfun=function(c) amap::Dist(c,
            method=distance), hclustfun=function(d) hclust(d, method=method),
            Rowv=NA, zlim=zLim, col=maigesPack:::greenRed(), scale="none", ...)
        
        else
            plot(hclust(dist(t(tmpMatrix), method=distance), method=method))
        
    }
    else {
        if(!doHeat)
            stop("With group 'B' you can't plot the hierarquical branch.")
        else
            heatmap(tmpMatrix, distfun=function(c) amap::Dist(c,
            method=distance), hclustfun=function(d) hclust(d, method=method),
            zlim=zLim, col=maigesPack:::greenRed(), scale="none", ...)
        
    }
}
