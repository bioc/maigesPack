## Generic function to do hierarquical cluster analysis
##
## Parameters: data      -> maiges (or maigesRaw) class object
##             group     -> Type of grouping: by rows ('R'), columns ('C') or
##                          both ('B')
##             distance  -> char string giving the type of distance to use from
##                          function Dist (lib amap)
##             method    -> char string specifying the linkage method for the
##                          hierarquical cluster
##             doHeat    -> logical indicanting to do or not the heatmap
##             sLabelID  -> Sample label id
##             gLabelID  -> Gene label id
##             rmSamples -> char list specifying samples to be removed
##             rmGenes   -> char list specifying genes to be removed
##             rmBad     -> remove bad spots?
##             geneGrp   -> char with name (or index) of the gene group to be
##                          ploted
##             path      -> char with name (or index) of the pathway to be
##                          ploted
##             ...       -> additional parameters for heatmap
##
## Gustavo H. Esteves
## 15/05/07
##
## Version: 1.0
##


hierM <- function(data, group=c("C", "R", "B")[1],
distance="correlation", method="complete", doHeat=TRUE, sLabelID="SAMPLE",
gLabelID="GeneName", rmGenes=NULL, rmSamples=NULL, rmBad=TRUE, geneGrp=NULL,
path=NULL, ...) {
    
    
    ## Some some basic tests
    if(!(group %in% c("C", "R", "B")))
        stop("Parameter 'group' must be 'C', 'R' or 'B'.")
    
    if(!is.null(geneGrp) & !is.null(path))
        stop("You must specify only geneGrp OR path (or none of them)")
    
    if(!(method %in% c("ward", "single", "complete", "average", "mcquitty",
    "median", "centroid")))
        stop("Parameter method must be 'ward', 'single', 'complete',
        'average', 'mcquitty', 'median' or 'centroid'.")
    
    if(!(distance %in% c("euclidean", "maximum", "manhattan", "canberra",
    "binary", "pearson", "correlation", "spearman")))
        stop("Parameter 'distance' must be 'euclidean', 'maximum', 'manhattan',
        'canberra', 'binary', 'pearson', 'correlation', 'spearman'.")
    
    
    ## Removing bad spots
    if(rmBad)
        data <- data[!data@BadSpots, ]
    
    
    ## Filtering genes samples and genes...
    if(!is.null(rmSamples)) {
        tmp <- getLabels(data, sLabelID)
        filter <- !is.element(tmp, rmSamples)
        data <- data[, filter]
    }
    if(!is.null(rmGenes)) {
        tmp <- getLabels(data, gLabelID, FALSE)
        filter <- !is.element(tmp, rmGenes)
        data <- data[filter, ]
    }
    
    
    ## Selecting the genes according to the gene groups
    if(!is.null(geneGrp)) {
        if(!is.numeric(geneGrp))
            geneGrp <- which(colnames(data@GeneGrps) == geneGrp)
        
        idxTmp <- data@GeneGrps[, geneGrp]
        
        if(sum(idxTmp) <= 1)
            stop("Group with less than one gene, it's impossible constructing
            a cluster.")
        
        tmpMatrix <- calcW(data[idxTmp, ])
        colnames(tmpMatrix) <- getLabels(data, sLabelID)
        tmp <- getLabels(data, gLabelID, FALSE)[idxTmp]
        tmp[data@BadSpots[idxTmp]] <- paste(tmp[data@BadSpots[idxTmp]], "(*)")
        rownames(tmpMatrix) <- tmp
        
    }
    ## Selecting the genes according to path
    else if(!is.null(path)) {
        if(!is.numeric(path)) {
            if(gLabelID != data@Paths$Glabel)
                stop("gLabelID is different from data@Paths$Glabel.")
            path <- which(names(data@Paths) == path)
        }

        idxTmp <- getLabels(data, data@Paths$Glabel, FALSE) %in%
        nodes(data@Paths[[path]])
        
        if(sum(idxTmp) <= 1)
            stop("Path with less than one gene, it's impossible constructing
            a cluster.")
        
        tmpMatrix <- calcW(data[idxTmp, ])
        colnames(tmpMatrix) <- getLabels(data, sLabelID)
        tmp <- getLabels(data, gLabelID, FALSE)[idxTmp]
        tmp[data@BadSpots[idxTmp]] <- paste(tmp[data@BadSpots[idxTmp]], "(*)")
        rownames(tmpMatrix) <- tmp
    }
    else {
        tmpMatrix <- calcW(data)
        colnames(tmpMatrix) <- getLabels(data, sLabelID)
        tmp <- getLabels(data, gLabelID, FALSE)
        tmp[data@BadSpots] <- paste(tmp[data@BadSpots], "(*)")
        rownames(tmpMatrix) <- tmp
    }
    
    
    ## Removing non finite genes
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
            plot(hclust(amap::Dist(t(tmpMatrix), method=distance),
            method=method))
        
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
