## Function to search by cliques (groups of genes) of classifiers,
## using exaustive search.
##
## Parameters: obj        -> Object of class maiges to do the calculations
##             sLabelID   -> Identification of the sample label to be used.
##             facToClass -> named list with 2 character vectors
##                           specifying the samples to be compared. If
##                           NULL the first 2 types of SlabelID are used
##             gNameID    -> Identification of gene name label ID.
##             geneGrp    -> char (or index) specifying the gene group to be
##                           tested (colnames of GeneGrps slot). If both
##                           geneGrp and path are NULL all genes are used.
##                           Defalts to 1 (first group).
##             path       -> char (or index) specifying the path to be tested
##                           (names of Paths slot). If both
##                           geneGrp and path are NULL all genes are used.
##                           Defalts to NULL.
##             nGenes     -> Number of genes in the classifier
##             kn         -> number of neighbours for the knn method
##
## Gustavo Esteves
## Adapted from Elier Cristo's functions
## 27/05/07
##
##


classifyKNN <- function(obj=NULL, sLabelID="Classification", facToClass=NULL,
gNameID="GeneName", geneGrp=1, path=NULL, nGenes=3, kn=5) {
    
    
    ## Testing some things
    if(nGenes < 2)
        stop("nGenes must be greater than 2.")
    if(!is.null(geneGrp) & !is.null(path))
        stop("You must specify only one of geneGrp and path, or leave
        both NULL.")
    
    
    ## Defining 1 additional function to do more iterations
    oneMoreIter <- function(tab, idxAnt, kn) {
        
        samp <- as.factor(colnames(tab))
        
        ng <- length(tab[, 1])
        nClique <- ncol(idxAnt)+1
        resCV <- NULL
        resSVD <- numeric()
        indexes <- NULL
        for (i in 1:nrow(idxAnt)) {
            g <- idxAnt[i, ]
            for (j in (1:ng)[-g]) {
                
                ## Doing first calculation
                if(i == 1 & j == (1:ng)[-g][1]) {
                    resCV <- c(resCV,
                    sum(class::knn.cv(train=as.data.frame(t(tab))[, c(g, j)],
                    cl=samp, k=kn) == samp))
                    
                    indexes <- rbind(indexes, c(g, j))
                }
                ## Verifying if the clique was tested already
                test <- 0
                for(k in 1:nrow(indexes)) {
                    if(sum(is.element(c(g, j), indexes[k, ])) == nClique)
                        test <- test+1
                }
                
                if(test == 0) {
                    resCV <- c(resCV,
                    sum(class::knn.cv(train=as.data.frame(t(tab))[, c(g, j)],
                    cl=samp, k=kn) == samp))
                    
                    indexes <- rbind(indexes, c(g, j))
                }
            }
        }
        
        idxGood <- sort(resCV, decreasing=TRUE, index.return=TRUE)$ix

        return(list(CV=resCV[idxGood], SVD=resSVD[idxGood],
        cliques=indexes[idxGood, ]))
        
    }
    
    
    
    ## Getting all labels for genes and samples
    allGenes <- getLabels(obj, gNameID, FALSE)
    allGenes[obj@BadSpots] <- paste(allGenes[obj@BadSpots], "(*)")
    samples <- getLabels(obj, sLabelID)
    
    if(is.null(facToClass)) {
        facToClass <- as.list(unique(samples)[1:2])
        names(facToClass) <- unique(samples)[1:2]
    }
    
    
    
    ## Getting samples from facToClass
    idxGrp1 <- samples %in% facToClass[[1]]
    idxGrp2 <- samples %in% facToClass[[2]]
    
    table <- calcW(obj)
    colnames(table)[idxGrp1] <- rep(names(facToClass)[1], sum(idxGrp1))
    colnames(table)[idxGrp2] <- rep(names(facToClass)[2], sum(idxGrp2))
    rownames(table) <- allGenes
    
    
    
    ## Getting genes from gene group if specified (or path)
    if(!is.null(geneGrp)) {
        if(!is.numeric(geneGrp))
            geneGrp <- which(colnames(obj@GeneGrps) == geneGrp)
        if(sum(obj@GeneGrps[,geneGrp]) <= nGenes)
            stop(paste("   There are less than", nGenes,
                       "genes in the group!!"))
        idxTmp <- obj@GeneGrps[, geneGrp]
        table <- table[idxTmp, ]
    }
    else if(!is.null(path)) {
        if(length(nodes(obj@Paths[[path]])) <= nGenes)
            stop(paste("   There are less than", nGenes,
                       "genes in the path!!"))
        idxTmp <- rownames(table) %in% nodes(obj@Paths[[path]])
        table <- table[idxTmp, ]
    }
    
    
    ## Removing samples that wer not used
    idx <- !is.na(colnames(table))
    table <- table[, idx]
    
    
    ## Doing the exaustive search for classifiers
    tmp <- matrix(1:nrow(table), nrow(table), 1)
    for(i in 1:(nGenes-1)) {
        classCliques <- oneMoreIter(table, tmp, kn)
        tmp <- classCliques$cliques
    }
    
    geneCliques <- NULL
    for(i in 1:dim(classCliques$cliques)[1])
        geneCliques <- rbind(geneCliques,
        rownames(table)[classCliques$cliques[i, ]])
    
    
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
    
    
    ## Defining an object to return
    result <- new("maigesClass", W=table, CV=classCliques$CV,
    SVD=classCliques$SVD, cliques=geneCliques, cliques.idx=classCliques$cliques,
    Date=date(), method="k-neighbours classifier, by exaustive search",
    V.info=vInfo)
    
    return(result)
    
}
