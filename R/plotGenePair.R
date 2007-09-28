## Generic function to do scatter plots between a pair of genes from relevance
## networks analysis
##
## Parameters: obj   -> object of maigesRelNet class
##             gene1 -> char string giving the first gene name
##             gene2 -> char string giving the first gene name
##             posL  -> x and y position for the legend
##             rCor  -> Is correlation robust? Defaults to TRUE
##
## Gustavo H. Esteves
## 15/05/07
##
## Version: 1.0
##


plotGenePair <- function(obj, gene1, gene2, posL=NULL, rCor=TRUE) {
    
    
    ## Getting table of Ws
    wTable <- obj@W
    
    
    ## Geting indexes of genes and sample types
    idxG1 <- rownames(wTable) == gene1
    idxG2 <- rownames(wTable) == gene2
    idxS1 <- colnames(wTable) %in% obj@types[1]
    idxS2 <- colnames(wTable) %in% obj@types[2]
    
    
    ## Calculating correlation values
    if(rCor)
        corr <- list(robustCorr(wTable[idxG1, idxS1], wTable[idxG2, idxS1]),
        robustCorr(wTable[idxG1, idxS2], wTable[idxG2, idxS2]))
    
    else
        corr <- list(list(Rcorr=obj@Corr1[idxG1, idxG2]),
        list(Rcorr=obj@Corr2[idxG1, idxG2]))
    
    ## Geting vectors of interest and calculating the regression line for sample
    ## type 1
    x <- wTable[idxG1, idxS1]
    y <- wTable[idxG2, idxS1]
    if(rCor)
        reg1 <- lm(y[-corr[[1]]$idx] ~ x[-corr[[1]]$idx])$coefficients
    else
        reg1 <- lm(y ~ x)$coefficients
    
    
    ## Geting vectors of interest and calculating the regression line for sample
    ## type 2
    x <- wTable[idxG1, idxS2]
    y <- wTable[idxG2, idxS2]
    if(rCor)
        reg2 <- lm(y[-corr[[2]]$idx] ~ x[-corr[[2]]$idx])$coefficients
    else
        reg2 <- lm(y ~ x)$coefficients
    
    
    ## Getting limits for the plot  
    limitsX <- range(wTable[idxG1, ])
    limitsY <- range(wTable[idxG2, ])
    limitsY <- c(limitsY[1], limitsY[2]+((limitsY[2]-limitsY[1])*0.2))
    
    
    ## Geting positins to print corr values
    if(is.null(posL))
        posL <- c(limitsX[1], limitsY[2])
    
    
    ## Pasting the correlation values in the labels
    labels <- obj@types
    labels[1] <- paste(labels[1], ", R =", round(corr[[1]]$Rcor, 3), sep="")
    labels[2] <- paste(labels[2], ", R =", round(corr[[2]]$Rcor, 3), sep="")
    
    
    
    
    ## Ploting the scatter plot
    plot(wTable[idxG1, idxS1], wTable[idxG2, idxS1], pch=19, xlim=limitsX,
    ylim=limitsY, xlab=gene1, ylab=gene2)
    
    points(wTable[idxG1, idxS2], wTable[idxG2, idxS2], pch=22)
    abline(reg1, lty=1)
    abline(reg2, lty=2)
    legend(posL[1], posL[2], labels, pch=c(19, 22))
    
}
