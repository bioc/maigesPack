## Function to calculate p-values by bootstrap using a function implemented
## in C row on microarray's data... This function calculate boot p.values
## for t statistic
##
## Parameters: x    -> numerical matrix to be t bootstraped
##             k    -> number of resamplings
##             obs1 -> indexes of first group
##             obs2 -> indexes of second group
##             ...  -> additional parameters for t test
##
## Gustavo H. Esteves
## 14/05/07
##
##


bootstrapT <- function(x, k=20000, obs1, obs2, ...) {
    
    
    ## Defining an additional function
    Test <- function(x, obs1, obs2, ...) {
        tmpOut <- t.test(x[obs1], x[obs2], ...)
        return(unname(tmpOut$statistic))
    }
    
    
    ## Writing indexes as integer if they are logical
    if(is.logical(obs1))
        obs1 <- which(obs1)
    if(is.logical(obs2))
        obs2 <- which(obs2)
    
    nRows <- length(x[, 1])
    
    if(nRows == 1) {
        tableObs1 <- as.matrix(t(x[, obs1]))
        tableObs2 <- as.matrix(t(x[, obs2]))
    }
    else {
        tableObs1 <- x[, obs1]
        tableObs2 <- x[, obs2]
    }
    
    
    ## geting additional parameters (for function t.test)
    vEq <- FALSE
    add <- list(...)
    if("var.equal" %in% names(add))
        vEq <- add$var.equal
    
    
    ## calculating the original means of the data.
    meanObs1 <- apply(tableObs1, 1, mean, na.rm=TRUE);
    meanObs2 <- apply(tableObs2, 1, mean, na.rm=TRUE);
    
    
    ## Calculating the original T statistic.
    tObs <- apply(x, 1, Test, obs1, obs2, ...)
    
    
    ## Doing the bootstraps by compiled C code.
    pBootT <- .C("bootT", as.double(as.vector(t(x))), as.integer(k),
    as.integer(obs1), as.integer(obs2), as.integer(nRows),
    as.integer(length(obs1)), as.integer(length(obs2)),
    as.integer(vEq), as.double(rep(1, nRows)), #DUP=FALSE,
    PACKAGE="maigesPack")[[9]]
    
    
    ## Returning the object
    return(cbind(meanDif=meanObs1-meanObs2, Statistic=tObs,  P.value=pBootT))
    
}
