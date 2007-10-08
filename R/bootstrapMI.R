## Function to calculate robust correlation coeficients (minimizing the
## influences of one outlier). This strategy uses an idea similar to
## leave-one-out of cross-validation, where we make a decision based on the
## "rule of three numbers".
##
## Parameters: x    -> matrix to calculate the bootstraped p-values of
##                     correlations (or maximum boot cor) between all
##                     pairs of rows from x. Return a square matrix with
##                     all p-values correlations (or maximum).
##                     If x is a matrix, the function calculates a matrix
##                     of p-values or maximuns, else
##                     x must be a vector and y must be specified
##             y    -> optional numerical vector. Must be specified if x
##                     is a vector. If x is a matrix, y is ignored
##             bRep -> number of bootstraps samples (by permutation) to
##                     generate
##             ret  -> of value to return? Must be 'p-value' ou 'max'
##
## Gustavo H. Esteves
## 14/05/07
##
##


bootstrapMI <- function(x, y=NULL, bRep, ret="p-value") {
    
    
    ## Doing basic tests
    if(!(ret %in% c("p-value", "max")))
        stop("Parameter 'ret' must be 'p-value' or 'max'.")  
    
    
    if(length(dim(x)) == 2) {
        
        ## Getting the dimensions of the matrix
        ns <- dim(x)
        
        ## Doing the MI calculations
        tmp <- array(NA, c(dim(x)[1], dim(x)[1], bRep))
        orig <- MI(x)
        for(j in 1:bRep) {
            if(ret == "max")
                tmp[, , j] <- MI(t(apply(x, 1, sample)))
            else if(ret == "p-value")
                tmp[, , j] <- MI(t(apply(x, 1, sample))) >= orig
        }
        
        if(ret == "max")
            ## Calculating maximum values
            res <- apply(tmp, c(1, 2), max)
        else if(ret == "p-value")
            ## Calculating p-value
            res <- apply(tmp, c(1, 2), sum)/bRep
        
        return(res)
    }
    else if(is.vector(x) & is.vector(y)) {
        
        tmp <- NULL
        orig <- MI(x, y)
        for(j in 1:bRep)
            tmp <- c(tmp, MI(sample(x), sample(y)))
        
        ## Calculating p-value
        if(ret == "p-value")
            res <- sum(tmp >= orig)/bRep
        else if(ret == "max")
            res <- max(tmp)
        
        return(res)
        
    }
    else
        stop("x must be a numerical matrix or vector, case where you must
        specify another vector y of same length.")
    
}
