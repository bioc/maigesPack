## Function to calculate robust correlation coeficients (minimizing the
## influences of one outlier). This strategy uses an idea similar to
## leave-one-out of cross-validation, where we make a decision based on the
## "rule of three numbers".
##
## Parameters: x -> matrix to calculate the correlations between all
##                  pairs of rows from x. Return a square matrix with
##                  all correlations. If x is a matrix, the function
##                  calculates a matrix of robust correlations, else
##                  x must be a vector and y must be specified as
##                  another vector of same length and the correlation
##                  between both are calculate
##             y -> optional numeric vector
##
## Gustavo H. Esteves
## 15/05/07
##
## Version: 1.0
##


robustCorr <- function(x, y=NULL) {
    
    
    if(length(dim(x)) == 2) {
        
        ## Getting the dimensions of the matrix
        ns <- dim(x)
        
        ## Doing the cor calculations from C compiled code
        tmp <- .C("robustCorr", as.double(as.vector(t(x))), as.integer(ns[1]),
        as.integer(ns[2]), as.double(rep(1, (ns[1])^2)), as.integer(rep(0,
        (ns[1])^2)), DUP=FALSE, PACKAGE="maigesPack")
        
        
        ## Defining the result matrix
        tmp1 <- matrix(tmp[[4]], ns[1], ns[1], TRUE)
        tmp2 <- matrix(tmp[[5]], ns[1], ns[1], TRUE)
        if(!is.null(rownames(x))) {
            rownames(tmp1) <- rownames(x)
            colnames(tmp1) <- rownames(x)
            rownames(tmp2) <- rownames(x)
            colnames(tmp2) <- rownames(x)
        }
        
        return(list(Rcor=tmp1, idx=tmp2))
        
    }
    else if(is.vector(x) & !is.null(y)) {
        
        ## Doing simple tests
        if(!is.vector(y))
            stop("y must be a vector of same length as x!")
        if(length(x) != length(y))
            stop("x and y must have same length!")
        
        ## calculating the correlation value
        calc <- rbind(x, y)
        ns <- dim(calc)
        
        tmp <- .C("robustCorr", as.double(as.vector(t(calc))),
        as.integer(ns[1]), as.integer(ns[2]), as.double(rep(1, (ns[1])^2)),
        as.integer(rep(0, (ns[1])^2)), DUP=FALSE, PACKAGE="maigesPack")
        
        
        ## Defining the result
        tmp1 <- matrix(tmp[[4]], ns[1], ns[1], TRUE)
        tmp2 <- matrix(tmp[[5]], ns[1], ns[1], TRUE)
        
        return(list(Rcor=tmp1[1,2], idx=tmp2[1,2]))
        
    }
    else
        stop("x must be a numerical matrix or vector, case where you must
        specify another vector y of same length.")
    
}
