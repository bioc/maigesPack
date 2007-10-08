## Function to calculate mutual information from every row of a matrix, or
## between 2 vectors
##
## Parameters: x -> matrix to calculate the MI between all pairs of rows from
##                  x. Return a square matrix with all MI's. If x is a matrix,
##                  the function calculates a matrix of robust correlations,
##                  else x must be a vector and y must be specified as another
##                  vector of same length and the correlation between they both
##                  are calculated
##             y -> another vector, optional. Must be especified it x is a
##                  vector.
##             k -> number of the neighbour to be used, in the calculation of
##                  the MI value
##
## Gustavo H. Esteves
## 15/05/07
##
##


MI <- function(x, y=NULL, k=1) {
    
    
    ## Testing if a matrix was specified
    if(length(dim(x)) == 2) {
        ## Getting the dimensions of the matrix
        ns <- dim(x)
        
        ## Scaling the rows of the matrix
        x <- t(scale(t(x)))
        
        ## Doing the cor calculations from C compiled code
        tmp <- .C("Minfo", as.double(as.vector(t(x))), as.integer(ns[1]),
        as.integer(ns[2]), as.integer(k), as.double(rep(0, (ns[1])^2)),
        DUP=FALSE, PACKAGE="maigesPack")
        
        ## Defining the result matrix
        tmp <- matrix(tmp[[5]], ns[1], ns[1], TRUE)
        if(!is.null(rownames(x))) {
            rownames(tmp) <- rownames(x)
            colnames(tmp) <- rownames(x)
        }
        
        return(tmp)
    }
    else if(is.vector(x) & !is.null(y)) {
        
        ## Doing simple tests
        if(!is.vector(y))
            stop("y must be a vector of same length as x!")
        if(length(x) != length(y))
            stop("x and y must have same length!")
        
        ## calculating the correlation value
        calc <- rbind(x, y)
        
        ## Scaling the rows of the matrix
        calc <- t(scale(t(calc)))
        
        ns <- dim(calc)
        tmp <- .C("Minfo", as.double(as.vector(t(calc))), as.integer(ns[1]),
        as.integer(ns[2]), as.integer(k), as.double(rep(0, (ns[1])^2)),
        DUP=FALSE, PACKAGE="maigesPack")
        
        ## Defining the result
        tmp <- matrix(tmp[[5]], ns[1], ns[1], TRUE)
        
        return(tmp[1, 2])
        
    }
    else
        stop("x must be a numerical matrix or vector, case where you must
        specify another vector y of same length.")
    
}
