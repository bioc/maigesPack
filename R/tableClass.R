## Function to write tables of classifiers genes (cliques)
##
## Parameters: classComp -> object of type maigesClass
##             file      -> string giving the file name to be saved (without
##                          extension)
##             type      -> type of file to be saved. 'HTML' (default) or 'CSV'
##             nCliques  -> number of cliques to be saved
##
## Gustavo Esteves
## 21/11/07
##
##


tableClass <- function(classComp=NULL, file="./class_result", type=c("HTML",
"CSV")[1], nCliques=NULL) {

    
    ## Doing a basic test
    if(!(type %in% c("HTML","CSV")))
        stop("type parameter must be 'HTML' or 'CSV'.")
    
    
    ## Defining all cliques if show is null
    if(is.null(nCliques))
        nCliques <- dim(classComp@cliques)[1]
    
    
    ## Verifying if the method used was lda (because of svd values)
    if(length(grep("lda", classComp@method, ignore.case=TRUE)) > 0)
        tmpDf <- data.frame(CV=classComp@CV, SVD=classComp@SVD,
        Clique=matrix(apply(classComp@cliques, 1, paste, collapse=", "),
        dim(classComp@cliques)[1], 1))
    
    else
        tmpDf <- data.frame(CV=classComp@CV,
        Clique=matrix(apply(classComp@cliques, 1, paste, collapse=", "),
        dim(classComp@cliques)[1],1))
    
    
    ## Saving the table
    if(type == "HTML") {
        
        R2HTML::HTML(tmpDf[1:nCliques, ], file=paste(file, ".html", sep=""),
                     innerBorder=1, append=FALSE)
        
    }
    else
        write.table(tmpDf[1:nCliques, ], file=paste(file, ".tsv", sep=""),
        sep="\t", quote=FALSE, row.names=FALSE)
    
}
