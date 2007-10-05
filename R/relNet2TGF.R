## Methods for relNet2TGF generic function
##
##
## Gustavo Esteves
## 15/05/07
##
##


## Defining a maigesRelNetM method
relNet2TGF.maigesRelNetM <- function(data, dir="./", filenames=c("group1.tgf",
"group2.tgf", "difPvalue.tgf"), pValue=0.05, ...) {
    
    ## Getting gene names as edges
    arestas <- rownames(data@DifP)
    for (i in filenames) {
        
        write.table(arestas, paste(dir, i, sep=""), col.names=FALSE,
        quote=FALSE)
        
        write.table("#", paste(dir, i, sep=""), col.names=FALSE,
        row.names=FALSE, append=TRUE, quote=FALSE)
        
    }
    
    for (j in 1:(length(arestas)-1)) {
        tmp <- upper.tri(data@DifP)
        idx <- which(data@DifP[j, tmp[j, ]] <= pValue)
        if (length(idx) > 0)
            for (k in (idx+j)) {
                
                write.table(paste(j, k, round(data@Corr1[j, k], 4)), paste(dir,
                filenames[1], sep=""), col.names=FALSE, row.names=FALSE,
                append=TRUE, quote=FALSE)
                
                write.table(paste(j, k, round(data@Corr2[j, k], 4)), paste(dir,
                filenames[2], sep=""), col.names=FALSE, row.names=FALSE,
                append=TRUE, quote=FALSE)
                
                write.table(paste(j, k, round(data@DifP[j, k], 4)), paste(dir,
                filenames[3], sep=""), col.names=FALSE, row.names=FALSE,
                append=TRUE, quote=FALSE)
                
            }
    }
    
    
}


## Defining a maigesRelNetB method
relNet2TGF.maigesRelNetB <- function(data, dir="./", filename="group.tgf",
corC=NULL, pValue=0.05, ...) {
    
    ## cutCor ->  must be numeric (in [0,1]) or 'max' (butte's method)
    
    ## Getting gene names as edges
    arestas <- rownames(data@Corr)

    write.table(arestas, paste(dir, filename, sep=""), col.names=FALSE,
    quote=FALSE)
    
    write.table("#", paste(dir, filename, sep=""), col.names=FALSE,
    row.names=FALSE, append=TRUE, quote=FALSE)
    
    if (!is.null(corC)) {
        if(corC == "max")
            corC <- max(data@maxB[upper.tri(data@maxB)])
        for (j in 1:(length(arestas)-1)) {
            tmp <- upper.tri(data@Corr)
            idx <- which(abs(data@Corr[j, tmp[j, ]]) >= corC)
            if (length(idx) > 0)
                for (k in (idx+j))
                    
                    write.table(paste(j, k, round(data@Corr[j, k], 4)),
                    paste(dir, filename, sep=""), col.names=FALSE,
                    row.names=FALSE, append=TRUE, quote=FALSE)
            
        }
    }
    else
        for (j in 1:(length(arestas)-1)) {
            tmp <- upper.tri(data@Pval)
            idx <- which(data@Pval[j, tmp[j, ]] <= pValue)
            if (length(idx) > 0)
                for (k in (idx+j))
                    
                    write.table(paste(j, k, round(data@Corr[j,k], 4)),
                    paste(dir, filename, sep=""), col.names=FALSE,
                    row.names=FALSE, append=TRUE, quote=FALSE)
            
        }
    
}
