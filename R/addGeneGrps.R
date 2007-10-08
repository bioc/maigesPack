## Function to load gene groups into maigesPreRaw object
##
## Parameters: data    -> object of maigesRaw class
##             folder  -> char string specyfing the directory of genes
##                        groups
##             ext     -> extension of the files, defaults to '.txt'
##
## Gustavo Esteves
## 14/05/07
##
##


addGeneGrps <- function(data, folder="./", ext=".txt") {
    
    
    ## adjusting the 'folder' and 'ext' parameters
    tmp1 <- strsplit(folder, NULL)
    if (tmp1[[1]][length(tmp1[[1]])] != "/")
        folder <- paste(folder, "/", sep="")
    if (!is.null(ext)) {
        tmp2 <- strsplit(ext, NULL)
        if (tmp2[[1]][1] != ".")
            ext <- paste(".", ext, sep="")
    }
    
    
    ## Doing a test onto the object
    if(class(data) != "maigesPreRaw")
        stop("This function can't be applied to this class of objects!")
    
    
    ## Reading and testing the presence of the Gene groups
    groups <- dir(folder, pattern=paste(ext, "$", sep=""))
    nGroups <- length(groups) ## Number of groups
    for (i in groups) {
        
        name <- i
        
        if (!is.null(ext))
            i <- as.character(strsplit(i, ext, fixed=TRUE)[[1]][1])
        else
            i <- as.character(i)
        
        ## Changig punctuation and blanck characters from the name
        i <- sub("^[_]", "", i)
        i <- sub("[_]$", "", i)
        i <- gsub("[[:blank:]]+", "_", i)
        i <- gsub("[[:punct:]]+", "_", i)
        
        ## reading the files
        tmpGenes <- scan(paste(folder, name, sep=""), what="character",
        sep="\n", quiet=TRUE)
        
        ##tmpGenes <- as.vector(as.matrix(read.table(paste(folder, name,
        ##sep=""), sep="\t", skip=skip)))
        
        ## Removing blanks and repeats
        tmpGenes <- unique(tmpGenes)
        blk <- (tmpGenes == "")
        if(sum(blk) > 0)
            tmpGenes <- tmpGenes[-which(blk)]
        
        ## putting the genes read in the groupsGenes list
        if(!(i %in% names(data@GeneGrps)))
            data@GeneGrps[[i]] <- tmpGenes
        else
            warning(paste("The object already has a group named ", i,
            ". It was not added.", sep=""))
    }
    
    
    ## Atualizing the date slot and returning the object
    data@Date <- date()
    return(data)
    
}
