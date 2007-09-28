## Function to load pathway files in the TGF format
##
## Parameters: data   -> object of maigesRaw class
##             folder -> directory to load the TGF files
##             ext    -> extension of the files, defaults to '.tgf'
##
## Gustavo Esteves
## 14/05/07
##
## Version: 1.0
##


addPaths <- function(data, folder="./", ext=".tgf") {
    
    
    ## adjusting the 'folder' and 'ext' parameters
    tmp1 <- strsplit(folder, NULL)
    if (tmp1[[1]][length(tmp1[[1]])] != "/")
        folder <- paste(folder, "/", sep="")
    if (!is.null(ext)) {
        tmp2 <- strsplit(ext, NULL)
        if (tmp2[[1]][1] != ".")
            ext <- paste(".", ext, sep="")
    }
    
    
    ## Doing a test onto the class of object
    if(class(data) != "maigesPreRaw")
        stop("This function can't be applied to this class of objects!")
    
    
    ## Read files from dir
    paths <- dir(folder, pattern=paste(ext, "$", sep=""))
    for (i in paths) {
        
        if (!is.null(ext))
            name <- as.character(strsplit(i, ext, fixed=TRUE)[[1]][1])
        else
            name <- as.character(i)
        
        name <- sub("^[_]", "", name)
        name <- sub("[_]$", "", name)
        name <- gsub("[[:blank:]]+", "_", name)
        name <- gsub("[[:punct:]]+", "_", name)
        
        ## Reading TGF files
        fileTgf <- scan(paste(folder, i, sep=""), what="character", sep="\n",
        quiet=TRUE)
        
        div <- which(fileTgf == "#")
        
        ## Getting gene names
        tmpGenes <- NULL
        for(j in 1:(div-1)) {
            tmp <- sub("^[[:digit:]]+[[:blank:]]+", "", fileTgf[j])
            tmpGenes <- c(tmpGenes, tmp)
        }
        
        vertices <- unique(tmpGenes)
        arestas <- vector("list", length=length(vertices))
        names(arestas) <- vertices
        
        
        ## Catching interactions and removing repeated ones
        interactions <- unique(fileTgf[(div+1):length(fileTgf)])
        
        
        ## Getting interactions between genes
        nArestas <- 0
        for (k in 1:length(tmpGenes)) {
            arestasTmp <- NULL
            weightsTmp <- NULL
            for(j in 1:length(interactions)) {
                aux <- strsplit(interactions[j], "\\s+", perl=TRUE)
                if(aux[[1]][1] == k & aux[[1]][2] %in% (1:length(tmpGenes))) {
                    arestasTmp <- c(arestasTmp,
                    tmpGenes[as.numeric(aux[[1]][2])])
                    
                    weightsTmp <- c(weightsTmp, as.numeric(aux[[1]][3]))
                    nArestas <- nArestas+1
                }
            }
            if(is.null(arestasTmp))
                arestasTmp <- character()
            if(is.null(weightsTmp))
                weightsTmp <- numeric()
            
            arestas[[tmpGenes[k]]] <- list(edges=arestasTmp, weights=weightsTmp)
        }
        
        ## Defining the object with the graphs
        if(!(name %in% names(data@Paths)))
            data@Paths[[name]] <- new("graphNEL", vertices, edgeL=arestas,
            edgemode="directed")
        else
            warning(paste("The object already has a Path named ", name,
            ".It was not added.", sep=""))
        
    }
    
    
    ## Actualising the date and returning the object
    data@Date <- date()
    
    return(data)
    
}
