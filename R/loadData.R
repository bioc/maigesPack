## Function to load microarray data (general data).
##
## Parameter: file.conf -> character string giving the configuration file
##                         name. All other parameters must be specified in this
##                         file
##
## Gustavo H. Esteves
## 27/05/07
##
##


loadData <- function(fileConf=paste(R.home(),
"library/maiges/doc/gastro/load_gastro.conf", sep="/")) {
    

    ## Defining some variables
    sampleFile <- character()
    dataDir <- character()
    geneMap <- character()
    ext <- character()
    headers <- character()
    skip <- numeric()
    sep <- character()
    gridR <- numeric()
    gridC <- numeric()
    printTipR <- numeric()
    printTipC <- numeric()
    datasetId <- character()
    
    
    ## Doing a basic test
    if(!is.character(fileConf))
        stop("fileConf must be a character string indicanting the
        parameter file!!")

    ## Parsing the configuration file
    eval(parse(fileConf))
    
    ## Opening a file to generate an output file
    outfile <- file("load.out", "w")
    sink(file=outfile, type=c("message"))
    
    
    ## Doing basic tests
    if (is.null(sampleFile))
        stop("Specify the 'sampleFile' parameter in the configuration file.")
    
    if (is.null(geneMap))
        stop("Specify the 'geneMap' parameter in the configuration file.")
    
    
    
    ## adjusting the 'dataDir' and 'ext' parameters
    tmp <- strsplit(dataDir, NULL)
    if (tmp[[1]][length(tmp[[1]])] != "/")
        dataDir <- paste(dataDir, "/", sep="")
    if (!is.null(ext)) {
        tmp <- strsplit(ext, NULL)
        if (tmp[[1]][1] != ".")
            ext <- paste(".", ext, sep="")
    }
    
    
    ## Create a list of files to be read
    cat("\nLoading the sample file ", sampleFile, ".\n\n", sep="",
    file=outfile)
    
    sampleTypes <- read.table(sampleFile, sep="\t", header=TRUE,
    as.is=TRUE)
    
    
    ## Doing basic tests
    if(sum(c("File", "Ref") %in% names(sampleTypes)) != 2)
        stop("fileConf must must have the fields 'File' and 'Ref'.")
    if(sum(!(tolower(sampleTypes$Ref) %in% c("green", "red"))) > 0)
        stop("Ref field from 'fileConf' must must contain only 'green'
        or 'red'.")
    
    
    ## Picking the number of files and lines
    nFiles <- length(sampleTypes[,1])
    nLines <- gridR*gridC*printTipR*printTipC
    
    
    ## Defining data matrices and vectors to receive the data
    nData <- list()
    headers2 <- gsub("[[:blank:]]+", ".", headers) ## changing
    ##BLANKS by .
    
    headers2 <- gsub("[[:punct:]]+", ".", headers2) ## changing PUNCTUATION by .
    
    
    ## Loading numerical values from files
    for (i in 1:nFiles) {
        
        fileTmp <- paste(dataDir,
        as.character(sampleTypes$File[i]), ext, sep="")
        
        cat("Loading file ", as.character(sampleTypes$File[i]), ext, ".\n",
        sep="", file=outfile)
        
        tmpTable <- read.table(fileTmp, header=TRUE, skip=skip,
        sep=sep, as.is=TRUE, nrows=nLines, check.names=FALSE,
        comment.char="", quote="\"")
        
        for (k in 1:length(headers))
            nData[[headers2[k]]] <- cbind(nData[[headers2[k]]],
            tmpTable[, headers[k]]) 
        
    }
    
    
    ## Create a layout for the chip
    dataLayout <- list(gridR=gridR, gridC=gridC, spotR=printTipR,
    spotC=printTipC, Nspots=gridR*gridC*printTipR*printTipC)
    
    
    ## Create the object for Sample labels
    sampleInfo <- as.data.frame(sampleTypes)
    
    
    ## Create the object for Gene labels 
    cat("\nLoading the gene map file ", geneMap, ".\n", sep="",
    file=outfile)
    genemap <- read.table(geneMap, sep="\t", header=TRUE, as.is=TRUE, quote="")
    
    
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
    
    
    ## Creating an object of type  maigesPreRaw
    cat("\nCreating a PreRaw object...\n", file=outfile)
    
    PreRawObj <- new("maigesPreRaw", Data=nData, Layout=dataLayout,
    Slabels=sampleInfo, BadSpots=rep(FALSE, length(nData[[1]][,1])),
    Glabels=genemap, Notes=datasetId, Date=date(), V.info=vInfo)
    
    
    ## closing the output file
    cat("\n\nOk, Object created succefuslly!!\n", file=outfile)
    sink(NULL, type=c("message"))
    
    
    ## Returning an object of type  maigesPreRaw
    return(PreRawObj)
    
}
