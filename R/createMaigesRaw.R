## Function to convert maigesPreRaw objects into maigesRaw
##
## Parameter: PreRaw           -> object of class maigesPreRaw
##            green.dataField  -> string with the name of the object from Data
##                                slot to be used as foreground information from
##                                channel 1
##            greenB.dataField -> string with the name of the object from Data
##                                slot to be used as background information from
##                                channel 1
##            red.dataField    -> string with the name of the object from Data
##                                slot to be used as foreground information from
##                                channel 2
##            redB.dataField   -> string with the name of the object from Data
##                                slot to be used as background information from
##                                channel 2
##            flag.dataField   -> string with the name of the object from Data
##                                slot to be used as flags information
##            gLabelGrp        -> string with the gene label to match gene
##                                groups
##            gLabelPath       -> string with the gene label to match gene
##                                networks
##
##
##
## Gustavo H. Esteves
## 14/05/07
##
## Version: 1.0
##


createMaigesRaw <- function(PreRaw, greenDataField, greenBackDataField,
redDataField, redBackDataField, flagDataField, gLabelGrp, gLabelPath) {
    
    
    ## Defining the dimensions of the object
    nFiles <- dim(PreRaw)[2]
    nLines <- dim(PreRaw)[1]
    
    
    
    ##
    ## Defining the data matrices
    ## 
    
    ## Defining indexes for references and samples in both channels
    tmp <- c("green", "red") %in% tolower(PreRaw@Slabels$Ref)
    if(sum(tmp) == 2) {
        ch1.ref <- (tolower(PreRaw@Slabels$Ref) == "green")
        
        ## Defining data matrices and vectors to receive intensity data values
        ## and flags
        Sf <- matrix(nrow=nLines, ncol=nFiles)
        Sb <- matrix(nrow=nLines, ncol=nFiles)
        Rf <- matrix(nrow=nLines, ncol=nFiles)
        Rb <- matrix(nrow=nLines, ncol=nFiles)
        
        ## Loading sample intensity values
        Sf[, !ch1.ref] <- PreRaw@Data[[greenDataField]][, !ch1.ref]
        Sf[, ch1.ref] <- PreRaw@Data[[redDataField]][, ch1.ref]
        ## loading sample background values
        Sb[, !ch1.ref] <- PreRaw@Data[[greenBackDataField]][, !ch1.ref]
        Sb[, ch1.ref] <- PreRaw@Data[[redBackDataField]][, ch1.ref]
        
        ## loading reference intensity values
        Rf[, ch1.ref] <- PreRaw@Data[[greenDataField]][, ch1.ref]
        Rf[, !ch1.ref] <- PreRaw@Data[[redDataField]][, !ch1.ref]
        ## loading reference background values
        Rb[, ch1.ref] <- PreRaw@Data[[greenBackDataField]][, ch1.ref]
        Rb[, !ch1.ref] <- PreRaw@Data[[redBackDataField]][, !ch1.ref]
        
        ## Defining the vector of dyes
        Rdye <- tolower(PreRaw@Slabels$Ref)
        Sdye <- vector(mode="character", length=dim(PreRaw)[2])
        Sdye[ch1.ref] <- Rdye[!ch1.ref]
        Sdye[!ch1.ref] <- Rdye[ch1.ref]
    }
    else if(c("green", "red")[tmp] == "green") {
        
        ## Loading sample intensity values
        Sf <- PreRaw@Data[[redDataField]]
        ## loading sample background values
        Sb <- PreRaw@Data[[redBackDataField]]
        
        ## loading reference intensity values
        Rf <- PreRaw@Data[[greenDataField]]
        ## loading reference background values
        Rb <- PreRaw@Data[[greenBackDataField]]
        
        ## Defining the vector of dyes
        Rdye <- tolower(PreRaw@Slabels$Ref)
        Sdye <- rep("red", length=dim(PreRaw)[2])
    }
    else if(c("green", "red")[tmp] == "red") {
        
        ## Loading sample intensity values
        Sf <- PreRaw@Data[[greenDataField]]
        ## loading sample background values
        Sb <- PreRaw@Data[[greenBackDataField]]
        
        ## loading reference intensity values
        Rf <- PreRaw@Data[[redDataField]]
        ## loading reference background values
        Rb <- PreRaw@Data[[redBackDataField]]
        
        ## Defining the vector of dyes
        Rdye <- tolower(PreRaw@Slabels$Ref)
        Sdye <- rep("green", length=dim(PreRaw)[2])
    }
    else
        stop("Something strange with your 'Ref' slot from sample label.")
    
    
    
    ##
    ## Matching gene groups from PreRaw to define the indexes matrix
    ## for the maigesRaw object
    ## 
    
    ## Getting gene labels from data object
    if(length(PreRaw@GeneGrps) > 0) {
        dataGenes <- PreRaw@Glabels[[gLabelGrp]]
        
        groupsGenes <- NULL ## Defining a NULL object to put the gene indexes
        for (i in 1:length(PreRaw@GeneGrps)) {
            idxTmp <- (dataGenes %in% PreRaw@GeneGrps[[i]])
            groupsGenes <- cbind(groupsGenes, idxTmp)
        }
        colnames(groupsGenes) <- names(PreRaw@GeneGrps)
    }
    else
        groupsGenes <- matrix(numeric(), nrow=0, ncol=0)
    
    
    
    ##
    ## Matching gene networks from PreRaw to define the graphs
    ## for the maigesRaw object
    ## 
    
    ## Getting gene labels from data object
    if(length(PreRaw@Paths) > 0) {

        dataGenes <- PreRaw@Glabels[[gLabelPath]]
        
        ## Defining a list to put the gene indexes
        pathGenes <- list(Glabel=gLabelPath)
        for (i in names(PreRaw@Paths)) {
            
            ## Getting gene names
            tmpGenes <- nodes(PreRaw@Paths[[i]])
            tmpGenes <- gsub("[[:blank:]]+", "", tmpGenes)
            tmpGenes <- toupper(gsub("[[:punct:]]+", "", tmpGenes))
            
            ## Indexing with dataGenes and removing all nodes not present into
            ## dataset
            idxTmp <- !(tmpGenes %in% toupper(dataGenes))
            if(sum(idxTmp) > 0) {
                tmpGenes <- tmpGenes[!idxTmp]
                newGraph <- removeNode(nodes(PreRaw@Paths[[i]])[idxTmp],
                PreRaw@Paths[[i]])
            }
            else
                newGraph <- PreRaw@Paths[[i]]
            
            
            ## Changing node labels from graphs
            nodes(newGraph) <- tmpGenes
            
            if(numEdges(newGraph) > 0)
                pathGenes[[i]] <- newGraph
            
        }
    }
    else
        pathGenes <- list()
    
    
    
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
    
    
    
    ##
    ## Creating and returning the maigesRaw object
    ## 
    
    RawObject <- new("maigesRaw", Sf=Sf, Sb=Sb, Sdye=Sdye, Rf=Rf, Rb=Rb,
    Rdye=Rdye, Flag=PreRaw@Data[[flagDataField]], UseSpots=matrix(TRUE,
    nrow(Sf), ncol(Sf)), GeneGrps=groupsGenes, Paths=pathGenes,
    BadSpots=PreRaw@BadSpots, Layout=PreRaw@Layout, Slabels=PreRaw@Slabels,
    Glabels=PreRaw@Glabels, Notes=PreRaw@Notes, Date=date(), V.info=vInfo)
    
    
    return(RawObject)
    
}
