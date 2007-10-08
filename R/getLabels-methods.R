## Function extract labels given a label name (for genes and samples)
##
## Parameters: obj     -> object to look for labels. Methods defined
##                        for maigesRaw, maiges, RGList, MAList,
##                        marrayRaw and marrayNorm
##             labelID -> string with label name to be searched
##             sLabel  -> logical indicating search in the sample
##                        labels, defaults to TRUE. If FALSE search
##                        for gene labels
##
## Gustavo H. Esteves
## 15/05/07
##
##


## Defining a default method
getLabels.default <- function(obj=NULL, labelID=NULL, sLabel=TRUE) {
    
    if(sLabel)
        labels <- eval(parse(text=paste("as.character(obj@Slabels$", labelID,
        ")", sep="")))
    
    else
        labels <- eval(parse(text=paste("as.character(obj@Glabels$", labelID,
        ")", sep="")))
    
    return(labels)
    
}


## Defining a maigesDE method
getLabels.maigesDE <- function(obj=NULL, labelID=NULL, sLabel=TRUE) {
    
    if(sLabel)
        labels <- eval(parse(text=paste("as.character(obj@SampleInfo$", labelID,
        ")", sep="")))
    
    else
        labels <- eval(parse(text=paste("as.character(obj@GeneInfo$", labelID,
        ")", sep="")))
    
    return(labels)
    
}


## Defining a maigesDEcluster method
getLabels.maigesDEcluster <- getLabels.maigesDE


## Defining a RGList method
getLabels.RGList <- function(obj=NULL, labelID=NULL, sLabel=TRUE) {
    
    if(sLabel)
        labels <- eval(parse(text=paste("as.character(obj$targets$", labelID,
        ")", sep="")))
    
    else
        labels <- eval(parse(text=paste("as.character(obj$genes$", labelID, ")",
        sep="")))
    
    return(labels)
    
}


## Defining a MAList method
getLabels.MAList <- function(obj=NULL, labelID=NULL, sLabel=TRUE) {
    
    if(sLabel)
        labels <- eval(parse(text=paste("as.character(obj$targets$", labelID,
        ")", sep="")))
    
    else
        labels <- eval(parse(text=paste("as.character(obj$genes$", labelID, ")",
        sep="")))
    
    return(labels)
    
}


## Defining a marrayRaw method
getLabels.marrayRaw <- function(obj=NULL, labelID=NULL, sLabel=TRUE) {
    
    if(sLabel)
        labels <- eval(parse(text=paste("as.character(obj@maTargets@maInfo$",
        labelID, ")", sep="")))
    
    else
        labels <- eval(parse(text=paste("as.character(obj@maGnames@maInfo$",
        labelID, ")", sep="")))
    
    return(labels)
    
}


## Defining a marrayNorm method
getLabels.marrayNorm <- function(obj=NULL, labelID=NULL, sLabel=TRUE) {
    
    if(sLabel)
        labels <- eval(parse(text=paste("as.character(obj@maTargets@maInfo$",
        labelID, ")", sep="")))
    
    else
        labels <- eval(parse(text=paste("as.character(obj@maGnames@maInfo$",
        labelID, ")", sep=""))) 
    
    return(labels)
    
}
