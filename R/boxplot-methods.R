## Define methods for boxplot generic function
##
## Gustavo H. Esteves
## 05/06/07
##
## Version: 1.1
##


## For maigesRaw class
boxplot.maigesRaw <- function(x, ...) {
    
    ## Testting if yvar was specified
    add <- list(...)
    if(length(add) == 0) {
        tmp <- as(x, "marrayRaw")
        tmp <- as(tmp, "marrayNorm")
        ## indexing ref labelled with green
        idx <- tolower(getLabels(x, "Ref")) == "red"
        if(sum(idx) > 0)
            maM(tmp)[, idx] <- -maM(tmp)[, idx]
        ## Testing if the object is indexed
        if(dim(x)[1] < x@Layout$Nspots) {
            Mvalues <- maM(tmp)
            Avalues <- maA(tmp)
            maM(tmp) <- matrix(NA, x@Layout$Nspots, dim(x)[2])
            maA(tmp) <- matrix(NA, x@Layout$Nspots, dim(x)[2])
            indexes <- as.numeric(rownames(x@Glabels))
            for (i in 1:dim(x)[1]) {
                maM(tmp)[indexes[i],] <- Mvalues[i,]
                maA(tmp)[indexes[i],] <- Avalues[i,]
            }
        }
        maBoxplot(tmp, ylab="W")
    }
    else {
        tmp <- as(x, "marrayRaw")
        ## Testing if the object is indexed
        if(dim(x)[1] < x@Layout$Nspots) {
            Mvalues <- maM(tmp)
            Avalues <- maA(tmp)
            maM(tmp) <- matrix(NA, x@Layout$Nspots, dim(x)[2])
            maA(tmp) <- matrix(NA, x@Layout$Nspots, dim(x)[2])
            indexes <- as.numeric(rownames(x@Glabels))
            for (i in 1:dim(x)[1]) {
                maM(tmp)[indexes[i],] <- Mvalues[i,]
                maA(tmp)[indexes[i],] <- Avalues[i,]
            }
        }
        maBoxplot(tmp, ...)
    }
    
}


## For maiges class
boxplot.maiges <- function(x, name=NULL, gLabelID=NULL, sLabelID=NULL,
gSamples=NULL, ...) {

    if(is.null(name)) {
        
        ## Testting if yvar was specified
        add <- list(...)
        if(length(add) == 0) {
            tmp <- as(x, "marrayNorm")
            ## indexing ref labelled with green
            idx <- tolower(getLabels(x, "Ref")) == "red"
            if(sum(idx) > 0)
                maM(tmp)[, idx] <- -maM(tmp)[, idx]
            ## Testing if the object is indexed
            if(dim(x)[1] < x@Layout$Nspots) {
                Mvalues <- maM(tmp)
                Avalues <- maA(tmp)
                maM(tmp) <- matrix(NA, x@Layout$Nspots, dim(x)[2])
                maA(tmp) <- matrix(NA, x@Layout$Nspots, dim(x)[2])
                indexes <- as.numeric(rownames(x@Glabels))
                for (i in 1:dim(x)[1]) {
                    maM(tmp)[indexes[i],] <- Mvalues[i,]
                    maA(tmp)[indexes[i],] <- Avalues[i,]
                }
            }
            maBoxplot(tmp, ylab="W")
        }
        else {
            tmp <- as(x, "marrayNorm")
            ## Testing if the object is indexed
            if(dim(x)[1] < x@Layout$Nspots) {
                Mvalues <- maM(tmp)
                Avalues <- maA(tmp)
                maM(tmp) <- matrix(NA, x@Layout$Nspots, dim(x)[2])
                maA(tmp) <- matrix(NA, x@Layout$Nspots, dim(x)[2])
                indexes <- as.numeric(rownames(x@Glabels))
                for (i in 1:dim(x)[1]) {
                    maM(tmp)[indexes[i],] <- Mvalues[i,]
                    maA(tmp)[indexes[i],] <- Avalues[i,]
                }
            }
            maBoxplot(tmp, ...)
        }
    }
    else {
        if(is.null(gLabelID) | is.null(sLabelID))
            stop("gLabelID and sLabelID must be specified together with name.")
        genes <- getLabels(x, gLabelID, FALSE)
        samples <- getLabels(x, sLabelID)
        x <- x[, !is.na(samples)]
        samples <- samples[!is.na(samples)]

        ## Using gSamples if it isn't NULL
        if(!is.null(gSamples))
            for(i in 1:length(gSamples))
                samples[samples %in% gSamples[[i]]] <- names(gSamples)[i]


        toPlot <- list()
        for (i in unique(gSamples))
            toPlot[[i]] <- calcW(x[genes == name, samples == i])


        boxplot(toPlot, ylab=name, col="grey", main=paste("Boxplot for", sLabelID,
        "sample label"))

    }
}


## For maigesANOVA class
boxplot.maigesANOVA <- boxplot.maiges


## For maigesDEcluster class (display boxplots by genes)
boxplot.maigesDEcluster <- function(x, name=NULL, gLabelID=NULL, sLabelID=NULL,
gSamples=NULL, ...) {

    if(is.null(gLabelID) | is.null(sLabelID))
        stop("gLabelID and sLabelID must be specified together with name.")
    genes <- getLabels(x, gLabelID, FALSE)
    if(is.null(name))
        name <- genes[1]
    samples <- getLabels(x, sLabelID)


    ## Using gSamples if it isn't NULL
    if(!is.null(gSamples))
        for(i in 1:length(gSamples))
            samples[samples %in% gSamples[[i]]] <- names(gSamples)[i]


    toPlot <- list()
    for (i in unique(samples))
        toPlot[[i]] <- x@W[genes == name, samples == i]

    boxplot(toPlot, ylab=name, col="grey", main=paste("Boxplot for", sLabelID,
    "sample label"))
    ##boxplot(toPlot, ylab=name, col=2:(length(unique(samples))+1),
    ##main=paste("Boxplot for", sLabelID, "sample label"))

}
