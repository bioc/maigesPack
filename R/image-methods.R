## Define methods for image generic function
##
## Gustavo H. Esteves
## 27/05/07
##
## Version: 1.1
##


## For maigesRaw class (display images of slides also using marray)
image.maigesRaw <- function(x, ...) {
    

  
    ## Testting if xvar was specified
    add <- list(...)
    if(length(add) == 0) {
        tmp <- as(x, "marrayRaw")
        tmp <- as(tmp, "marrayNorm")
        ## indexing ref labelled with green
        idx <- tolower(getLabels(x, "Ref")) == "red"
        if(sum(idx) > 0)
            maM(tmp)[, idx] <- -maM(tmp)[, idx]
        title <- paste(maLabels(maTargets(tmp))[1], ": image of W", sep="")
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
        maImage(tmp, main=title, subset=indexes)
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
        maImage(tmp, subset=indexes, ...)
        
    }
}


## For maiges class (display images for slides, also using marray)
image.maiges <- function(x, ...) {
    
    ## Testting if xvar was specified
    add <- list(...)
    if(length(add) == 0) {
        tmp <- as(x, "marrayNorm")
        ## indexing ref labelled with green
        idx <- tolower(getLabels(x, "Ref")) == "red"
        if(sum(idx) > 0)
            maM(tmp)[, idx] <- -maM(tmp)[, idx]
        title <- paste(maLabels(maTargets(tmp))[1], ": image of W", sep="")
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
        maImage(tmp, subset=indexes, main=title)
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
        maImage(tmp, subset=indexes, ...)
    }
}


## For maigesANOVA class
image.maigesANOVA <- image.maiges


## For maigesRelNetM class (heatmaps of the correlation coefficients)
image.maigesRelNetM <- function(x=NULL, names=NULL, ...) {

    if(is.null(names))
        names <- c(x@types, "Significance of differences")

    ## Defining the initial display
    par(mfrow=c(1,3))

    ## Doing the first heatmap

    limite <- max(abs(range(x@Corr1)))
    limite <- c(-limite, limite)
    
    par(mar=c(5, 9, 9, 1))

    image(1:nrow(x@Corr1), 1:ncol(x@Corr1), t(x@Corr1[nrow(x@Corr1):1,]),
    col=maigesPack:::greenRed(), axes=FALSE, xlab=names[1], ylab="",
    zlim=limite, ...)

    axis(2, 1:nrow(x@Corr1), las=2, labels=rev(rownames(x@Corr1)),
    cex.axis=1.6, col.axis=1)

    axis(3, 1:ncol(x@Corr1), las=2, labels=colnames(x@Corr1),
    cex.axis=1.6, col.axis=1)

    ## Doing the second heatmap

    limite <- max(abs(range(x@Corr2)))
    limite <- c(-limite, limite)
    
    par(mar=c(5, 9, 9, 1))

    image(1:nrow(x@Corr2), 1:ncol(x@Corr2), t(x@Corr2[nrow(x@Corr2):1,]),
    col=maigesPack:::greenRed(), axes=FALSE, xlab=names[2], ylab="",
    zlim=limite, ...)

    axis(2, 1:nrow(x@Corr2), las=2, labels=rev(rownames(x@Corr2)),
    cex.axis=1.6, col.axis=1)

    axis(3, 1:ncol(x@Corr2), las=2, labels=colnames(x@Corr2),
    cex.axis=1.6, col.axis=1)

    ## Doing the third heatmap

    limite <- max(abs(range(-log(x@DifP))))
    limite <- c(0, limite)
    
    par(mar = c(5, 9, 9, 1))

    image(1:nrow(x@DifP), 1:ncol(x@DifP), t(-log(x@DifP)[nrow(x@DifP):1,]),
    col=maigesPack:::blackBlue(), axes=FALSE, xlab=names[3],
    ylab="", zlim=limite, ...)

    axis(2, 1:nrow(x@DifP), las=2, labels=rev(rownames(x@DifP)),
    cex.axis=1.6, col.axis=1)

    axis(3, 1:ncol(x@DifP), las=2, labels=colnames(x@DifP),
    cex.axis=1.6, col.axis=1)

}


## For maigesRelNetB class (heatmaps of the correlation coefficients)
image.maigesRelNetB <- function(x=NULL, name=NULL, ...) {

    if(is.null(name))
        name <- x@type

    ## Doing the heatmap

    limite <- max(abs(range(x@Corr)))
    limite <- c(-limite, limite)
    
    par(mar=c(5, 9, 9, 1))

    image(1:nrow(x@Corr), 1:ncol(x@Corr), t(x@Corr[nrow(x@Corr):1,]),
    col=maigesPack:::greenRed(), axes=FALSE, xlab=name, ylab="", zlim=limite,
    ...)

    axis(2, 1:nrow(x@Corr), las=2, labels=rev(rownames(x@Corr)),
    cex.axis=1.6, col.axis=1)

    axis(3, 1:ncol(x@Corr), las=2, labels=colnames(x@Corr),
    cex.axis=1.6, col.axis=1)

}


## For maigesActMod class (same result as plot.maigesActMod)
image.maigesActMod <- function(x, type=c("S", "C")[2],
keepEmpty=FALSE, ...) {

    ## Making some basic initial tests...
    if(is.null(x))
        stop("You MUST specify an object generated by activeMod function.")
    if(!is.element(type, c("S", "C")))
        stop("You must be 'C' or 'S'.")

    if(type == "S") {
        if(keepEmpty)
            table <- x@modBySamp
        else {
            idx <- apply(x@modBySamp != 0, 2, sum, na.rm=TRUE) != 0
            if(sum(idx) < 2)
                stop("Less than 2 elements present significant results!")
            table <- x@modBySamp[, idx]
        }
        limite <- max(abs(range(table, na.rm=TRUE)))
        limite <- c(-limite, limite)
        idx1 <- order(rownames(table))
        idx2 <- order.dendrogram(as.dendrogram(hclust(dist(t(table)))))

        heatmap(table[idx1, idx2], scale="none", col=maigesPack:::greenRed(),
        zlim=limite, Rowv=NA, Colv=NA, ...)
    }
    else if(type == "C") {
        if(keepEmpty)
            table <- x@modByCond
        else {
            idx <- apply(x@modByCond != 0, 2, sum, na.rm=TRUE) != 0
            if(sum(idx) < 2)
                stop("Less than 2 elements present significant results!")
            table <- x@modByCond[, idx]
        }

        limite <- max(abs(range(table, na.rm=TRUE)))
        limite <- c(-limite, limite)
        idx1 <- order(rownames(table))
        idx2 <- order.dendrogram(as.dendrogram(hclust(dist(t(table)))))

        heatmap(table[idx1, idx2], scale="none", col=maigesPack:::greenRed(),
        zlim=limite, Rowv=NA, Colv=NA, ...)
    }
}


## For maigesActNet class (heatmap of the significative results)
image.maigesActNet <- function(x, type=c("score","p-value")[1], ...) {

    ## Making some basic initial tests...
    if(is.null(x))
        stop("You MUST specify an object generated by activeNet function.")
    if(!is.element(type, c("score", "p-value")))
        stop("You must be 'score' or 'p-value'.")

    if(type == "score") {

        limite <- max(x@scores, na.rm=TRUE)
        limite <- c(0, limite)
        idx <- order.dendrogram(as.dendrogram(hclust(dist(t(x@scores)))))

        heatmap(x@scores[, idx], scale="none", col=maigesPack:::blackBlue(),
        zlim=limite, Rowv=NA, Colv=NA, ...)
    }
    else {
        
        limite = max(-log10(x@Pvalues), na.rm=TRUE)
        limite = c(0, limite)
        
        idx <- order.dendrogram(as.dendrogram(hclust(dist(t(x@Pvalues)))))
        
        heatmap(-log10(x@Pvalues)[, idx], scale="none",
        col=maigesPack:::blackBlue(), zlim=limite, Rowv=NA, Colv=NA, ...)
    }
}
