## Define methods for coerce (as) generic function
##
## Gustavo H. Esteves
## 11/05/07
##
##


## Convert maigesRaw into RGList object
setAs("maigesRaw", "RGList", function(from, to) {

    y <- new(to)

    tmpf <- matrix(NA, nrow=dim(from)[1], ncol=dim(from)[2])
    tmpb <- matrix(NA, nrow=dim(from)[1], ncol=dim(from)[2])
    idxR1 <- (tolower(from@Sdye) == "red")
    idxR2 <- (tolower(from@Rdye) == "red")
    tmpf[, idxR1] <- from@Sf[, idxR1]
    tmpb[, idxR1] <- from@Sb[, idxR1]
    tmpf[, idxR2] <- from@Rf[, idxR2]
    tmpb[, idxR2] <- from@Rb[, idxR2]
    y$R <- tmpf
    y$Rb <- tmpb

    tmpf <- matrix(NA, nrow=dim(from)[1], ncol=dim(from)[2])
    tmpb <- matrix(NA, nrow=dim(from)[1], ncol=dim(from)[2])
    idxG1 <- (tolower(from@Sdye) == "green")
    idxG2 <- (tolower(from@Rdye) == "green")
    tmpf[, idxG1] <- from@Sf[, idxG1]
    tmpb[, idxG1] <- from@Sb[, idxG1]
    tmpf[, idxG2] <- from@Rf[, idxG2]
    tmpb[, idxG2] <- from@Rb[, idxG2]
    y$G <- tmpf
    y$Gb <- tmpb

    y$weights <- from@Flag
    y$printer$ngrid.r <- from@Layout$gridR
    y$printer$ngrid.c <- from@Layout$gridC
    y$printer$nspot.r <- from@Layout$spotR
    y$printer$nspot.c <- from@Layout$spotC
    y$printer$notes <- from@Layout$Nspots
    y$genes <- from@Glabels
    y$genes$Labels <- from@Glabels[[1]]
    attr(y$genes,"notes") <- NULL
    y$genes$Sub <- NULL
    y$genes$Plate <- NULL
    y$genes$Controls <- NULL
    y$targets <- from@Slabels
    y$targets$Labels <- from@Slabels[[1]]
    y$notes <- "Converted from maigesRaw class"

    return(y)

})


## Convert RGList into maigesRaw object
setAs("RGList", "maigesRaw", function(from, to) {

    y <- new(to)

    y@Sf <- from$R
    y@Sb <- from$Rb
    y@Sdye <- rep("red", dim(from)[2])

    y@Rf <- from$G
    y@Rb <- from$Gb
    y@Rdye <- rep("green", dim(from)[2])

    y@Flag <- from$weights
    y@Layout$gridR <- from$printer$ngrid.r
    y@Layout$gridC <- from$printer$ngrid.c
    y@Layout$spotR <- from$printer$nspot.r
    y@Layout$spotC <- from$printer$nspot.c
    y@Layout$Nspots <- (from$printer$ngrid.r*from$printer$ngrid.c*
    from$printer$nspot.r*from$printer$nspot.c)
    y@Glabels <- from$genes
    y@Slabels <- from$targets
    y@BadSpots <- rep(FALSE, dim(from)[2])
    y@UseSpots <- matrix(TRUE, dim(from)[1], dim(from)[2])
    y@Notes <- "Converted from RGList class"

    return(y)

})


## Convert maigesRaw into marrayRaw object
setAs("maigesRaw", "marrayRaw", function(from, to) {

    y <- new(to)

    tmpf <- matrix(NA, nrow=dim(from)[1], ncol=dim(from)[2])
    tmpb <- matrix(NA, nrow=dim(from)[1], ncol=dim(from)[2])
    idxR1 <- (tolower(from@Sdye) == "red")
    idxR2 <- (tolower(from@Rdye) == "red")
    tmpf[, idxR1] <- from@Sf[, idxR1]
    tmpb[, idxR1] <- from@Sb[, idxR1]
    tmpf[, idxR2] <- from@Rf[, idxR2]
    tmpb[, idxR2] <- from@Rb[, idxR2]
    y@maRf <- tmpf
    y@maRb <- tmpb

    tmpf <- matrix(NA, nrow=dim(from)[1], ncol=dim(from)[2])
    tmpb <- matrix(NA, nrow=dim(from)[1], ncol=dim(from)[2])
    idxG1 <- (tolower(from@Sdye) == "green")
    idxG2 <- (tolower(from@Rdye) == "green")
    tmpf[, idxG1] <- from@Sf[, idxG1]
    tmpb[, idxG1] <- from@Sb[, idxG1]
    tmpf[, idxG2] <- from@Rf[, idxG2]
    tmpb[, idxG2] <- from@Rb[, idxG2]
    y@maGf <- tmpf
    y@maGb <- tmpb

    y@maW <- from@Flag
    y@maLayout@maNgr <- from@Layout$gridR
    y@maLayout@maNgc <- from@Layout$gridC
    y@maLayout@maNsr <- from@Layout$spotR
    y@maLayout@maNsc <- from@Layout$spotC
    y@maLayout@maNspots <- from@Layout$Nspots
    y@maLayout@maSub <- rep(TRUE, dim(from)[1])
    y@maLayout@maPlate <- factor(NA)
    y@maLayout@maControls <- factor(NA)
    y@maGnames@maInfo <- from@Glabels
    y@maGnames@maLabels <- as.character(from@Glabels[[1]])
    y@maTargets@maInfo <- from@Slabels
    y@maTargets@maLabels <- as.character(from@Slabels[[1]])
    y@maNotes <- "Converted from maigesRaw object"

    return(y)

})


## Convert marrayRaw into maigesRaw object
setAs("marrayRaw", "maigesRaw", function(from, to) {

    y <- new(to)

    y@Sf <- from@maRf
    y@Sb <- from@maRb
    y@Sdye <- rep("red", dim(from)[2])

    y@Rf <- from@maGf
    y@Rb <- from@maGb
    y@Rdye <- rep("green", dim(from)[2])

    y@Flag <- from@maW
    y@Layout$gridR <- from@maLayout@maNgr
    y@Layout$gridC <- from@maLayout@maNgc
    y@Layout$spotR <- from@maLayout@maNsr
    y@Layout$spotC <- from@maLayout@maNsc
    y@Layout$Nspots <- from@maLayout@maNspots
    y@Glabels <- from@maGnames@maInfo
    y@Slabels <- from@maTargets@maInfo
    y@BadSpots <- rep(FALSE, dim(from)[2])
    y@UseSpots <- matrix(TRUE, dim(from)[1], dim(from)[2])
    y@Notes <- "Converted from marrayRaw class"

    return(y)

})


## Convert maiges into MAList object
setAs("maiges", "MAList", function(from, to) {

    y <- new(to)
    y$M <- from@W
    ## indexing ref labelled with green
    idx <- tolower(getLabels(from, "Ref")) == "red"
    if(sum(idx) > 0)
        y$M[, idx] <- (-1)*y$M[, idx]
    y$A <- from@A
    y$weights <- matrix(numeric())
    y$printer$ngrid.r <- from@Layout$gridR
    y$printer$ngrid.c <- from@Layout$gridC
    y$printer$nspot.r <- from@Layout$spotR
    y$printer$nspot.c <- from@Layout$spotC
    y$printer$notes <- from@Layout$Nspots
    y$genes <- from@Glabels
    y$genes$Labels <- from@Glabels[[1]]
    attr(y$genes,"notes") <- NULL
    y$genes$Sub <- NULL
    y$genes$Plate <- NULL
    y$genes$Controls <- NULL
    y$targets <- from@Slabels
    y$targets$Labels <- from@Slabels[[1]]
    y$notes <- "Converted from maiges class"

    return(y)

})


## Convert MAList into maiges object
setAs("MAList", "maiges", function(from, to) {

    y <- new(to)

    y@W <- from$M
    y@A <- from$A

    y@Layout$gridR <- from$printer$ngrid.r
    y@Layout$gridC <- from$printer$ngrid.c
    y@Layout$spotR <- from$printer$nspot.r
    y@Layout$spotC <- from$printer$nspot.c
    y@Layout$Nspots <- (from$printer$ngrid.r*from$printer$ngrid.c*
    from$printer$nspot.r*from$printer$nspot.c)
    y@Glabels <- from$genes
    y@Slabels <- from$targets
    y@BadSpots <- rep(FALSE, dim(from)[2])
    y@Notes <- "Converted from MAList class"

    return(y)

})


## Convert maiges into marrayNorm
setAs("maiges", "marrayNorm", function(from, to) {

    y <- new(to)

    y@maM <- from@W
    ## indexing ref labelled with green
    idx <- tolower(getLabels(from, "Ref")) == "red"
    if(sum(idx) > 0)
        y@maM[, idx] <- (-1)*y@maM[, idx]
    y@maA <- from@A
    y@maW <- matrix(numeric())
    y@maLayout@maNgr <- from@Layout$gridR
    y@maLayout@maNgc <- from@Layout$gridC
    y@maLayout@maNsr <- from@Layout$spotR
    y@maLayout@maNsc <- from@Layout$spotC
    y@maLayout@maNspots <- from@Layout$Nspots
    y@maLayout@maSub <- rep(TRUE, dim(from)[1])
    y@maLayout@maPlate <- factor(NA)
    y@maLayout@maControls <- factor(NA)
    y@maGnames@maInfo <- from@Glabels
    y@maGnames@maLabels <- as.character(from@Glabels[[1]])
    y@maTargets@maInfo <- from@Slabels
    y@maTargets@maLabels <- as.character(from@Slabels[[1]])
    y@maNotes <- "Converted from maiges object"

    return(y)

})


## Convert marrayNorm into maiges object
setAs("marrayNorm", "maiges", function(from, to) {

    y <- new(to)

    y@W <- from@maM
    y@A <- from@maA

    y@Layout$gridR <- from@maLayout@maNgr
    y@Layout$gridC <- from@maLayout@maNgc
    y@Layout$spotR <- from@maLayout@maNsr
    y@Layout$spotC <- from@maLayout@maNsc
    y@Layout$Nspots <- from@maLayout@maNspots
    y@Glabels <- from@maGnames@maInfo
    y@Slabels <- from@maTargets@maInfo
    y@BadSpots <- rep(FALSE, dim(from)[2])
    y@Notes <- "Converted from marrayNorm class"

    return(y)

})
