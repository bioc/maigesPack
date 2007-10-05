## Define methods for indexing generic methods
##
## Gustavo H. Esteves
## 26/08/07
##
##


## For maigesRaw class
setMethod("[", "maigesPreRaw", structure(function(x, i, j, ..., drop=FALSE) {
    newx <- x
    if (missing(j))
        j <- TRUE
    if (missing(i))
        i <- TRUE
    if (length(x@Data) != 0)
        for (k in names(x@Data))
            newx@Data[[k]] <- x@Data[[k]][i, j, drop=FALSE]
    gastro@Layout <- gastro@Layout
    gastro@Glabels <- gastro@Glabels
    gastro@Slabels <- gastro@Slabels
    gastro@Notes <- gastro@Notes

    if (length(x@BadSpots) != 0)
        newx@BadSpots <- x@BadSpots[i]

    
    newx@Glabels <- x@Glabels[i, ]
    newx@Slabels <- x@Slabels[j, ]
    newx@Date <- date()

    return(newx)
}, ## Close the function definition
class=structure("MethodDefinition", package="methods"),
target=structure("maigesPreRaw", .Names="x", class=structure("signature",
package="methods")), defined=structure("maigesPreRaw", .Names="x",
class=structure("signature", package="methods")))) ## Close setMethod


## For maigesRaw class
setMethod("[", "maigesRaw", structure(function(x, i, j, ..., drop=FALSE) {
    newx <- x
    if (missing(j))
        j <- TRUE
    if (missing(i))
        i <- TRUE
    if (length(x@Sf) != 0)
        newx@Sf <- x@Sf[i, j, drop=FALSE]
    if (length(x@Sb) != 0)
        newx@Sb <- x@Sb[i, j, drop=FALSE]
    if (length(x@Sdye) != 0)
        newx@Sdye <- x@Sdye[j]
    if (length(x@Rf) != 0)
        newx@Rf <- x@Rf[i, j, drop=FALSE]
    if (length(x@Rb) != 0)
        newx@Rb <- x@Rb[i, j, drop=FALSE]
    if (length(x@Rdye) != 0)
        newx@Rdye <- x@Rdye[j]
    if (length(x@Flag) != 0)
        newx@Flag <- x@Flag[i, j, drop=FALSE]
    if (length(x@UseSpots) != 0)
        newx@UseSpots <- x@UseSpots[i, j, drop=FALSE]
    if (length(x@BadSpots) != 0)
        newx@BadSpots <- x@BadSpots[i]
    if (length(x@GeneGrps) != 0)
        newx@GeneGrps <- x@GeneGrps[i, , drop=FALSE]

    if (length(x@Paths) != 0) {
        gUsed <- getLabels(x, x@Paths$Glabel, FALSE)[i]
        for(k in 2:length(x@Paths)) {
            idx <- (nodes(x@Paths[[k]]) %in% gUsed)
            newx@Paths[[k]] <- subGraph(nodes(x@Paths[[k]])[idx], x@Paths[[k]])
        }
    }

    newx@Glabels <- x@Glabels[i, ]
    newx@Slabels <- x@Slabels[j, ]
    newx@Date <- date()

    return(newx)
}, ## Close the function definition
class=structure("MethodDefinition", package="methods"),
target=structure("maigesRaw", .Names="x", class=structure("signature",
package="methods")), defined=structure("maigesRaw", .Names="x",
class=structure("signature", package="methods")))) ## Close setMethod


## For maiges class
setMethod("[", "maiges", structure(function(x, i, j, ..., drop=FALSE) {
    newx <- x
    if (missing(j))
        j <- TRUE
    if (missing(i))
        i <- TRUE
    if (length(x@W) != 0)
        newx@W <- x@W[i, j, drop=FALSE]
    if (length(x@A) != 0)
        newx@A <- x@A[i, j, drop=FALSE]
    if (length(x@SD) != 0)
        newx@SD <- x@SD[i, j, drop=FALSE]
    if (length(x@IC1) != 0)
        newx@IC1 <- x@IC1[i, j, drop=FALSE]
    if (length(x@IC2) != 0)
        newx@IC2 <- x@IC2[i, j, drop=FALSE]
    if (length(x@BadSpots) != 0)
        newx@BadSpots <- x@BadSpots[i]
    if (length(x@UseSpots) != 0)
        newx@UseSpots <- x@UseSpots[i, j, drop=FALSE]
    if (length(x@GeneGrps) != 0)
        newx@GeneGrps <- x@GeneGrps[i, , drop=FALSE]

    if (length(x@Paths) != 0) {
        gUsed <- getLabels(x, x@Paths$Glabel, FALSE)[i]
        for(k in 2:length(x@Paths)) {
            idx <- nodes(x@Paths[[k]]) %in% gUsed
            newx@Paths[[k]] <- subGraph(nodes(x@Paths[[k]])[idx], x@Paths[[k]])
        }
    }

    newx@Glabels <- x@Glabels[i, ]
    newx@Slabels <- x@Slabels[j, ]
    newx@Date <- date()

    return(newx)
}, ## Close function
class=structure("MethodDefinition", package="methods"),
target=structure("maiges", .Names="x", class=structure("signature",
package="methods")), defined=structure("maiges", .Names="x",
class=structure("signature", package="methods")))) ## close setMethod


## For maigesANOVA class
setMethod("[", "maigesANOVA", structure(function(x, i, j, ..., drop=FALSE) {
    newx <- x
    if (missing(j))
        j <- TRUE
    if (missing(i))
        i <- TRUE
    if (length(x@W) != 0)
        newx@W <- x@W[i, j, drop=FALSE]
    if (length(x@A) != 0)
        newx@A <- x@A[i, j, drop=FALSE]
    if (length(x@SD) != 0)
        newx@SD <- x@SD[i, j, drop=FALSE]
    if (length(x@IC1) != 0)
        newx@IC1 <- x@IC1[i, j, drop=FALSE]
    if (length(x@IC2) != 0)
        newx@IC2 <- x@IC2[i, j, drop=FALSE]
    if (length(x@BadSpots) != 0)
        newx@BadSpots <- x@BadSpots[i]
    if (length(x@GeneGrps) != 0)
        newx@GeneGrps = x@GeneGrps[i, , drop=FALSE]
    if (length(x@Dmatrix) != 0)
        newx@Dmatrix <- x@Dmatrix[j, , drop=FALSE]

    if (length(x@Paths) != 0) {
        gUsed <- getLabels(x, x@Paths$Glabel, FALSE)[i]
        for(k in 2:length(x@Paths)) {
            idx <- nodes(x@Paths[[k]]) %in% gUsed
            newx@Paths[[k]] <- subGraph(nodes(x@Paths[[k]])[idx], x@Paths[[k]])
        }
    }

    newx@Glabels <- x@Glabels[i, ]
    newx@Slabels <- x@Slabels[j, ]
    newx@Date <- date()

    return(newx)
}, ## Close function
class=structure("MethodDefinition", package="methods"),
target=structure("maiges", .Names="x", class=structure("signature",
package="methods")), defined=structure("maiges", .Names="x",
class=structure("signature", package="methods")))) ## close setMethod
