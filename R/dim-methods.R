## Define methods for dim generic function
##
## Gustavo H. Esteves
## 10/05/07
##
## Version: 1.0
##


## For maigesPreRaw class
dim.maigesPreRaw <- function(x)
    c(length(x@Glabels[,1]), length(x@Slabels[,1]))


## For maigesRaw class
dim.maigesRaw <- function(x)
    dim(x@Sf)


## For maiges class
dim.maiges <- function(x)
    dim(x@W)


## For maigesANOVA class
dim.maigesANOVA <- function(x)
    dim(x@W)

