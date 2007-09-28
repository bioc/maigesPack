## Definition of generic functions
##
## Gustavo H. Esteves
## 10/05/07
##
## Version: 1.0
##


## Generic to extract A values from microarray objects
calcA <- function(object, ...)
    UseMethod("calcA")


## Generic to extract W values from microarray objects
calcW <- function(object, ...)
    UseMethod("calcW")


## Generic to pick information about labels from microarray objects
getLabels <- function(obj=NULL, labelID=NULL, sLabel=TRUE)
    UseMethod("getLabels")


## Generic to convert relevance networks in TGF files
relNet2TGF <- function(...)
    UseMethod("relNet2TGF")
