## Define Classes for maigesPack package
##
## Gustavo H. Esteves
## 17/05/07
##
## Version: 1.0
##



## Defining maigesPreRaw class: Class for describing Object with raw intensity
## values and information about genes and samples used in the data. Here, user
## can put all fields  of the original data files.
##
setClass("maigesPreRaw", representation(Data="list", GeneGrps="list",
Paths="list", Layout="list", Glabels="data.frame", Slabels="data.frame",
BadSpots="logical", Notes="character", Date="character", V.info="list"))


## Defining maigesRaw class: Class for describing Object with raw intensity
## values and information about genes and samples used in the data
##
setClass("maigesRaw", representation(Sf="matrix", Sb="matrix", Sdye="vector",
Rf="matrix", Rb="matrix", Rdye="vector", Flag="matrix", BadSpots="logical",
UseSpots="matrix", GeneGrps="matrix", Paths="list", Layout="list",
Glabels="data.frame", Slabels="data.frame", Notes="character", Date="character",
V.info="list"))


## Defining maiges class: Class for describing Object with (normalized) values
## and information about genes and samples used in the data
##
setClass("maiges", representation(W="matrix", A="matrix", SD="matrix",
IC1="matrix", IC2="matrix", BadSpots="logical", UseSpots="matrix",
GeneGrps="matrix", Paths="list", Layout="list", Glabels="data.frame",
Slabels="data.frame", Notes="character", Date="character", V.info="list"))


## Defining maigesANOVA class: Class extending maiges class containing design
## and contrat matrices
##
setClass("maigesANOVA", representation("maiges", Dmatrix="matrix",
Cmatrix="matrix"))


## Defining maigesDE class: Class for describing results of differential
## expression analyses
##
setClass("maigesDE", representation(GeneInfo="data.frame",
SampleInfo="data.frame", fold="matrix", stat="matrix", p.value="matrix",
factors="character", test="character", Date="character", V.info="list"))


## Defining maigesDEcluster class: Class for describing results of differential
## expression analyses. Extends the previous one, containing the W matrix
##
setClass("maigesDEcluster", representation("maigesDE", W="matrix"))


## Defining maigesClass class: Class for describing results of discrimination
## analysis.
##
setClass("maigesClass", representation(W="matrix", CV="vector", SVD="vector",
cliques="matrix", cliques.idx="matrix", method="character", Date="character",
V.info="list"))


## Defining maigesRelNetM class: Class for describing results of Relevance
## networks analysis, using our model adapted the Butte's one
##
setClass("maigesRelNetM", representation(W="matrix", Corr1="matrix",
Corr2="matrix", DifP="matrix", types="vector", Slabel="character",
Date="character", V.info="list"))


## Defining maigesRelNetB class: Class for describing results of the original
## Butte's method for constructing relevance networks
##
setClass("maigesRelNetB", representation(W="matrix", Corr="matrix",
Pval="matrix", maxB="matrix", type="character", Slabel="character",
Date="character", V.info="list"))


## Defining maigesActMod class: Class for describing results of functional
## classification of gene groups (Segal's method)
##
setClass("maigesActMod", representation(modBySamp="matrix", modByCond="matrix",
globalScore="list", tissueScore="list", Date="character", V.info="list"))


## Defining maigesActNet class: Class for describing results of our statistical
## model for functional classification of gene regulatory networks
##
setClass("maigesActNet", representation(scores="matrix", Pvalues="matrix",
Date="character", V.info="list"))


## Guaranting the loading of methods package to use classes and methods defined
## in it
##
.onLoad <- function(lib, pkg) require("methods")
