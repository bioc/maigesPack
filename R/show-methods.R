## Define methods for show generic functions
##
## Gustavo H. Esteves
## 10/05/07
##
##

## For maigesPreRaw class
setMethod("show", "maigesPreRaw", function(object) print(object))


## For maigesRaw class
setMethod("show", "maigesRaw", function(object) print(object))


## For maiges class
setMethod("show", "maiges", function(object) print(object))


## For maigesANOVA class
setMethod("show", "maigesANOVA", function(object) print(object))


## For maigesDE class
setMethod("show", "maigesDE", function(object) print(object))


## For maigesDEcluster class
setMethod("show", "maigesDEcluster", function(object) print(object))


## For maigesClass class
setMethod("show", "maigesClass", function(object) print(object))


## For maigesRelNetM class
setMethod("show", "maigesRelNetM", function(object) print(object))


## For maigesRelNetB class
setMethod("show", "maigesRelNetB", function(object) print(object))


## For maigesActMod class
setMethod("show", "maigesActMod", function(object) print(object))


## For maigesActNet class
setMethod("show", "maigesActNet", function(object) print(object))

