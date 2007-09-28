## Function to write MEV files to be loaded by TIGR MeV, from a maiges
## object. This save TDMS file
##
## Parameters: data     -> object of class maiges
##             sLabelID -> Sample label ID top be saved
##             file     -> file name, with 'txt' extension.
##
## Gustavo Henrique Esteves
## 14/05/07
##
## Version: 1.0
##


createTDMS <- function(data=NULL, sLabelID=names(data@Slabels)[1],
file="data.txt") {
    
    
    ## Doing a check
    if(is.null(data))
        stop("The data object MUST be specified!!!")
    
    
    ## Picking data and generating the table
    chips <- cbind(data@Glabels, calcW(data))
    colnames(chips) <- c(names(data@Glabels), getLabels(data, sLabelID))
    
    
    ## Writing file
    write.table(chips, file=file, sep="\t", quote=FALSE, row.names=FALSE)
    
}
