## Function to write tables of DE genes.
##
## Parameters: deComp    -> object pf maigesDE class
##             dir       -> char specifying the directory to save the files
##             filenames -> char vector with file names to be saved, if NULL
##                          default names are given.
##             dataID    -> char giving an identifier for the dataset
##             type      -> type of file to be saved. May be 'HTML' (default)
##                          or 'CSV'
##             geneID    -> char giving the Id of gene name label in the
##                          dataset
##             hsID      -> char giving the Id of unigene (unique cluster)
##                          label in the dataset
##             gbID      -> char giving the Id of Genbank entry label in the
##                          dataset
##             annotID   -> char giving the Id of gene annotation label in the
##                          dataset
##             logFold   -> logical. Must fold be in log2?
##             genes     -> character vector specifying the genes to be saved
##             adjP      -> type of p-value adjustment. May be "Bonferroni",
##                          "Holm", "Hochberg", "SidakSS", "SidakSD", "BH",
##                          "BY" or "none". Defaults to "none"
##             sort      -> what field to sort, may be "p.value" (default),
##                          "fold" or "statistic".
##             nDEgenes  -> number of differentially genes to be save in the
##                          file. Defaults to NULL for all genes. If a number
##                          in (0,1), a cuttoff p.value is used.
##
## Gustavo Esteves
## 26/06/12
##
##


tablesDE <- function(deComp=NULL, dir="./", filenames=NULL, dataID="Someone's",
type=c("HTML","CSV")[1], geneID="GeneName", hsID="ClusterId", gbID="GeneId",
annotID="Annot", genes=NULL, logFold=TRUE, adjP="none", sort="p.value",
nDEgenes=NULL) {
  
  
  ## Defining a new repositoty to NCBI gene name
  repofun <- function(ids){
    out <- paste("http://www.ncbi.nlm.nih.gov/gene?term=", ids, sep = "")
    out
  }
  require(annotate)
  setRepository("gene", repofun)
  
    
    ## Doing some initial tests...
    if(is.null(deComp))
        stop("You MUST specify the deComp, a maigesDE object.")
    
    if(!is.element(type, c("HTML","CSV")))
        stop("type must be 'HTML' or 'CSV'.")
    
    if(is.null(geneID))
        stop("You must specify geneID parameter.")
    
    if(!is.element(adjP, c("Bonferroni", "Holm", "Hochberg", "SidakSS",
    "SidakSD", "BH", "BY", "none")))
        stop("Incorrect value for parameter adjP see help pages for
        or 'mt.rawp2adjp'")
    
    if(!is.element(sort, c("p.value", "fold", "statistic")))
        stop("Sort must be 'p.value', 'fold' or 'statistic'.")
    
    
    
    ## tranforming nDEgenes in the number of genes if it is null
    if(is.null(nDEgenes))
        nDEgenes <- nrow(deComp@W)
    
    
    ## Adjusting dir name for function won't crash
    tmp <- strsplit(dir, NULL)
    if (tmp[[1]][length(tmp[[1]])] != "/")
        dir <- paste(dir, "/", sep="")
    
    deCompOrd <- list()
    if(!is.null(colnames(deComp@p.value)))
        namesDEcomp <- colnames(deComp@p.value)
    else
        namesDEcomp <- ""
    
    nComp <- length(namesDEcomp)
    
    
    ## Getting gene information
    geneInfo <- NULL
    geneInfo <- cbind(geneInfo, Name=getLabels(deComp, geneID, FALSE))
    if(!is.null(gbID))
        geneInfo <- cbind(geneInfo, Genbank=getLabels(deComp, gbID, FALSE))
    
    if(!is.null(hsID))
        geneInfo <- cbind(geneInfo, Unigene=getLabels(deComp, hsID, FALSE))
    
    if(!is.null(annotID))
        geneInfo <- cbind(geneInfo, Annotation=getLabels(deComp, annotID,
        FALSE))
    
    
    nInfo <- length(geneInfo[1, ])
    
    
    ## Adjusting p-values if specified by the user.
    if(adjP != "none") {
        for(i in 1:nComp) {
            tmp1 <- multtest::mt.rawp2adjp(deComp@p.value[, i], proc=adjP)
            deComp@p.value[, i] <- tmp1$adjp[order(tmp1$index), 2]
        }
    }
    
    
    ## log-transforming (or not) fold and p-values
    if(!logFold)
        for(i in 1:nComp)
            deComp@fold[, i] <- 2^deComp@fold[, i]
    
    
    ## ordering data...
    if(sort == "p.value") {
        if(!is.null(colnames(deComp@p.value)))
            for(i in 1:nComp) {
                auxIndex <- sort(deComp@p.value[, i], index.return=TRUE)$ix

                deCompOrd[[i]] <- cbind(geneInfo[auxIndex, ],
                Fold=deComp@fold[auxIndex, i],
                Statistic=deComp@stat[auxIndex, i],
                P.value=deComp@p.value[auxIndex, i])
                
                if(type == "HTML")
                    deCompOrd[[i]][, (nInfo+1):(nInfo+3)] <- 
                    round(as.numeric(deCompOrd[[i]][, (nInfo+1):(nInfo+3)]), 4)
                
            }
        else {

            auxIndex <- sort(deComp@p.value, index.return=TRUE)$ix
            deCompOrd[[1]] <- cbind(geneInfo[auxIndex, ],
            Statistic=deComp@stat[auxIndex, 1],
            P.value=deComp@p.value[auxIndex, 1])
            
            if(type == "HTML")
                deCompOrd[[1]][, (nInfo+1):(nInfo+2)] <-
                round(as.numeric(deCompOrd[[1]][, (nInfo+1):(nInfo+2)]), 4)
            
        }
    }
    else if (sort == "fold") {
        if(!is.null(colnames(deComp@fold)))
            for(i in 1:nComp) {

                auxIndex <- sort(abs(deComp@fold[, i]), index.return=TRUE,
                decreasing=TRUE)$ix
                
                deCompOrd[[i]] <- cbind(geneInfo[auxIndex, ],
                Fold=deComp@fold[auxIndex, i], Statistic=deComp@stat[auxIndex,i],
                P.value=deComp@p.value[auxIndex, i])
                
                if(type == "HTML")
                    deCompOrd[[i]][, (nInfo+1):(nInfo+3)] <- 
                    round(as.numeric(deCompOrd[[i]][, (nInfo+1):(nInfo+3)]), 4)
                
            }
        else {
            auxIndex <- sort(abs(deComp@fold), index.return=TRUE,
            decreasing=TRUE)$ix
            
            deCompOrd[[1]] <- cbind(geneInfo[auxIndex, ],
            Statistic=deComp@stat[auxIndex, 1],
            P.value=deComp@p.value[auxIndex, 1])
            
            if(type == "HTML")
                deCompOrd[[1]][, (nInfo+1):(nInfo+2)] <- 
                round(as.numeric(deCompOrd[[1]][, (nInfo+1):(nInfo+2)]), 4)
            
        }
    }
    else if (sort == "statistic") {
        if(!is.null(colnames(deComp@stat)))
            for(i in 1:nComp) {
                auxIndex <- sort(abs(deComp@stat[, i]), index.return=TRUE,
                decreasing=TRUE)$ix
                
                deCompOrd[[i]] <- cbind(geneInfo[auxIndex, ],
                Fold=deComp@fold[auxIndex, i],
                Statistic=deComp@stat[auxIndex, i],
                P.value=deComp@p.value[auxIndex, i])
                
                if(type == "HTML")
                    deCompOrd[[i]][, (nInfo+1):(nInfo+3)] <- 
                    round(as.numeric(deCompOrd[[i]][, (nInfo+1):(nInfo+3)]), 4)
                
            }
        else {
            auxIndex <- sort(abs(deComp@stat), index.return=TRUE,
            decreasing=TRUE)$ix
            
            deCompOrd[[1]] <- cbind(geneInfo[auxIndex, ],
            Statistic=deComp@stat[auxIndex, 1],
            P.value=deComp@p.value[auxIndex, 1])
            
            if(type == "HTML")
                deCompOrd[[1]][, (nInfo+1):(nInfo+2)] <- 
                round(as.numeric(deCompOrd[[1]][, (nInfo+1):(nInfo+2)]), 4)
            
        }
    }
    
    
    
    ##  Construct HTML (or CSV) tables of diferentially expressed genes
    head <- colnames(deCompOrd[[1]])
    
    nn <- length(head)
    if(!is.null(gbID) & !is.null(hsID))
        repository=list("gene","gb","ug")
    else if(!is.null(gbID))
        repository=list("gene","gb")
    else if(!is.null(hsID))
        repository=list("gene","ug")
    else
        repository=list("gene")
    
    for(i in 1:nComp) {
        
        ## Getting indexes of the gene list specified or more DE genes
        if(!is.null(genes))
            idx <- which(deCompOrd[[i]][, 1] %in% genes)
        else if(nDEgenes > 0 & nDEgenes < 1) {
            tmp <- colnames(deCompOrd[[i]]) == "P.value"
            idx <- which(as.numeric(deCompOrd[[i]][, tmp]) <= nDEgenes)
        }
        else {
            tmp1 <- colnames(deCompOrd[[i]]) == "P.value"
            tmp2 <- sort(deCompOrd[[i]][, tmp1])[nDEgenes]
            idx <- which(deCompOrd[[i]][, tmp1] <= tmp2)
        }
        
        if(type == "HTML") {

            require("annotate")
            
            if (is.null(filenames)) {
                
                filename <- paste(dir, deComp@test, "_", namesDEcomp[i], "_",
                format(Sys.Date(), "%d%m%Y"), ".html", sep="")
                
            }
            else
                filename <- paste(dir, filenames[i])
            
            title <- paste(namesDEcomp[i], ". Method: ", deComp@test, sep="")
            
            if(!is.null(gbID) & !is.null(hsID))
                genelist<-as.list(as.data.frame(matrix(deCompOrd[[i]][idx, 1:3],
                length(idx), 3)))
            
            else if(!is.null(gbID) | !is.null(hsID))
                genelist<-as.list(as.data.frame(matrix(deCompOrd[[i]][idx, 1:2],
                length(idx), 2)))
            
            else
                genelist<-as.list(as.data.frame(matrix(deCompOrd[[i]][idx, 1],
                length(idx), 1)))
            
            othernames <- as.list(as.data.frame(deCompOrd[[i]][idx,
            (length(genelist)+1):nn]))
            
            annotate::htmlpage(genelist, filename, title, othernames, head,
            repository=repository)
            
        }
        else {
            if(is.null(filenames)) {
                
                filename <- paste(dir, namesDEcomp[i], "_", format(Sys.Date(),
               "%d%m%Y"), ".csv", sep="")
                
            }
            else
                filename <- paste(dir, filenames[i])
            
            write.table(deCompOrd[[i]][idx,], file=filename, row.names=FALSE,
            quote=FALSE, sep="\t")
            
        }
    }
    
}
