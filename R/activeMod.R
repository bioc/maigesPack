## Function to search by modules (groups) of genes with similar profiles of
## expression
##
## Parameters: data      -> object of maiges class
##             gNameID   -> Identification of gene name label ID.
##             sLabelID  -> Label identification to be used
##             samples   -> a list with character vectors specifying
##                          the groups that must be compared
##             usePaths  -> logical specifying if the Paths must also
##                          be used, defaults to FALSE
##             cutExp    -> real number specifying the cuttoff for
##                          expression levels (to discretize the expression)
##             cutPhiper -> p-value cutoff to select significant gene groups
##             adjP      -> type of p-value adjustment. May be
##                          "Bonferroni", "Holm", "Hochberg",
##                          "SidakSS", "SidakSD", "BH",
##                          "BY" or "none". Defaults to "none"
##
##
##
##
## Gustavo H. Esteves
## 27/05/07
##
##


activeMod <- function(data=NULL, gNameID="GeneName", samples=NULL,
usePaths=FALSE, sLabelID="Classification", adjP="none", cutExp=1,
cutPhiper=0.05) {


    ## Making some basic initial tests...
    if(is.null(data))
        stop("You MUST specify a valid 'data' object.")
    if(cutExp < 0)
        stop("You must specify a cutExp parameter greater than 0.")
    if(!is.element(adjP, c("Bonferroni", "Holm", "Hochberg", "SidakSS",
    "SidakSD", "BH", "BY", "none")))
        stop("Incorrect value for par. adjP see help pages for 'mt.rawp2adjp'")


    ## Geting gene and sample labels
    if(is.null(samples))
        samples <- getLabels(data, sLabelID)
    else {
        tmp <- getLabels(data, sLabelID)
        idxF <- NULL
        for (i in 1:length(samples)) {
            idxTmp <- which(is.element(tmp, samples[[i]]))
            tmp[idxTmp] <- names(samples)[i]
            idxF <- c(idxF, idxTmp)
        }
        samples <- tmp[idxF]

        ## Picking only the samples specified
        data <- data[, idxF]

    }


    ## Joining paths to GeneGrps, if specified
    if(usePaths) {
        oldNames <- colnames(data@GeneGrps)
        for(i in 2:length(data@Paths)) {
            idxPath <- getLabels(data, data@Paths[[1]], FALSE) %in%
            nodes(data@Paths[[i]])

            data@GeneGrps <- cbind(data@GeneGrps, idxPath)
        }
        colnames(data@GeneGrps) <- c(oldNames, names(data@Paths)[-1])
    }


    ## Removing NAs from table
    goodIdx <- (apply(is.na(calcW(data)), 1, sum) == 0)
    data <- data[goodIdx, ]


    ## Removing eventual groups wiht less than one gene...
    idxTmp <- (apply(data@GeneGrps, 2, sum) >= 1)
    data@GeneGrps <- data@GeneGrps[, idxTmp]


    ## Geting colnames of GeneGrps slot as modules
    modules <- colnames(data@GeneGrps)


    table <- calcW(data)
    genes <- getLabels(data, gNameID, FALSE)
    genes[data@BadSpots] <- paste(genes[data@BadSpots], "(*)")

    nGenes <- length(table[, 1]) ## Number of Genes
    nArrays <- length(table[1, ]) ## Number of Arrays
    nModules <- length(modules) ## Number of modules (gene groups)

    for (i in 1:nGenes)
        table[i, ] <- table[i, ]-mean(table[i, ])


    ## Defining the number of genes induced and repressed in each array...
    ## defining the matrix of induced (repressed) genes
    tableIR <- matrix(0, nrow=nGenes, ncol=nArrays)
    for (i in 1:nGenes) {
        idxTmp1 <- table[i, ] <= -cutExp
        idxTmp2 <- table[i, ] >= cutExp
        tableIR[i, idxTmp1] <- -1
        tableIR[i, idxTmp2] <- 1
    }


    ## Counting the number of genes induced in each module for each
    ## array
    arrayVSmoduleIR <- matrix(0, nrow=nArrays, ncol=nModules)
    for (j in 1:nArrays) {
        inducedA <- sum(tableIR[, j] == 1) ## induced at array j
        repressedA <- sum(tableIR[, j] == -1) ## repressed at array j
        pValueI <- NULL
        pValueR <- NULL
        inducedG <- NULL
        repressedG <- NULL
        nG <- NULL
        for (i in 1:nModules) {
            nG <- c(nG, sum(data@GeneGrps[, i])) ## number of genes at ith group
            idxTmp <- data@GeneGrps[, i]

            ## induced at array j restricted at group i
            inducedG <- c(inducedG, sum(tableIR[idxTmp,j] == 1))

            ## repr at array j restricted at group i
            repressedG <- c(repressedG, sum(tableIR[idxTmp,j] == -1))

            pValueI <- c(pValueI, dhyper(inducedG[i], inducedA,
            (nGenes-inducedA), nG[i])+phyper(inducedG[i], inducedA,
            (nGenes-inducedA), nG[i], lower.tail=FALSE))

            pValueR <- c(pValueR, dhyper(repressedG[i], repressedA,
            (nGenes-repressedA), nG[i])+phyper(repressedG[i], repressedA,
            (nGenes-repressedA), nG[i], lower.tail=FALSE))

        }

        ## Adjusting p-values if specified by the user.
        if(adjP != "none") {
            tmp1 <- multtest::mt.rawp2adjp(pValueI, proc=adjP)
            tmp2 <- multtest::mt.rawp2adjp(pValueR, proc=adjP)
            pValueI <- tmp1$adjp[order(tmp1$index), 2]
            pValueR <- tmp2$adjp[order(tmp2$index), 2]
        }
        for (i in 1:nModules) {

            if ((pValueI[i] <= cutPhiper) & (pValueR[i] > cutPhiper) &
            (inducedG[i]/nG[i] > inducedA/nGenes))
                arrayVSmoduleIR[j, i] <- inducedG[i]/nG[i]

            if ((pValueR[i] <= cutPhiper) & (pValueI[i] > cutPhiper) &
            (repressedG[i]/nG[i] > repressedA/nGenes))
                arrayVSmoduleIR[j, i] <- -repressedG[i]/nG[i]

            if ((pValueI[i] <= cutPhiper) & (pValueR[i] <= cutPhiper) &
            (inducedG[i]/nG[i] > inducedA/nGenes) &
            (repressedG[i]/nG[i] > repressedA/nGenes))
                arrayVSmoduleIR[j, i] <- NA

        }
    }
    rownames(arrayVSmoduleIR) <- samples
    colnames(arrayVSmoduleIR) <- modules


    ## Counting the number of arrays by tissue type that present alterations...
    conditions <- unique(samples)
    nConditions <- length(conditions)

    ## Calculating changes by chi-square distribution...
    ##tableFinal <- NULL
    ##for (j in 1:nModules) {
    ##chi <- 0
    ##inducedAll <- sum(arrayVSmoduleIR[, j] > 0)
    ##repressedAll <- sum(arrayVSmoduleIR[, j] < 0)
    ##neutralAll <- sum(arrayVSmoduleIR[, j] == 0)
    ##for (i in 1:nConditions) {
    ##    idxTmp <- is.element(samples, conditions[i])
    ##    nC <- sum(idxTmp)
    ##    inducedC.obs <- sum(arrayVSmoduleIR[idxTmp, j] > 0)
    ##    repressedC.obs <- sum(arrayVSmoduleIR[idxTmp, j] < 0)
    ##    neutralC.obs <- sum(arrayVSmoduleIR[idxTmp, j] == 0)
    ##    inducedC.exp <- (inducedAll/nArrays)*nC
    ##    repressedC.exp <- (repressedAll/nArrays)*nC
    ##    neutralC.exp <- (neutralAll/nArrays)*nC
    ##    if(inducedC.exp != 0 & repressedC.exp != 0 & neutralC.exp != 0)
    ##        chi <- chi+((inducedC.obs-inducedC.exp)^2)/inducedC.exp+
    ##        ((repressedC.obs-repressedC.exp)^2)/repressedC.exp+
    ##        ((neutralC.obs-neutralC.exp)^2)/neutralC.exp
    ##}
    ##print(chi)
    ##Pvalue <- pchisq(chi, df=(2*(nConditions-1)), lower.tail=FALSE)
    ##tableFinal <- c(tableFinal, Pvalue)


    ## Calculating the changes by hipergeometric distribution...
    tableFinal <- matrix(0, nrow=nConditions, ncol=nModules)
    for (i in 1:nConditions) {
        idxTmp <- is.element(samples, conditions[i])
        nC <- sum(idxTmp)
        inducedAll <- NULL
        repressedAll <- NULL
        inducedC <- NULL
        repressedC <- NULL
        pValueI <- NULL
        pValueR <- NULL
        for (j in 1:nModules) {
            inducedAll <- c(inducedAll, sum(arrayVSmoduleIR[, j] > 0,
            na.rm=TRUE))

            repressedAll <- c(repressedAll, sum(arrayVSmoduleIR[, j] < 0,
            na.rm=TRUE))

            inducedC <- c(inducedC, sum(arrayVSmoduleIR[idxTmp, j] > 0,
            na.rm=TRUE))

            repressedC <- c(repressedC, sum(arrayVSmoduleIR[idxTmp, j] < 0,
            na.rm=TRUE))

            pValueI <- c(pValueI, dhyper(inducedC[j], inducedAll[j],
            (nArrays-inducedAll[j]), nC)+phyper(inducedC[j], inducedAll[j],
            (nArrays-inducedAll[j]), nC, lower.tail=FALSE))

            pValueR <- c(pValueR, dhyper(repressedC[j], repressedAll[j],
            (nArrays-repressedAll[j]), nC)+phyper(repressedC[j],
            repressedAll[j], (nArrays-repressedAll[j]), nC, lower.tail=FALSE))

        }

        ## Adjusting p-values if specified by the user.
        if(adjP != "none") {
            tmp1 <- multtest::mt.rawp2adjp(pValueI, proc=adjP)
            tmp2 <- multtest::mt.rawp2adjp(pValueR, proc=adjP)
            pValueI <- tmp1$adjp[order(tmp1$index), 2]
            pValueR <- tmp2$adjp[order(tmp2$index), 2]
        }
        for (j in 1:nModules) {
            if ((pValueI[j] <= cutPhiper) & (pValueR[j] > cutPhiper) &
            (inducedC[j]/nC > inducedAll[j]/nArrays))
                tableFinal[i, j] <- inducedC[j]/nC

            if ((pValueR[j] <= cutPhiper) & (pValueI[j] > cutPhiper) &
            (repressedC[j]/nC > repressedAll[j]/nArrays))
                tableFinal[i, j] <- -repressedC[j]/nC

            if ((pValueI[j] <= cutPhiper) & (pValueR[j] <= cutPhiper) &
            (inducedC[j]/nC > inducedAll[j]/nArrays) &
            (repressedC[j]/nC > repressedAll[j]/nArrays))
                tableFinal[i, j] <- NA

        }
    }
    rownames(tableFinal) <- conditions
    colnames(tableFinal) <- modules



    ##
    ## from here, we calculate the score of the genes...
    ##

    ## Calculating the (global) score of each gene (for all tissue types)...
    score <- list()
    for (j in 1:nModules) {
        scoreTable <- NULL
        arrayI <- unname(arrayVSmoduleIR[, j] > 0)
        arrayR <- unname(arrayVSmoduleIR[, j] < 0)
        for (n in which(data@GeneGrps[, j])) {
            idxTmp <- n
            scoreTmp <- 0
            mu <- 0
            sigma2 <- 0
            for (i in 1:nArrays) {
                if(is.finite(arrayI[i]) & arrayI[i] &
                (mean(tableIR[idxTmp,i])== 1))
                    scoreTmp <- scoreTmp-log(abs(arrayVSmoduleIR[i, j]))

                if(is.finite(arrayR[i]) & arrayR[i] & (mean(tableIR[idxTmp,i])
                == -1))
                    scoreTmp <- scoreTmp-log(abs(arrayVSmoduleIR[i, j]))

                if(is.finite(arrayI[i]) & is.finite(arrayR[i]))
                    if(arrayI[i] | arrayR[i]) {
                        mu <- mu-abs(arrayVSmoduleIR[i, j])*
                        (log(abs(arrayVSmoduleIR[i, j])))

                        sigma2 <- sigma2+abs(arrayVSmoduleIR[i, j])*
                        (1-abs(arrayVSmoduleIR[i, j]))*
                        ((log(abs(arrayVSmoduleIR[i, j])))^2)
                    }
            }
            if(sigma2 != 0)
                pValue <- pnorm(((scoreTmp-mu)/sqrt(sigma2)),
                lower.tail=FALSE)
            else
                pValue <- pnorm(0, lower.tail=FALSE)

            scoreTable <- rbind(scoreTable, c(Score=scoreTmp,
            P.valor=pValue))
        }
        score[[j]] <- scoreTable
        rownames(score[[j]]) <- genes[data@GeneGrps[, j]]
    }
    names(score) <- modules


    ## Calculating the score of each gene (by module) independently for each
    ## tissue type...
    scoreByTissue <- list()
    for (j in 1:nModules) {
        listTmp <- array(0, c(sum(data@GeneGrps[, j]), 2,
        length(unique(samples))), list(genes[data@GeneGrps[, j]],
        c("Score","P.value"), unique(samples)))

        arrayI <- unname(arrayVSmoduleIR[, j] > 0)
        arrayR <- unname(arrayVSmoduleIR[, j] < 0)

        for (k in 1:length(unique(samples))) {
            idxTmp.A <- which(is.element(samples, unique(samples)[k]))
            scoreTable <- NULL
            for (n in which(data@GeneGrps[, j])) {
                idxTmp <- n
                scoreTmp <- 0
                mu <- 0
                sigma2 <- 0
                for (i in idxTmp.A) {
                    if(is.finite(arrayI[i]) & arrayI[i] &
                    (mean(tableIR[idxTmp, i]) == 1))
                        scoreTmp <- scoreTmp-log(abs(arrayVSmoduleIR[i, j]))

                    if(is.finite(arrayR[i]) & arrayR[i] &
                    (mean(tableIR[idxTmp, i]) == -1))
                        scoreTmp <- scoreTmp-log(abs(arrayVSmoduleIR[i, j]))

                    if(is.finite(arrayI[i]) & is.finite(arrayR[i]))
                        if(arrayI[i] | arrayR[i]) {
                            mu <- mu-abs(arrayVSmoduleIR[i,j])*
                            (log(abs(arrayVSmoduleIR[i,j])))

                            sigma2 <- sigma2+abs(arrayVSmoduleIR[i, j])*
                            (1-abs(arrayVSmoduleIR[i, j]))*
                            ((log(abs(arrayVSmoduleIR[i,j])))^2)

                        }
                }
                if(sigma2 != 0)
                    pValue <- pnorm(((scoreTmp-mu)/sqrt(sigma2)),
                    lower.tail=FALSE)
                else
                    pValue <- pnorm(0, lower.tail=FALSE)

                scoreTable <- rbind(scoreTable, c(scoreTmp, pValue))
            }

            listTmp[, , k] <- scoreTable
        }

        scoreByTissue[[j]] <- listTmp
    }

    names(scoreByTissue) <- modules


    ## Picking R and packages version information
    tmp <- sessionInfo()
    vInfo <- list()
    vInfo$R.version <- tmp$R.version$version.string
    vInfo$BasePacks <- tmp$basePkgs
    tmp1 <- NULL
    for (i in 1:length(tmp$otherPkgs))
        tmp1 <- c(tmp1, paste(tmp$otherPkgs[[i]]$Package, "version",
        tmp$otherPkgs[[i]]$Version))

    vInfo$AddPacks <- tmp1


    ## Defining the object to return
    result <- new("maigesActMod", modBySamp=arrayVSmoduleIR,
    modByCond=tableFinal, globalScore=score, Date=date(),
    tissueScore=scoreByTissue, V.info=vInfo)

    return(result)

}
