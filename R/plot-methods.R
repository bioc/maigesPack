## Define methods for plot generic function
##
## Gustavo H. Esteves
## 27/05/07
##
## Version: 1.1
##


## For maigesRaw class
plot.maigesRaw <- function(x, bkgSub="subtract", z=NULL, legend.func=NULL,
ylab="W", ...) {

    ## Correcting background
    tmp <- backgroundCorrect(as(x, "RGList"), bkgSub)
    tmp <- as(tmp, "marrayRaw")
    tmp <- as(tmp, "marrayNorm")

    ## indexing ref labelled with green
    idx <- tolower(getLabels(x, "Ref")) == "red"
    if(sum(idx) > 0)
        maM(tmp)[, idx] <- -maM(tmp)[, idx]

    maLabels(maTargets(tmp)) <- as.character(maInfo(maTargets(tmp))[[1]])
    maPlot(tmp, z=z, legend.func=legend.func, ylab=ylab, ...)
}


## For maiges class
plot.maiges <- function(x, z=NULL, legend.func=NULL, ylab="W", ...) {

    ## Converting to marrayNorm class (from marray package)
    tmp <- as(x, "marrayNorm")

    ## indexing ref labelled with green
    idx <- tolower(getLabels(x, "Ref")) == "red"
    if(sum(idx) > 0)
        maM(tmp)[, idx] <- -maM(tmp)[, idx]

    maPlot(tmp, z=z, legend.func=legend.func, ylab=ylab, ...)

}


## For maigesANOVA class
plot.maigesANOVA <- plot.maiges


## For maigesDE class (shows volcano plot)
plot.maigesDE <- function(x, adjP="none", idx=1, ...) {

    ## Adjusting p-values if specified by the user.
    if(adjP != "none") {
        tmp1 <- multtest::mt.rawp2adjp(x@p.value[, idx], proc=adjP)
        x@p.value[, idx] <- tmp1$adjp[order(tmp1$index), 2]
    }

    if(sum(dim(x@fold)) == 0)
        stop("I can't plot volcano with fold slot empty.")
    else
        tmp1 <- x@fold[, idx]

    tmp2 <- -log10(x@p.value[, idx])
    if(!is.null(colnames(x@stat)[idx]))
        plotName <- colnames(x@stat)[idx]
    else
        plotName <- x@test

    plot(tmp1, tmp2, main=paste("Volcano plot -", plotName), xlab="log2(fold)",
    ylab=paste("-log10(P.value)"), ...)

}


## For maigesDEcluster class
plot.maigesDEcluster <- function(x, adjP="none", idx=1, ...) {

    ## Adjusting p-values if specified by the user.
    if(adjP != "none") {
        tmp1 <- multtest::mt.rawp2adjp(x@p.value[, idx], proc=adjP)
        x@p.value[, idx] <- tmp1$adjp[order(tmp1$index), 2]
    }

    if(sum(dim(x@fold)) == 0)
        tmp1 <- apply(x@W, 1, mean, na.rm=TRUE)
    else
        tmp1 <- x@fold[, idx]

    tmp2 <- -log10(x@p.value[, idx])
    if(!is.null(colnames(x@stat)[idx]))
        plotName <- colnames(x@stat)[idx]
    else
        plotName <- x@test

    plot(tmp1, tmp2, main=paste("Volcano plot -", plotName), xlab="log2(fold)",
    ylab=paste("-log10(P.value)"), ...)

}


## For maigesClass class (bi or tri-dimensional graphs of classifiers)
plot.maigesClass <- function(x, idx=1, ...) {

    ## Getting sample types
    types <- colnames(x@W)
    uTypes <- unique(types)

    if(dim(x@cliques)[2] == 2) {
        tmp1 <- x@W[x@cliques.idx[idx, 1], ]
        tmp2 <- x@W[x@cliques.idx[idx, 2], ]
        xlimite <- c(2*min(tmp1)-quantile(tmp1, 0.35), max(tmp1))
        ylimite <- c(min(tmp2), 2*max(tmp2)-quantile(tmp2, 0.65))

        toLegend <- c(2*min(tmp1)-quantile(tmp1, 0.3),
        2*max(tmp2)-quantile(tmp2, 0.7))

        plot(tmp1[types == uTypes[1]], tmp2[types == uTypes[1]], ylim=ylimite,
        xlim=xlimite, main="Classification graphic", xlab=x@cliques[idx,1],
        ylab=x@cliques[idx,2], pch=22, bg="red", col="red", ...)

        points(tmp1[types == uTypes[2]], tmp2[types == uTypes[2]], pch=19,
        col="green")

        legend(toLegend[1], toLegend[2], uTypes, pch=c(22,19),
        pt.bg=c("red","green"), col=c("red", "green"))

    }
    else if(dim(x@cliques)[2] == 3) {

        ## Loading required library
        require("rgl")

        tmp1 <- x@W[x@cliques.idx[idx, 1], ]
        tmp2 <- x@W[x@cliques.idx[idx, 2], ]
        tmp3 <- x@W[x@cliques.idx[idx, 3], ]

        labels <- c(x@cliques[idx, 1], x@cliques[idx, 2], x@cliques[idx, 3])

        rgl::open3d()

        rgl::text3d(tmp1[types == uTypes[1]], tmp2[types == uTypes[1]],
        tmp3[types == uTypes[1]], text=types[types == uTypes[1]], col="red")

        rgl::text3d(tmp1[types != uTypes[1]], tmp2[types != uTypes[1]],
        tmp3[types != uTypes[1]], text=types[types != uTypes[1]], col="green")

        rgl::decorate3d(xlab=labels[1], ylab=labels[2], zlab=labels[3],
        box=FALSE)

        rgl::bbox3d(color=c("#767676","black"), emission="#767676",
        specular="#767676", shininess=10, alpha=0.9, marklen=50)
    }
    else
        stop("I can't plot graphics with more than 3 genes in the classifiers!")

}


## For maigesRelNetM class (circular graph with significative iterations)
plot.maigesRelNetM <- function(x=NULL, cutPval=0.05, names=NULL, ...) {

    ## Giving default names if it is NULL
    if(is.null(names))
        names <- c(x@types, "Significance of differences")


    ## Defining additional functions
    graphPath <- function(data=NULL, cuttoffPvalue) {
        N <- length(rownames(data@Corr1))
        corObj1 <- data@Corr1
        corObj2 <- data@Corr2
        Diff <- data@DifP
        ## Construct graphs centered on the 1st table
        vertices <- rownames(corObj1)
        arestas1 <- vector("list", length=length(vertices))
        names(arestas1) <- vertices
        arestas3 <- arestas2 <- arestas1
        for (i in 1:length(vertices)) {
            idx <- Diff[i, ] <= cuttoffPvalue
            arestas1[[i]] <- list(edges=vertices[idx],
            weights=unname(corObj1[i, idx]))

            arestas2[[i]] <- list(edges=vertices[idx],
            weights=unname(corObj2[i, idx]))

            arestas3[[i]] <- list(edges=vertices[idx],
            weights=unname(Diff[i, idx]))
        }

        ## Defining the object with the graphs
        Path <- list(Type1=new("graphNEL", vertices, edgeL=arestas1),
        Type2=new("graphNEL", vertices, edgeL=arestas2),
        Dif=new("graphNEL", vertices, edgeL=arestas3))

        return(Path)
    }

    plot.dots <- function(xy, v, n) {
        text(xy[, 1], xy[, 2], v, cex=1.2, col="blue")
    }

    draw.edges <- function(coor, vertices, edges, alpha, sorting) {

        main <- edges[[1]]
        a <- coor[which(vertices == main), ]
        if(sorting == "A")
            strength <- order(abs(edges[[3]]))
        else
            strength <- order(abs(edges[[3]]), decreasing=TRUE)

        for (k in 1:length(edges[[2]])) {
            b <- coor[edges[[2]][k], ]
            ba <- b-a
            ba <- ba/sqrt(sum(ba*ba))
            x <- a+ba*alpha
            y <- b-ba*alpha
            if (edges[[3]][k] > 0 )
                color <- "red"
            else
                color <- "green"
            a1 <- c(x[1], y[1])
            a2 <- c(x[2], y[2])

            lines(a1, a2, lwd=strength[k], lty=1, col=color)
        }
    }

    def.coor <- function(ce, k, h, w) {
        if (k == 1)
            return(ce)
        else if (k == 2) {
            r1 <- c(ce[1], ce[1])
            r2 <- c(ce[2]+h*0.3, ce[2]-h*0.3)
        }
        else if (k == 3) {
            r1 <- c(ce[1], ce[1], ce[1])
            r2 <- c(ce[2]+h*0.25, ce[2], ce[2]-h*0.25)
        }
        else if (k == 4) {
            r1 <- c(ce[1]-w*0.3, ce[1]+w*0.3, ce[1]+w*0.3, ce[1]-w*0.3)
            r2 <- c(ce[2]-h*0.3, ce[2]-h*0.3, ce[2]+h*0.3, ce[2]+h*0.3)
        }
        else {
            a <- 1
            z <- seq(a, a+2*pi, len=k+1)
            z <- z[-1]
            r1 <- ce[1]+w/2.5*cos(z)
            r2 <- ce[2]+h/2.5*sin(z)
        }
        cbind(r1, r2)
    }

    plotGraph <- function (graph, coor=NULL, alpha=2.5, main="Some Graph",
    sorting) {

        ## Define the nodes and the length of the graph (number of nodes)
        v <- nodes(graph)
        n <- length(v)

        ## Create an empty draw device
        plot(c(0, 100), c(0, 100), type="n", axes=FALSE, xlab="", ylab="",
        main=main)

        ## Define an initial center to the graph
        center <- matrix(c(50, 50), ncol=2)

        ## Create an object of the locations for each node
        if (is.null(coor))
            coor <- def.coor(center, n, 100, 100)

        ## Plot each node on the device
        plot.dots(coor, v, n)

        ## Locate the nodes to plot on the figure
        for (i in 1:n) {
            if(length(edgeL(graph)[[i]]$edges) > 0) {
                elo <- list(v[i], edgeL(graph)[[i]]$edges,
                as.numeric(edgeWeights(graph, i)[[1]]))
                draw.edges(coor, v, elo, alpha, sorting)
            }
        }
        colnames(coor) <- c("x", "y")
        return(invisible(coor))
    }


    ## Ploting the graphs
    graph <- graphPath(x, cutPval)

    par(mfrow=c(1,3))
    tmp <- plotGraph(graph[[1]], main=names[1], sorting="A")
    plotGraph(graph[[2]], coor=tmp, main=names[2], sorting="A")
    plotGraph(graph[[3]], coor=tmp, main=names[3], sorting="D")

}


## For maigesRelNetB class (circular graph with significative iterations)
plot.maigesRelNetB <- function(x=NULL, cutPval=0.05, cutCor=NULL, name=NULL,
...) {

    ## Some tests
    if(!is.null(cutPval) & !is.null(cutCor))
        stop("You must specify only, cutPval or cutCor.")

    if(is.null(name))
        name <- x@type

    ## Defining additional functions
    graphPath <- function(data=NULL, cuttoffCor=NULL, cuttoffP=NULL) {

        N <- length(rownames(data@Corr))
        corObj <- data@Corr
        ## Construct graphs centered on the 1st table
        vertices <- rownames(corObj)
        arestas <- vector("list", length=length(vertices))
        names(arestas) <- vertices
        for (i in 1:length(vertices)) {
            if(!is.null(cuttoffCor))
                idx <- abs(corObj[i, ]) >= cuttoffCor
            if(!is.null(cuttoffP))
                idx <- abs(data@Pval[i, ]) <= cuttoffP
            arestas[[i]] <- list(edges=vertices[idx], weights=unname(corObj[i, idx]))
        }

        ## Defining the object with the graphs
        Path <- new("graphNEL", vertices, edgeL=arestas)
        return(Path)
    }

    plot.dots <- function(xy, v, n) {
        text(xy[, 1], xy[, 2], v, cex=1.2, col="blue")
    }

    draw.edges <- function(coor, vertices, edges, alpha, sorting) {

        main <- edges[[1]]
        a <- coor[which(vertices == main), ]
        if(sorting == "A")
            strength <- order(abs(edges[[3]]))
        else
            strength <- order(abs(edges[[3]]), decreasing=TRUE)

        for (k in 1:length(edges[[2]])) {
            b <- coor[edges[[2]][k], ]
            ba <- b-a
            ba <- ba/sqrt(sum(ba*ba))
            x <- a+ba*alpha
            y <- b-ba*alpha
            if (edges[[3]][k] > 0 )
                color <- "red"
            else
                color <- "green"
            a1 <- c(x[1], y[1])
            a2 <- c(x[2], y[2])
            lines(a1, a2, lwd=strength[k], lty=1, col=color)
        }
    }

    def.coor <- function(ce, k, h, w) {
        if (k == 1)
            return(ce)
        else if (k == 2) {
            r1 <- c(ce[1], ce[1])
            r2 <- c(ce[2]+h*0.3, ce[2]-h*0.3)
        }
        else if (k == 3) {
            r1 <- c(ce[1], ce[1], ce[1])
            r2 <- c(ce[2]+h*0.25, ce[2], ce[2]-h*0.25)
        }
        else if (k == 4) {
            r1 <- c(ce[1]-w*0.3, ce[1]+w*0.3, ce[1]+w*0.3, ce[1]-w*0.3)
            r2 <- c(ce[2]-h*0.3, ce[2]-h*0.3, ce[2]+h*0.3, ce[2]+h*0.3)
        }
        else {
            a <- 1
            z <- seq(a, a+2*pi, len=k+1)
            z <- z[-1]
            r1 <- ce[1]+w/2.5*cos(z)
            r2 <- ce[2]+h/2.5*sin(z)
        }
        cbind(r1, r2)
    }

    plotGraph <- function (graph, coor=NULL, alpha=2.5, main="Some Graph",
    sorting) {

        ## Define the nodes and the length of the graph (number of nodes)
        v <- nodes(graph)
        n <- length(v)

        ## Create an empty draw device
        plot(c(0, 100), c(0, 100), type="n", axes=FALSE, xlab="", ylab="",
        main=main)

        ## Define an initial center to the graph
        center <- matrix(c(50, 50), ncol=2)

        ## Create an object of the locations for each node
        if (is.null(coor)) {
            coor <- def.coor(center, n, 100, 100)
        }
        ## Plot each node on the device
        plot.dots(coor, v, n)
        ## Locate the nodes to plot on the figure
        for (i in 1:n) {
            if(length(edgeL(graph)[[i]]$edges) > 0) {
                elo <- list(v[i], edgeL(graph)[[i]]$edges,
                as.numeric(edgeWeights(graph, i)[[1]]))
                ##if(length(elo[[2]]) > 0)
                draw.edges(coor, v, elo, alpha, sorting)
            }
        }
        colnames(coor) <- c("x", "y")
        return(invisible(coor))
    }


    if(!is.null(cutCor)) {
        ## geting maximum bootstrap value to use in cutCor
        if((is.character(cutCor)) & (cutCor == "max"))
            cutCor <- max(x@maxB[upper.tri(x@maxB)])
        else if((cutCor < 0) | (cutCor >= 1))
            stop("cutCor must a number in [0,1) or 'max'.")

        ## Ploting the graphs
        graph <- graphPath(x, cutCor, NULL)
        plotGraph(graph, main=name, sorting="A")
    }
    else {
        ## Ploting the graphs
        graph <- graphPath(x, NULL, cutPval)
        plotGraph(graph, main=name, sorting="A")
    }
}


## For maigesActMod class (heatmap of the significative results)
plot.maigesActMod <- function(x, type=c("S", "C")[2], keepEmpty=FALSE, ...) {

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
plot.maigesActNet <- function(x, type=c("score","p-value")[1], ...) {

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
