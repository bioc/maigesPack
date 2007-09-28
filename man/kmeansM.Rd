\name{kmeansM}
\alias{kmeansM}

\title{
  Function to do k-means cluster analysis
}

\description{
  This is a function to do k-means clustering
  analysis for objects of classes \code{\link{maiges}},
  \code{\link{maigesRaw}} and \code{\link{maigesANOVA}}. Use the
  function \code{\link{kmeansMde}} for objects of class
  \code{\link{maigesDEcluster}}.
}

\usage{
kmeansM(data, group=c("C", "R")[1], distance="correlation",
        method="complete", sampleT=NULL, doHier=FALSE, sLabelID="SAMPLE",
        gLabelID="GeneName", rmGenes=NULL, rmSamples=NULL, rmBad=TRUE,
        geneGrp=NULL, path=NULL, \dots)
}

\arguments{
  \item{data}{object of class \code{\link{maigesRaw}},
    \code{\link{maiges}} or \code{\link{maigesANOVA}}.}
  \item{group}{character string giving the type of grouping: by rows
    'R' or columns 'C' (default).}
  \item{distance}{char string giving the type of distance to use. Here we
    use the function \code{\link[amap:dist]{Dist}} and the possible values
    are 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary',
    'pearson', 'correlation' (default) and 'spearman'.}
  \item{method}{char string specifying the linkage method for the
    hierarchical cluster. Possible values are 'ward', 'single',
    'complete' (default), 'average', 'mcquitty', 'median' or 'centroid'}
  \item{sampleT}{list with 2 vectors. The first one specify the first
    letter of different sample types to be coloured by distinct colours,
    that are given in the second vector. If NULL (default) no colour is used.}
  \item{doHier}{logical indicating if you want to do the hierarchical
    branch in the opposite dimension of clustering. Defaults to FALSE.}
  \item{sLabelID}{character string specifying the sample label ID to be
    used to label the samples.}
  \item{gLabelID}{character string specifying the gene label ID to be
    used to label the genes.}
  \item{rmGenes}{char list specifying genes to be removed.}
  \item{rmSamples}{char list specifying samples to be removed.}
  \item{rmBad}{logical indicating to remove or not bad spots (slot
    \code{BadSpots} in objects of class \code{\link{maiges}},
    \code{\link{maigesRaw}} or \code{\link{maigesANOVA}}).}
  \item{geneGrp}{numerical or character specifying the gene group to be
    clustered. This is given by the columns of the slot \code{GeneGrps}
    in objects of classes \code{\link{maiges}}, \code{\link{maigesRaw}}
    and \code{\link{maigesANOVA}}.}
  \item{path}{numerical or character specifying the gene network to be
    clustered. This is given by the items of the slot \code{Paths}
    in objects of classes \code{\link{maiges}}, \code{\link{maigesRaw}}
    and \code{\link{maigesANOVA}}.}
  \item{\dots}{additional parameters for \code{\link[amap]{Kmeans}} function.}
}

\details{
  This function implements the k-means clustering method for
  objects of microarray data defined in this package. The method uses
  the function \code{\link[amap]{Kmeans}} from package \emph{amap}.
}

\value{
  This function display the heatmaps and return invisibly a list
  resulted from the function \code{\link[amap]{Kmeans}}.
}

\seealso{
  \code{\link[amap]{Kmeans}} from package
  \emph{amap}. \code{\link{somM}} and \code{\link{hierM}} for displaying SOM and
  hierarchical clusters, respectively.
}

\examples{
## Loading the dataset
data(gastro)

## Doing a K-means cluster with 2 groups using all genes, for maigesRaw class
kmeansM(gastro.raw, rmGenes=c("BLANK","DAP","LYS","PHE", "Q_GENE","THR","TRP"),
        sLabelID="Sample", gLabelID="Name", centers=2)

## The same as above, but for maigesNorm class
kmeansM(gastro.norm, rmGenes=c("BLANK","DAP","LYS","PHE", "Q_GENE","THR","TRP"),
        sLabelID="Sample", gLabelID="Name", centers=2)

## Another example with 3 groups
kmeansM(gastro.norm, rmGenes=c("BLANK","DAP","LYS","PHE", "Q_GENE","THR","TRP"),
        sLabelID="Sample", gLabelID="Name", centers=3)

## If you want to use euclidean distance to group genes (or spots) with
## 4 groups
kmeansM(gastro.summ, rmGenes=c("BLANK","DAP","LYS","PHE", "Q_GENE","THR","TRP"),
        sLabelID="Sample", gLabelID="Name", centers=4, group="R", distance="euclidean")
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{hplot}
