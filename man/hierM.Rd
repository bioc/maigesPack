\name{hierM}
\alias{hierM}

\title{
  Function to do hierarchical cluster analysis
}

\description{
  This is a function to do hierarchical clustering
  analysis for objects of classes \code{\link{maiges}},
  \code{\link{maigesRaw}} and \code{\link{maigesANOVA}}. Use the
  function \code{\link{hierMde}} for objects of class
  \code{\link{maigesDEcluster}}.
}

\usage{
hierM(data, group=c("C", "R", "B")[1], distance="correlation",
      method="complete", doHeat=TRUE, sLabelID="SAMPLE",
      gLabelID="GeneName", rmGenes=NULL, rmSamples=NULL,
      rmBad=TRUE, geneGrp=NULL, path=NULL, \dots)
}

\arguments{
  \item{data}{object of class \code{\link{maigesRaw}}, \code{\link{maiges}},
    \code{\link{maigesANOVA}} or \code{\link{maigesDEcluster}}.}
  \item{group}{character string giving the type of grouping: by rows
    'R', columns 'C' (default) or both 'B'.}
  \item{distance}{char string giving the type of distance to use. Here we
    use the function \code{\link[amap:dist]{Dist}} and the possible values
    are 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary',
    'pearson', 'correlation' (default) and 'spearman'.}
  \item{method}{char string specifying the linkage method for the
    hierarchical cluster. Possible values are 'ward', 'single',
    'complete' (default), 'average', 'mcquitty', 'median' or 'centroid'}
  \item{doHeat}{logical indicating to do or not the heatmap. If FALSE,
    only the dendrogram is displayed.}
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
  \item{\dots}{additional parameters for \code{\link[stats]{heatmap}} function.}
}

\details{
  This function implements the hierarchical clustering method for
  objects of microarray data defined in this package. The default
  function for hierarchical clustering is the
  \code{\link[stats]{hclust}}.
}

\value{
  This function display the heatmaps and don't return any object or value.
}

\seealso{
  \code{\link{somM}} and \code{\link{kmeansM}} for displaying SOM and
  k-means clusters, respectively.
}

\examples{
## Loading the dataset
data(gastro)

## Doing a hierarchical cluster using all genes, for maigesRaw class
hierM(gastro.raw, rmGenes=c("BLANK","DAP","LYS","PHE", "Q_GENE","THR","TRP"),
      sLabelID="Sample", gLabelID="Name", doHeat=FALSE)

## Doing a hierarchical cluster using all genes, for maigesNorm class
hierM(gastro.norm, rmGenes=c("BLANK","DAP","LYS","PHE", "Q_GENE","THR","TRP"),
      sLabelID="Sample", gLabelID="Name", doHeat=FALSE)

## If you want to show the heatmap do
hierM(gastro.norm, rmGenes=c("BLANK","DAP","LYS","PHE", "Q_GENE","THR","TRP"),
      sLabelID="Sample", gLabelID="Name", doHeat=TRUE)

## If you want to show the hierarchical branch in both margins do
hierM(gastro.summ, rmGenes=c("BLANK","DAP","LYS","PHE", "Q_GENE","THR","TRP"),
      sLabelID="Sample", gLabelID="Name", doHeat=TRUE, group="B")

## If you want to use euclidean distance only into rows (spots or genes)
hierM(gastro.summ, rmGenes=c("BLANK","DAP","LYS","PHE", "Q_GENE","THR","TRP"),
      sLabelID="Sample", gLabelID="Name", doHeat=FALSE, group="R", distance="euclidean")
}


\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{hplot}
