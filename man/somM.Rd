\name{somM}
\alias{somM}

\title{
  Function to do SOM cluster analysis
}

\description{
  This is a function to do SOM (Self Organising Maps) clustering
  analysis for objects of classes \code{\link{maiges}},
  \code{\link{maigesRaw}} and \code{\link{maigesANOVA}}. Use the
  function \code{\link{somMde}} for objects of class
  \code{\link{maigesDEcluster}}.
}

\usage{
somM(data, group=c("C", "R")[1], distance="correlation",
     method="complete", sampleT=NULL, doHier=FALSE, sLabelID="SAMPLE",
     gLabelID="GeneName", rmGenes=NULL, rmSamples=NULL, rmBad=TRUE,
     geneGrp=NULL, path=NULL, \dots)
}

\arguments{
  \item{data}{object of class \code{\link{maigesRaw}}, \code{\link{maiges}},
    or \code{\link{maigesANOVA}}.}
  \item{group}{character string giving the type of grouping: by rows
    'R' or columns 'C' (default).}
  \item{distance}{char string giving the type of distance to use. Only two
    options are available here: 'euclidean' and 'correlation' (default).}
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
  \item{\dots}{additional parameters for \code{\link[som]{som}} function.}
}

\details{
  This function implements the SOM clustering method for
  objects of microarray data defined in this package. The method uses
  the function \code{\link[som]{som}} from package \emph{som}.
}

\value{
  This function display the heatmaps and return invisibly an object
  of class \code{som} resulted from the function \code{\link[som]{som}}.
}

\seealso{
  \code{\link[som]{som}} from package \emph{som}.
  \code{\link{kmeansM}} and \code{\link{hierM}} for displaying k-means and
  hierarchical clusters, respectively.
}

\examples{
## Loading the dataset
data(gastro)

## Doing a SOM cluster with 2 groups using all genes, for maigesRaw class
somM(gastro.raw, rmGenes=c("BLANK","DAP","LYS","PHE", "Q_GENE","THR","TRP"),
     sLabelID="Sample", gLabelID="Name", xdim=2, ydim=1, topol="rect")

## Doing a SOM cluster with 3 groups using all genes, for maigesNorm class
somM(gastro.norm, rmGenes=c("BLANK","DAP","LYS","PHE", "Q_GENE","THR","TRP"),
     sLabelID="Sample", gLabelID="Name", xdim=3, ydim=1, topol="rect")

## Another example with 4 groups
somM(gastro.norm, rmGenes=c("BLANK","DAP","LYS","PHE", "Q_GENE","THR","TRP"),
     sLabelID="Sample", gLabelID="Name", xdim=2, ydim=2, topol="rect")

## If you want to use euclidean distance to group genes (or spots), with
## 3 groups
somM(gastro.summ, rmGenes=c("BLANK","DAP","LYS","PHE", "Q_GENE","THR","TRP"),
     sLabelID="Sample", gLabelID="Name", group="R",
     distance="euclidean", xdim=3, ydim=1, topol="rect")
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{hplot}
