\name{somMde}
\alias{somMde}

\title{
  Function to do SOM cluster analysis
}

\description{
  This is a function to do SOM (Self Organising Maps) clustering
  analysis for objects of classes \code{\link{maigesDEcluster}}.
}

\usage{
somMde(data, group=c("C", "R")[1], distance="correlation",
       method="complete", sampleT=NULL, doHier=FALSE, sLabelID="SAMPLE",
       gLabelID="GeneName", idxTest=1, adjP="none", nDEgenes=0.05, \dots)

}

\arguments{
  \item{data}{object of class \code{\link{maigesDEcluster}}.}
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
  \item{idxTest}{numerical index of the test to be used to sort the
    genes when clustering objects of class \code{\link{maigesDEcluster}}.}
  \item{adjP}{string specifying the method of p-value adjustment. May be
    'none', 'Bonferroni', 'Holm', 'Hochberg', 'SidakSS', 'SidakSD',
    'BH', 'BY'.}
  \item{nDEgenes}{number of DE genes to be selected. If a real number
    in (0,1) all genes with p.value <= \code{nDEgenes} will be
    used. If an integer, the \code{nDEgenes} genes with smaller
    p-values will be used.}
  \item{\dots}{additional parameters for \code{\link[som]{som}} function.}
}

\details{
  This function implements the SOM clustering method for
  objects resulted from differential expression analysis. The method uses
  the function \code{\link[som]{som}} from package \emph{som}. For
  the adjustment of p-values in the selection of genes differentially
  expressed, we use the function \code{\link[multtest]{mt.rawp2adjp}}
  from package \emph{multtest}.
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

## Doing bootstrap from t statistic test fot 'Type' sample label, k=1000
## specifies one thousand bootstraps
gastro.ttest = deGenes2by2Ttest(gastro.summ, sLabelID="Type")

## SOM cluster with 2 groups adjusting p-values by FDR, and showing all genes
## with p-value < 0.05
somMde(gastro.ttest, sLabelID="Type", adjP="BH", nDEgenes=0.05,
       xdim=2, ydim=1, topol="rect")

## SOM cluster with 4 groups adjusting p-values by FDR, and showing all genes
## with p-value < 0.05
somMde(gastro.ttest, sLabelID="Type", adjP="BH", nDEgenes=0.05,
       xdim=2, ydim=2, topol="rect")
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{hplot}
