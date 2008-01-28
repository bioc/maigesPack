\name{kmeansMde}
\alias{kmeansMde}

\title{
  Function to do k-means cluster analysis
}

\description{
  This is a function to do k-means clustering
  analysis for objects of class \code{\link{maigesDEcluster}}.
}

\usage{
kmeansMde(data, group=c("C", "R")[1], distance="correlation",
          method="complete", sampleT=NULL, doHier=FALSE, sLabelID="SAMPLE",
          gLabelID="GeneName", idxTest=1, adjP="none", nDEgenes=0.05, \dots)
}

\arguments{
  \item{data}{object of class \code{\link{maigesDEcluster}}.}
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
  \item{idxTest}{numerical index of the test to be used to sort the
    genes when clustering objects of class \code{\link{maigesDEcluster}}.}
  \item{adjP}{string specifying the method of p-value adjustment. May be
    'none', 'Bonferroni', 'Holm', 'Hochberg', 'SidakSS', 'SidakSD',
    'BH', 'BY'.}
  \item{nDEgenes}{number of DE genes to be selected. If a real number
    in (0,1) all genes with p.value <= \code{nDEgenes} will be
    used. If an integer, the \code{nDEgenes} genes with smaller
    p-values will be used.}
  \item{\dots}{additional parameters for \code{\link[amap]{Kmeans}} function.}
}

\details{
  This function implements the k-means clustering method for
  objects resulted from differential analysis. The method uses
  the function \code{\link[amap]{Kmeans}} from package \emph{amap}. For
  the adjustment of p-values in the selection of genes differentially
  expressed, we use the function \code{\link[multtest]{mt.rawp2adjp}}
  from package \emph{multtest}.
}

\value{
  This function display the heatmaps and return invisibly a list
  resulted from the function \code{\link[amap]{Kmeans}}.
}

\seealso{
  \code{\link[amap]{Kmeans}} from package
  \emph{amap}. \code{\link[multtest]{mt.rawp2adjp}} from package
  \emph{multtest}. \code{\link{somM}} and \code{\link{hierM}} for displaying SOM and
  hierarchical clusters, respectively.
}

\examples{
## Loading the dataset
data(gastro)

## Doing bootstrap from t statistic test fot 'Type' sample label, k=1000
## specifies one thousand bootstraps
gastro.ttest = deGenes2by2Ttest(gastro.summ, sLabelID="Type")

## K-means cluster with 2 groups adjusting p-values by FDR, and showing all genes
## with p-value < 0.05
kmeansMde(gastro.ttest, sLabelID="Type", adjP="BH", nDEgenes=0.05, centers=2)

## K-means cluster with 3 groups adjusting p-values by FDR, and showing all genes
## with p-value < 0.05
kmeansMde(gastro.ttest, sLabelID="Type", adjP="BH", nDEgenes=0.05, centers=3)

dev.off()
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{hplot}
