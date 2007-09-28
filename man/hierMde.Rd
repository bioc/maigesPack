\name{hierMde}
\alias{hierMde}

\title{
  Function to do hierarchical cluster analysis
}

\description{
  This is a function to do hierarchical clustering
  analysis for objects of classe \code{\link{maigesDEcluster}}.
}

\usage{
hierMde(data, group=c("C", "R", "B")[1], distance="correlation",
        method="complete", doHeat=TRUE, sLabelID="SAMPLE",
        gLabelID="GeneName", idxTest=1, adjP="BH",
        nDEgenes=0.05, \dots)

}

\arguments{
  \item{data}{object of class \code{\link{maigesDEcluster}}.}
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
  \item{idxTest}{numerical index of the test to be used to sort the
    genes when clustering objects of class \code{\link{maigesDEcluster}}.}
  \item{adjP}{string specifying the method of p-value adjustment. May be
    'none', 'Bonferroni', 'Holm', 'Hochberg', 'SidakSS', 'SidakSD',
    'BH', 'BY'.}
  \item{nDEgenes}{number of DE genes to be selected. If a real number
    in (0,1) all genes with p.value <= \code{nDEgenes} will be
    used. If an integer, the \code{nDEgenes} genes with smaller
    p-values will be used.}
  \item{\dots}{additional parameters for \code{\link[stats]{heatmap}} function.}
}

\details{
  This function implements the hierarchical clustering method for
  objects resulted from differential expression analysis. The default
  function for hierarchical clustering is the
  \code{\link[stats]{hclust}}. For the adjustment of p-values in the
  selection of genes differentially expressed, we use the function
  \code{\link[multtest]{mt.rawp2adjp}} from package \emph{multtest}.
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

## Doing bootstrap from t statistic test fot 'Type' sample label, k=1000
## specifies one thousand bootstraps
gastro.ttest = deGenes2by2Ttest(gastro.summ, sLabelID="Type")

## Hierarchical cluster adjusting p-values by FDR, and showing all genes
## with p-value < 0.05
hierMde(gastro.ttest, sLabelID="Type", adjP="BH", nDEgenes=0.05)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{hplot}
