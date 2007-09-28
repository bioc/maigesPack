\name{activeNet}
\alias{activeNet}

\title{Functional classification of gene networks}

\description{
  This function calculate a statistic for each gene network in each
  biological condition that measure the profile of activation of the
  network in that condition. Also the function measures the significance
  of the results.
}

\usage{
activeNet(data=NULL, samples=NULL, sLabelID="Classification",
          type="Rpearson", bRep=1000, alternative = "greater",
          adjP="none")
}

\arguments{
  \item{data}{object of class \code{\link{maiges}} to be used to
    functionally classify gene networks stored in \code{Paths} slot.}
  \item{sLabelID}{character string specifying identification of sample
    label to be used. }
  \item{samples}{a list with character vectors specifying the groups
    that must be compared.}
  \item{type}{character string giving the type of correlation to be
    calculated. May be 'Rpearson' (default), 'pearson', 'kendall',
    'spearman' or 'MI'.}
  \item{bRep}{integer number specifying the bootstraps to be done in the
    correlation test.}
  \item{alternative}{character string specifying the alternative
    hypotheses. May be 'greater' (default) to test the activity of the
    networks in accordance to the to the graph or 'less' to test the
    activity of the network antagonic to the graph.}
  \item{adjP}{character string giving the type of p-value
    adjustment. May be 'Bonferroni', 'Holm', 'Hochberg', 'SidakSS',
    'SidakSD', 'BH', 'BY' or 'none'. Defaults to 'none'. See function
    \code{\link[multtest]{mt.rawp2adjp}} in package multtest for more details.}
}

\value{
  The result of this function is an object of class \code{\link{maigesActNet}}.
}

\details{
  If the argument \code{samples} is NULL, all types defined by the
  sample label given by \code{sLabelID} are used. It is possible to use
  the \code{\link{plot.maigesActNet}} and
  \code{\link{image.maigesActNet}} methods to display the results of this
  analysis.
}

\seealso{
  \code{\link{activeNetScoreHTML}}, \code{\link{maigesActNet}},
  \code{\link{plot.maigesActNet}}, \code{\link{image.maigesActNet}}, \code{\link[multtest]{mt.rawp2adjp}}
}

\examples{
## Loading the dataset
data(gastro)

## Doing functional classification of gene networks for sample Label
## given by 'Tissue'
gastro.net = activeNet(gastro.summ, sLabelID="Tissue")
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
