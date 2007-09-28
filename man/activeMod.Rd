\name{activeMod}
\alias{activeMod}

\title{Functional classification of gene groups}

\description{
  This function calculate profiles of expression for groups of genes in
  each sample of the dataset and for each biological condition (group of
  samples).
}

\usage{
activeMod(data=NULL, gNameID="GeneName", samples=NULL, usePaths=FALSE,
          sLabelID="Classification", adjP="none", cutExp=1,
          cutPhiper=0.05)
}

\arguments{
  \item{data}{object of class \code{\link{maiges}} to be used to
    functionally classify gene groups stored in \code{GeneGrps} slot.}
  \item{gNameID}{character string specifying identification of gene
    label to be used. }
  \item{sLabelID}{idem to the previous argument for identification of
    sample label.}
  \item{samples}{a list with character vectors specifying the groups
    that must be compared.}
  \item{usePaths}{logical specifying if the pathways given in
    \code{Paths} slot must also be used, defaults to FALSE.}
  \item{cutExp}{real number specifying the cutoff for expression levels
    (to discretise the expression)}
  \item{adjP}{character string giving the type of p-value
    adjustment. May be 'Bonferroni', 'Holm', 'Hochberg', 'SidakSS',
    'SidakSD', 'BH', 'BY' or 'none'. Defaults to 'none'. See function
    \code{\link[multtest]{mt.rawp2adjp}} in package multtest for more details.}
  \item{cutPhiper}{p-value cutoff to select significant gene groups.}
}

\value{
  The result of this function is an object of class \code{\link{maigesActMod}}.
}

\details{
  If the argument \code{samples} is NULL, all types defined by the
  sample label given by \code{sLabelID} are used. It is possible to use
  the \code{\link{plot.maigesActMod}} and
  \code{\link{image.maigesActMod}} to display the results of this
  analysis. This function is based in the method proposed by
  Segal et al. (2004).
}

\references{
  Segal, E.; Friedman, N.; Koller, D. and Regev, A. A module map showing
  conditional activity of expression modules in cancer. \bold{Nature
    Genetics}, 36, 1090-1098, 2004. (\url{http://www.nature.com/ng/journal/v36/n10/abs/ng1434.html})
}

\seealso{
  \code{\link{activeModScoreHTML}}, \code{\link{maigesActMod}},
  \code{\link{plot.maigesActMod}}, \code{\link{image.maigesActMod}}, \code{\link[multtest]{mt.rawp2adjp}}
}

\examples{
## Loading a little dataset
data(gastro)

## Doing functional classification of gene groups for 'Tissue' sample label
gastro.mod = activeMod(gastro.summ, sLabelID="Tissue", cutExp=1,
  cutPhiper=0.05)

## Doing functional classification of gene groups together with the
## networks given by Paths slot for 'Tissue' sample label. Also we are
## using a cuttoff for p-value of hipergeometric test as 0.1
gastro.mod = activeMod(gastro.summ, sLabelID="Tissue", cutExp=1,
  cutPhiper=0.1, usePaths=TRUE)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
