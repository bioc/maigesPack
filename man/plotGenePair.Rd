\name{plotGenePair}
\alias{plotGenePair}

\title{
  Scatter plots for pair of genes
}

\description{
  This function displays scatter plots for pair of genes that presented
  altered correlation values in Relevance Network analysis.
}

\usage{
plotGenePair(obj, gene1, gene2, posL=NULL, rCor=TRUE)
}

\arguments{
  \item{obj}{object of class \code{\link{maigesRelNetM}}.}
  \item{gene1}{character string giving the first gene identification.}
  \item{gene2}{character string giving the first gene identification.}
  \item{posL}{numerical vector of length 2, specifying the x and y
    position of the legend.}
  \item{rCor}{logical specifying if the correlation are robust
    (calculated by the function \code{\link{robustCorr}}. Defaults to TRUE.}
}

\details{
  This function only picks the result of the \code{\link{relNetworkM}}
  and display scatter plots for a pair of genes giving the regression
  lines and the correlation values for the two biological groups tested.
}

\value{
  This function don't return any object.
}

\seealso{
  \code{\link{maigesRelNetM}}, \code{\link{robustCorr}}, \code{\link{relNetworkM}}.
}

\examples{
## Loading the dataset
data(gastro)

## Constructing the relevance network for sample
## 'Tissue' comparing 'Neso' and 'Aeso' for the 1st gene group
gastro.net = relNetworkM(gastro.summ, sLabelID="Tissue", 
  samples = list(Neso="Neso", Aeso="Aeso"), geneGrp=11,
  type="Rpearson")

## As the sample size is small, because we used a small fraction of the
## genes from the original dataset, this isn't so reliable.
plotGenePair(gastro.net, "KLK13", "EVPL")
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
