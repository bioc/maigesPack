\name{activeNetScoreHTML}
\alias{activeNetScoreHTML}

\title{
  Save HTML file with scores and p-values from functional gene
  networks classification
}

\description{
  This function takes an object of class \code{\link{maigesActNet}},
  that is generated using the function \code{\link{activeNet}} to do
  functional classification of gene groups, and save an HTML file with
  global score for genes separated by gene groups (modules).
}

\usage{
activeNetScoreHTML(mod=NULL, dir="./", fileSave="scores")
}

\arguments{
  \item{mod}{object of class \code{\link{maigesActNet}}, resulted from
    functional classification of gene networks.}
  \item{dir}{character string giving the folder to save the file.}
  \item{fileSave}{string giving the file name. You don't need to put
    the extension 'html', it will be put automatically.}
}

\value{
  This function generates an HTML file and don't return any R value or object.
}

\seealso{
  \code{\link{plot.maigesActNet}}, \code{\link{activeNet}}
}

\examples{
## Loading the dataset
data(gastro)

## Doing functional classification of gene networks for sample Label
## given by 'Tissue'
gastro.net = activeNet(gastro.summ, sLabelID="Tissue")

activeNetScoreHTML(gastro.net)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
