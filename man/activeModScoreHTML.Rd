\name{activeModScoreHTML}
\alias{activeModScoreHTML}

\title{
  Save HTML file with global gene scores from functional gene
  groups classification
}

\description{
  This function takes an object of class \code{\link{maigesActMod}},
  that is generated using the function \code{\link{activeMod}} to do
  functional classification of gene groups, and save an HTML file with
  global score for genes separated by gene groups (modules).
}

\usage{
activeModScoreHTML(mod=NULL, dir="./", fileSave="scores")
}

\arguments{
  \item{mod}{object of class \code{\link{maigesActMod}}, resulted from
    functional classification of gene groups.}
  \item{dir}{character string giving the folder to save the file.}
  \item{fileSave}{string giving the file name. You don't need to put
    the extension 'html', it will be put automatically.}
}

\value{
  This function generates an HTML file and don't return any R value or object.
}

\seealso{
  \code{\link{plot.maigesActMod}}, \code{\link{activeMod}}
}

\examples{
## Loading a little dataset
data(gastro)

## Doing functional classification of gene groups for 'Tissue' sample label
gastro.mod = activeMod(gastro.summ, sLabelID="Tissue", cutExp=1,
  cutPhiper=0.05)

activeModScoreHTML(gastro.mod)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
