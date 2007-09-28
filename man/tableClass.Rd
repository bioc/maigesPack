\name{tableClass}
\alias{tableClass}

\title{
  Save HTML or CSV tables of good classifiers (cliques)
}

\description{
  This function takes an object of class \code{\link{maigesClass}}
  resulting from classification analysis and write tables (in HTML or
  CSV format) a table of the cliques, or classifiers.
}

\usage{
tableClass(classComp=NULL, file="./class_result",
           type=c("HTML","CSV")[1], nCliques=NULL)
}

\arguments{
  \item{classComp}{object of class \code{\link{maigesClass}}.}
  \item{file}{character string giving the file name to be saved (without
    extension).}
  \item{type}{character string with the type of file to be saved. Must
    be 'HTML' (default) or 'CSV'.}
  \item{nCliques}{integer specifying the number of cliques to be
    saved. If NULL (default) all cliques in the object are saved.}
}

\value{
  This function don't return any object.
}

\seealso{
  \code{\link{maigesClass}}, \code{\link{classifyLDA}},
  \code{\link{classifySVM}}, \code{\link{classifyKNN}},
  \code{\link{classifyLDAsc}}, \code{\link{classifySVMsc}} and
  \code{\link{classifyKNNsc}}.
}

\examples{
## Loading the dataset
data(gastro)


## Doing LDA classifier with 3 genes for the 6th gene group comparing
## the 2 categories from 'Type' sample label.
gastro.class = classifyLDA(gastro.summ, sLabelID="Type",
  gNameID="GeneName", nGenes=3, geneGrp=6)

tableClass(gastro.class)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
