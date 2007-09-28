\name{deGenesANOVA}
\alias{deGenesANOVA}

\title{
  Function to do differential expression analysis, using ANOVA models
}

\description{
  This function takes an object of class \code{\link{maiges}} and do
  differential expression analysis for the genes onto dataset, comparing
  more than two samples using ANalysis Of VAriance (ANOVA)
  models.
}

\usage{
deGenesANOVA(data=NULL, eBayes=FALSE, retOrig=FALSE,
             retF=FALSE, doClust=TRUE, \dots)
}

\arguments{
  \item{data}{object of class \code{\link{maigesANOVA}}.}
  \item{eBayes}{logical indicating the use or not (default) of empirical
    Bayes statistics implemented in \emph{limma} package.}
  \item{retOrig}{logical indicating if the object of class
    \code{\link[limma:marraylm]{MArrayLM}} from \emph{limma} package must be
    returned. Defaults to FALSE.}
  \item{retF}{logical asking to return the results associated with the
    F test (TRUE) or with the individual contrasts (FALSE - default).}
  \item{doClust}{logical indicating if the object generated from this
    analysis will be used for cluster analysis. Defaults to TRUE.}
  \item{\dots}{additional parameters to function \code{\link[limma]{lmFit}}.}
}

\value{
  The result of this function is an object of class
  \code{\link[limma:marraylm]{MArrayLM}} when \code{retOrig} is TRUE. When it is
  FALSE, the result is an object of class \code{\link{maigesDE}} if
  \code{doClust} if FALSE or of class \code{\link{maigesDEcluster}} if
  it is TRUE.
}

\details{
  The object of class \code{\link{maigesANOVA}} of the argument
  \code{data} is created by the function \code{\link{designANOVA}}. This
  function calculate statistics and p-values using the function
  \code{\link[limma]{lmFit}} from package \emph{limma}.
}

\seealso{
  \code{\link{designANOVA}}, \code{\link[limma]{lmFit}}, \code{\link[limma:marraylm]{MArrayLM}}.
}

\examples{
## Loading the dataset
data(gastro)

## Constructing a maigesANOVA object for the 'Tissue' sample label using
## default model (simple linear model with intercept) and contrasts (all
## parameters are equal between themselves)
gastro.ANOVA = designANOVA(gastro.summ, factors="Tissue")

## Fitting the ANOVA model designed by the above command
gastro.ANOVAfit = deGenesANOVA(gastro.ANOVA, retF=TRUE)
gastro.ANOVAfit
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
