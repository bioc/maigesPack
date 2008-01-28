\name{boxplot}
\alias{boxplot.maigesRaw}
\alias{boxplot.maiges}
\alias{boxplot.maigesANOVA}
\alias{boxplot.maigesDEcluster}
\alias{boxplot}

\title{
  Method boxplot for objects defined in this package
}

\description{
  Generic function \code{\link[graphics]{boxplot}} to display boxplots
  of the data.
}

\usage{
\method{boxplot}{maigesRaw}(x, \dots)

\method{boxplot}{maiges}(x, name=NULL, gLabelID=NULL, sLabelID=NULL, gSamples=NULL, \dots)

\method{boxplot}{maigesANOVA}(x, name=NULL, gLabelID=NULL, sLabelID=NULL, gSamples=NULL, \dots)

\method{boxplot}{maigesDEcluster}(x, name=NULL, gLabelID=NULL, sLabelID=NULL, gSamples=NULL, \dots)
}

\arguments{
  \item{x}{an object of any class defined in this package}
  \item{name}{character string specifying a gene to do boxplot of their
    expression values along the types specified by \code{sLabelID}.}
  \item{gLabelID}{character value giving the name of gene label to be
    used to search for parameter \code{name}.}
  \item{sLabelID}{idem to \code{gLabelID}, specifying the name of sample
    label to be used to separate the gene observations.}
  \item{gSamples}{a named list containing character vectors defining
    groups of samples from \code{sLabelID}.}
  \item{\dots}{additional parameters to \code{\link[graphics]{boxplot}}
    method defined in \emph{graphics} package or to
    \code{\link[marray]{maBoxplot}} method defined in \emph{marray}
    package (for \code{maigesRaw}, \code{maiges} or \code{maigesANOVA}
    classes), in this case the additional parameters must not be named,
    because these names conflict with the \code{boxplot} generic
    function definition.}
}

\details{
  This method uses the function \code{\link[marray]{maBoxplot}} from
  marray package to show boxplots of the W values along all the slides of a
  dataset or along specific accessor methods to stratify the data in
  objects of class \code{\link{maigesRaw}}. For objects of classes
  \code{\link{maiges}} or \code{\link{maigesANOVA}} this is also done if
  the argument \code{name} is NULL, else the method shows boxplots for
  the expression values of the gene specified by \code{name} stratified by
  sample types from \code{sLabelID}. For objects of class
  \code{\link{maigesDEcluster}} only boxplots of genes are produced and
  the argument \code{name} may not be NULL.

  If you especify the \code{y} parameter (but not named), defined for
  the method \code{\link[marray]{maBoxplot}} in package \emph{marray},
  it will be displayed the M values instead of W.
}

\examples{
## Loading the dataset
data(gastro)

## To see the boxplots for W values in all chips
boxplot(gastro.raw)  ## maigesRaw class
boxplot(gastro.norm) ## maigesNorm class
boxplot(gastro.summ) ## summarized data (also maigesNorm class)

## To see the boxplots for W values in individual chips
## separating into print tips.
boxplot(gastro.raw[,1])  ## maigesRaw class, first chip
boxplot(gastro.norm[,8]) ## maigesNorm class, 8th chip
boxplot(gastro.summ[,19]) ## summarized data (also maigesNorm class), 19th chip


## Boxplot for individual genes into ANOVA model fitting
gastro.ANOVA = designANOVA(gastro.summ, factors="Tissue")
gastro.ANOVAfit = deGenesANOVA(gastro.ANOVA, retF=TRUE)

boxplot(gastro.ANOVAfit, name="KLK13", gLabelID="GeneName", 
sLabelID="Tissue")
}

\author{
  Gustavo H. Esteves <gesteves@vision.ime.usp.br>
}

\seealso{
  \code{\link[marray]{maBoxplot}} in the \emph{marray} package and
  \code{\link[graphics]{boxplot}} in the \emph{graphics} package.
}

\keyword{array}
