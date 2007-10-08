\name{image}
\alias{image.maigesRaw}
\alias{image.maiges}
\alias{image.maigesANOVA}
\alias{image.maigesRelNetB}
\alias{image.maigesRelNetM}
\alias{image.maigesActMod}
\alias{image.maigesActNet}
\alias{image}

\title{
  Method image for objects defined in this package
}

\description{
  Generic function \code{\link[graphics]{image}} to display colour maps
  of numerical values stored in objects defined in this package.
}

\usage{
\method{image}{maigesRaw}(x, \dots)

\method{image}{maiges}(x, \dots)

\method{image}{maigesANOVA}(x, \dots)

\method{image}{maigesRelNetB}(x=NULL, name=NULL, \dots)

\method{image}{maigesRelNetM}(x=NULL, names=NULL, \dots)

\method{image}{maigesActMod}(x, type=c("S","C")[2], keepEmpty=FALSE, ...)

\method{image}{maigesActNet}(x, type=c("score","p-value")[1], ...)
}

\arguments{
  \item{x}{an object of class \code{\link{maigesRaw}},
    \code{\link{maiges}}, \code{\link{maigesANOVA}},
    \code{\link{maigesRelNetB}}, \code{\link{maigesRelNetM}},
    \code{\link{maigesActMod}} or \code{\link{maigesActNet}} defined in
    this package.}
  \item{name}{character string giving a name for sample type tested to
    be plotted as a name in the method for class \code{\link{maigesRelNetB}}.}
  \item{names}{similar to the previous one, but it is a vector of length
    3 for class \code{\link{maigesRelNetM}}.}
  \item{type}{string specifying the type of colour map to be plotted. For
    class \code{\link{maigesActMod}} it must be 'S' or 'C' for samples
    or biological conditions, respectively. For class
    \code{\link{maigesActNet}} it must be 'score' or 'p-value' for the
    statistics or p-values of the tests, respectively.}
  \item{keepEmpty}{logical, if true the results of all gene groups are
    displayed, else only the gene groups that present at least one
    significant result are displayed.}
  \item{\dots}{additional arguments for the generic method
    \code{\link[graphics]{image}} from \emph{graphics} package or
    \code{\link[marray]{maImage}} defined in package \emph{marray} (for
    \code{maigesRaw}, \code{maiges} or \code{maigesANOVA} classes), in
    this case the additional parameters must not be named, because these
    names conflict with the \code{boxplot} generic function definition.}
}

\details{
  This method uses the function \code{\link[marray]{maImage}} from
  marray package to display colour maps of accessor methods defined into
  \emph{marray} package for objects of class \code{\link{maiges}} or
  \code{\link{maigesRaw}} and \code{\link{maigesANOVA}}.

  For objects of class
  \code{\link{maigesRelNetM}} the method displays 3 colour maps
  representing the correlation values for the two groups tested and the
  p-values of the tests.

  For class \code{\link{maigesRelNetB}} it
  displays the correlation values for the type tested.

  In objects of
  class \code{\link{maigesActMod}} it displays the fraction of genes
  induced or repressed for each gene group, by samples or biological
  type.

  Finally, for class \code{\link{maigesActNet}}, the method
  display the matrix of statistics or p-values of the tests.

  Pay attention, if you specify the parameter \emph{x} (but not named)
  for \code{\link[marray]{maImage}}) it will plot the M values instead
  of W (default). 
}

\examples{
## Loading the dataset
data(gastro)

## Doing image plots (using package marray), default method showing the
## W values (for 1st chip), after showing the A values (2nd chip) and
## red background (20th chip). As we are using as small piece of a large
## dataset, this command gives some inconvenient warnings...
\dontrun{
image(gastro.raw[,1])
image(gastro.raw[,2], "maA")
image(gastro.raw[,20], "maRb")

## Example for normalized objects (showing A values for the 5th chip).
image(gastro.norm[,5], "maA")
}



## Example for object of class maigesRelNetB

## Constructing the relevance network (Butte's method) for sample
## 'Tissue' equal to 'Neso' for the 1st gene group
gastro.net = relNetworkB(gastro.summ, sLabelID="Tissue", 
  samples="Neso", geneGrp=1, type="Rpearson")

image(gastro.net)




## Example for object of class maigesRelNetM

## Constructing the relevance network for sample
## 'Tissue' comparing 'Neso' and 'Aeso' for the 1st gene group
gastro.net = relNetworkM(gastro.summ, sLabelID="Tissue", 
  samples = list(Neso="Neso", Aeso="Aeso"), geneGrp=11,
  type="Rpearson")

image(gastro.net)
}

\author{
  Gustavo H. Esteves <gesteves@vision.ime.usp.br>
}

\seealso{
  \code{\link[graphics]{image}} in the \emph{graphics} package and
  \code{\link[marray]{maImage}} in \emph{marray} package.
}

\keyword{array}
