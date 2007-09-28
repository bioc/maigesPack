\name{getLabels}
\alias{getLabels.default}
\alias{getLabels.maigesDE}
\alias{getLabels.maigesDEcluster}
\alias{getLabels.RGList}
\alias{getLabels.MAList}
\alias{getLabels.marrayRaw}
\alias{getLabels.marrayNorm}
\alias{getLabels}

\title{
  Method getLabels to pick gene and sample labels
}

\description{
  Generic function \code{\link{getLabels}} to extract labels given an
  ID to genes or samples.
}

\usage{
getLabels(obj, labelID=NULL, sLabel=TRUE)

\method{getLabels}{default}(obj=NULL, labelID=NULL, sLabel=TRUE)

\method{getLabels}{maigesDE}(obj=NULL, labelID=NULL, sLabel=TRUE)

\method{getLabels}{maigesDEcluster}(obj=NULL, labelID=NULL, sLabel=TRUE)

\method{getLabels}{RGList}(obj=NULL, labelID=NULL, sLabel=TRUE)

\method{getLabels}{MAList}(obj=NULL, labelID=NULL, sLabel=TRUE)

\method{getLabels}{marrayRaw}(obj=NULL, labelID=NULL, sLabel=TRUE)

\method{getLabels}{marrayNorm}(obj=NULL, labelID=NULL, sLabel=TRUE)

}

\arguments{
  \item{obj}{object to look for labels. Methods defined for classes
    \code{\link{maigesRaw}}, \code{\link{maiges}},
    \code{\link{maigesDE}}, \code{\link{maigesDEcluster}}, \code{\link[limma:rglist]{RGList}},
    \code{\link[limma:malist]{MAList}}, \code{\link[marray:marrayRaw-class]{marrayRaw}} and
    \code{\link[marray:marrayNorm-class]{marrayNorm}}.}
  \item{labelID}{character string with label name to be searched.}
  \item{sLabel}{logical indicating search in the sample labels, defaults
    to TRUE. If FALSE search is done for gene labels.}
}

\details{
  The name of gene labels are done by the names of the slot \code{Glabels}
  in objects of classes \code{\link{maigesRaw}} or
  \code{\link{maiges}}, the slot \code{GeneInfo} in objects of classes
  \code{\link{maigesDE}} or \code{\link{maigesDEcluster}}, the slot
  \code{genes} in objects of classes \code{\link[limma:rglist]{RGList}} or
  \code{\link[limma:malist]{MAList}} and the slot \code{maGnames@maInfo} in
  objects of classes \code{\link[marray:marrayRaw-class]{marrayRaw}} or
  \code{\link[marray:marrayNorm-class]{marrayNorm}}.
  Equivalently, the name of sample labels are done by the names of the
  slots \code{Slabels}, \code{SampleInfo}, \code{targets} and \code{maTargets@maInfo}.
}

\examples{
## Loading dataset
data(gastro)

## Getting the 'Tissue' label for samples in maigesRaw class object
getLabels(gastro.raw, "Tissue", sLabel=TRUE)

## Getting the 'Type' label for samples in maigesNorm class object
getLabels(gastro.summ, "Type", sLabel=TRUE)

## Getting the 'GeneName' label for genes (spots) in maigesRaw class object
getLabels(gastro.raw, "GeneName", sLabel=FALSE)

## Getting the 'Annot' label for samples in maigesNorm class object
getLabels(gastro.summ, "Annot", sLabel=FALSE)
}

\author{
  Gustavo H. Esteves <gesteves@vision.ime.usp.br>
}

\keyword{array}
