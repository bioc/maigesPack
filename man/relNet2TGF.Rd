\name{relNet2TGF}
\alias{relNet2TGF.maigesRelNetB}
\alias{relNet2TGF.maigesRelNetM}
\alias{relNet2TGF}

\title{
  Transform Relevance Network analysis in TGF output
}

\description{
  This function receive an object of class \code{\link{maigesRelNetB}}
  or \code{\link{maigesRelNetM}} and write TGF files with the relevance
  networks produced.
}

\usage{

relNet2TGF(\dots)

\method{relNet2TGF}{maigesRelNetB}(data, dir = "./",
  filename="group.tgf", corC=NULL, pValue=0.05, \dots)

\method{relNet2TGF}{maigesRelNetM}(data, dir = "./",
  filenames=c("group1.tgf", "group2.tgf", "difPvalue.tgf"),
  pValue=0.05, \dots)

}

\arguments{
  \item{data}{object of class \code{\link{maigesRelNetB}} or \code{\link{maigesRelNetM}}.}
  \item{dir}{character string specifying the folder to save the TGF files.}
  \item{filename}{character string specifying the file name, for objects
    of class \code{\link{maigesRelNetB}}.}
  \item{filenames}{character vector of length 3 with the file names to
    be saved, for objects of class \code{\link{maigesRelNetM}}.}
  \item{corC}{numeric in [0,1] specifying the cutoff for selecting
    absolute correlation. May also be 'max' to select the maximum
    correlation values in a permutation bootstrap strategy, as proposed
    by Butte et al. (2000).}
  \item{pValue}{numeric in [0,1] specifying the cutoff for selecting
    correlation values by p-values.}
  \item{\dots}{additional parameters.}
}

\details{
  This function only picks the result of the \code{\link{relNetworkB}}
  or \code{\link{relNetworkM}} and display write TGF files. This files
  are interesting to be used with \emph{Yed} graph visualisation and
  editing tool, wrote in Java (\url{http://www.yworks.com/en/products_yed_about.htm}). 
}

\value{
  This function don't return any object.
}

\references{
  Butte, A.J.; Tamayo, P.; Slonim, D.; Golub, T.R. and Kohane,
  I.S. Discovering functional relationships between RNA expression and
  chemotherapeutic susceptibility using relevance networks, \bold{PNAS},
  97, 12182-12186, 2000 (\url{http://www.pnas.org/cgi/content/full/97/22/12182})
}

\seealso{
  \code{\link{relNetworkB}}, \code{\link{relNetworkM}},
  \code{\link{maigesRelNetB}}, \code{\link{maigesRelNetM}}.
}

\examples{
## Loading the dataset
data(gastro)

## Constructing the relevance network for sample
## 'Tissue' comparing 'Neso' and 'Aeso' for the 1st gene group
## The same is also true for objects of class maigesRelNetB
gastro.net = relNetworkM(gastro.summ, sLabelID="Tissue", 
  samples = list(Neso="Neso", Aeso="Aeso"), geneGrp=11,
  type="Rpearson")

relNet2TGF(gastro.net)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
