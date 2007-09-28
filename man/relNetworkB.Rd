\name{relNetworkB}
\alias{relNetworkB}

\title{
  Relevance Network analysis
}

\description{
  Function to construct Relevance Networks for one biological type
  (Butte's Relevance Network).
}

\usage{
relNetworkB(data=NULL, gLabelID="GeneName", sLabelID="Classification",
            geneGrp=NULL, path=NULL, samples=NULL,
            type="Rpearson", bRep=1000, \dots)
}

\arguments{
  \item{data}{object of class \code{\link{maiges}}.}
  \item{gLabelID}{character string giving the identification of gene label ID.}
  \item{sLabelID}{character string giving the identification of sample label ID.}
  \item{geneGrp}{character string (or numeric index) specifying the gene
    group to calculate the correlation values between them. If
    NULL (together with path) all genes are used.}
  \item{path}{character string (or numeric index) specifying the gene network
    to calculate the correlation values between them. If NULL (together
    with geneGrp) all genes are used.}
  \item{samples}{a character vector specifying the group to be compared.}
  \item{type}{type of correlation to be calculated. May be 'Rpearson'
    (default), 'pearson', 'kendall', 'spearman' or 'MI'.}
  \item{bRep}{integer specifying the number of bootstrap permutation to
    calculate the significance of correlation values.}
  \item{\dots}{additional parameters for functions
    \code{\link{robustCorr}} or \code{\link[stats]{cor}}.}
}

\value{
  The result of this function is an object of class \code{\link{maigesRelNetB}}.
}

\details{
  This method uses the function \code{\link[stats]{cor}} to calculate
  the usual correlation values, \code{\link{robustCorr}} to calculate
  a robust correlation using an idea similar to the leave-one-out or
  \code{\link{MI}} to calculate mutual information values.
}

\seealso{
  \code{\link[stats]{cor}}, \code{\link{robustCorr}}, \code{\link{MI}}
  \code{\link{maigesRelNetB}}, \code{\link{plot.maigesRelNetB}},
  \code{\link{image.maigesRelNetB}}.
}

\references{
  Butte, A.J. and Kohane, I.S. Unsupervised Knowledge discovery in
  medical databases using relevance networks. In Proc. AMIA Symp.,
  711-715, 1999 (\url{http://www.amia.org/pubs/symposia/D005550.HTM})

  Butte, A.J.; Tamayo, P.; Slonim, D.; Golub, T.R. and Kohane,
  I.S. Discovering functional relationships between RNA expression and
  chemotherapeutic susceptibility using relevance networks, \bold{PNAS},
  97, 12182-12186, 2000
  (\url{http://www.pnas.org/cgi/content/full/97/22/12182})

  Butte, A.J. and Kohane, I.S. Mutual information relevance networks:
  functional genomic clustering using pairwise entropy measurements. In
  Pacific Symposium on Biocomputing, 5, 415-426, 2000
  (\url{http://psb.stanford.edu/psb-online/proceedings/psb00/})
}

\examples{
## Loading the dataset
data(gastro)

## Constructing the relevance network (Butte's method) for sample
## 'Tissue' equal to 'Neso' for the 1st gene group
gastro.net = relNetworkB(gastro.summ, sLabelID="Tissue", 
  samples="Neso", geneGrp=1, type="Rpearson")

## Constructing the relevance network (Butte's method) for sample
## 'Type' equal to 'Col' for the 1st gene group using the conventional
## pearson correlation
gastro.net = relNetworkB(gastro.summ, sLabelID="Type", 
  samples="Col", geneGrp=1, type="pearson")
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
