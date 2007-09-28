\name{relNetworkM}
\alias{relNetworkM}

\title{
  Relevance Network analysis
}

\description{
  Function to construct Relevance Networks comparing two distinct
  biological types.
}

\usage{
relNetworkM(data=NULL, gLabelID="GeneName", sLabelID="Classification",
            geneGrp=NULL, path=NULL, samples=NULL,
            type="Rpearson", \dots)
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
  \item{samples}{a named list with two character vectors specifying the
    two groups that must be compared.}
  \item{type}{type of correlation to be calculated. May be 'Rpearson'
    %(default), 'pearson', 'kendall', 'spearman' or 'MI'.}
    (default), 'pearson', 'kendall' or 'spearman'.}
  \item{\dots}{additional parameters for functions
    \code{\link{robustCorr}} or \code{\link[stats]{cor}}.}
}

\value{
  The result of this function is an object of class \code{\link{maigesRelNetM}}.
}

\details{
  This method uses the function \code{\link[stats]{cor}} to calculate
  %the usual correlation values, \code{\link{robustCorr}} to calculate
  the usual correlation values or \code{\link{robustCorr}} to calculate
  %a robust correlation using an idea similar to the leave-one-out or
  a robust correlation using an idea similar to the leave-one-out.
  %\code{\link{MI}} to calculate the mutual information values. The
  The
  correlation values are calculated for pairs of genes in the two groups
  specified by the argument \code{samples}, then a Fisher's Z
  transformation are done to calculate the significance for the
  difference between the two correlation values, this is implemented in
  the function \code{\link{compCorr}}. This method was first used in the
  work from Gomes et al. (2005).
}

\seealso{
  %\code{\link[stats]{cor}}, \code{\link{robustCorr}}, \code{\link{MI}}
  \code{\link[stats]{cor}}, \code{\link{robustCorr}}
  \code{\link{compCorr}}, \code{\link{maigesRelNetM}},
  \code{\link{plot.maigesRelNetM}},
  \code{\link{image.maigesRelNetM}}.
}

\references{
  Gomes, L.I.; Esteves, G.H.; Carvalho, A.F.; Cristo, E.B.; Hirata Jr.,
  R.; Martins, W.K.; Marques, S.M.; Camargo, L.P.; Brentani, H.;
  Pelosof, A.; Zitron, C.; Sallum, R.A.; Montagnini, A.; Soares, F.A.;
  Neves, E.J. & Reis, L.F. Expression Profile of Malignant and
  Nonmalignant Lesions of Esophagus and Stomach: Differential Activity
  of Functional Modules Related to Inflammation and Lipid Metabolism,
  \bold{Cancer Research}, 65, 7127-7136, 2005
  (\url{http://cancerres.aacrjournals.org/cgi/content/abstract/65/16/7127})
}

\examples{
## Loading the dataset
data(gastro)

## Constructing the relevance network for sample
## 'Tissue' comparing 'Neso' and 'Aeso' for the 1st gene group
gastro.net = relNetworkM(gastro.summ, sLabelID="Tissue", 
  samples = list(Neso="Neso", Aeso="Aeso"), geneGrp=11,
  type="Rpearson")
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
