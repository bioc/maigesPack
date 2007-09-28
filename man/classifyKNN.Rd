\name{classifyKNN}
\alias{classifyKNN}

\title{
  Function to do discrimination analysis
}

\description{
  Function to search by groups of few genes, also called cliques, that
  can discriminate (or classify) between two distinct biological
  sample types, using the k nearest neighbourhood method. This function
  uses exhaustive search.
}

\usage{
classifyKNN(obj=NULL, sLabelID="Classification", facToClass=NULL,
            gNameID="GeneName", geneGrp=1, path=NULL, nGenes=3, kn=5)
}

\arguments{
  \item{obj}{object of class \code{\link{maiges}} to search the classifiers.}
  \item{sLabelID}{character string with the identification of the sample
    label to be used.}
  \item{facToClass}{named list with 2 character vectors specifying the
    samples to be compared. If NULL (default) the first 2 types of
    sLabelID are used.}
  \item{gNameID}{character string with the identification of gene label ID.}
  \item{geneGrp}{character or integer specifying the gene group to be
    tested (\code{colnames} of \code{GeneGrps} slot). If both \code{geneGrp} and
    \code{path} are NULL all genes are used. Defaults to 1 (first group).}
  \item{path}{character or integer specifying the gene network to be
    tested (\code{names} of \code{Paths} slot). If both \code{geneGrp}
    and \code{path} are NULL all genes are used. Defaults to NULL.}
  \item{nGenes}{integer specifying the number of genes in the clique, or
    classifier.}
  \item{kn}{number of neighbours for the \emph{knn} method.}
}

\value{
  The result of this function is an object of class \code{\link{maigesClass}}.
}

\details{
  Pay attention with the arguments \code{geneGrp} and \code{path}, if
  both of them is NULL an exhaustive search for all dataset will be done,
  and this search may be extremely computational intensive, which may
  result in a process during some weeks or months depending on the
  number of genes in your dataset.
  
  If you want to construct classifiers from a group of several genes,
  the \emph{search and choose} (SC) method may be an interesting option. It is
  implemented in the function \code{\link{classifyKNNsc}}.
  This function uses the function \code{\link[class]{knn.cv}} from
  package \emph{class} to construct k-nearest neighbour classifiers. It
  possible to use functions \code{\link{classifyLDA}} or
  \code{\link{classifySVM}} to construct classifiers using Fisher's
  linear discriminant analysis or support vector machines methods, respectively.
}

\seealso{
  \code{\link[class]{knn.cv}}, \code{\link{classifyKNNsc}},
  \code{\link{classifyLDA}}, \code{\link{classifySVM}}.
}

\examples{
## Loading the dataset
data(gastro)

## Doing KNN classifier with 2 genes for the 6th gene group comparing
## the 2 categories from 'Type' sample label.
gastro.class = classifyKNN(gastro.summ, sLabelID="Type",
  gNameID="GeneName", nGenes=2, geneGrp=6)
gastro.class

## To do classifier with 3 genes for the 6th gene group comparing
## normal vs adenocarcinomas from 'Tissue' sample label
gastro.class = classifyKNN(gastro.summ, sLabelID="Tissue",
  gNameID="GeneName", nGenes=3, geneGrp=6,
  facToClass=list(Norm=c("Neso","Nest"), Ade=c("Aeso","Aest")))
}

\author{
  Elier B. Cristo, adapted by Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
