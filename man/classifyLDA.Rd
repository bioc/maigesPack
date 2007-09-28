\name{classifyLDA}
\alias{classifyLDA}

\title{
  Function to do discrimination analysis
}

\description{
  Function to search by groups of few genes, also called cliques, that
  can discriminate (or classify) between two distinct biological
  sample types, using the Fisher's linear discriminat analysis. This
  function uses exhaustive search.
}

\usage{
classifyLDA(obj=NULL, sLabelID="Classification", facToClass=NULL,
            gNameID="GeneName", geneGrp=1, path=NULL, nGenes=3,
            sortBy="cv")
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
  \item{sortBy}{character string with field to sort the result. May
    be 'cv' (default) or 'svd' for cross validation by leave-one-out or
    the singular value decomposition, respectively.}
}

\value{
  The result of this function is an object of class \code{\link{maigesClass}}.
}

\details{
  Pay attention with the arguments \code{geneGrp} and \code{path}, if
  both of them is NULL an exhaustive search for all dataset will be done,
  and this search may be extremely computational intensive, which may
  result in a process running during some weeks or months depending on the
  number of genes in your dataset.
  
  If you want to construct classifiers from a group of several genes,
  the \emph{search and choose} (SC) method may be an interesting option. It is
  implemented in the function \code{\link{classifyLDAsc}}.
  This function uses the function \code{\link[MASS]{lda}} from package
  \emph{MASS} to search by classifiers using Fisher's linear
  discriminant analysis. The functions \code{\link{classifySVM}} and
  \code{\link{classifyKNN}} were also dedined to construct classifiers
  by support vector machines ans k-neighbours, respectively.
}

\seealso{
  \code{\link[MASS]{lda}}, \code{\link{classifySVM}},
  \code{\link{classifyKNN}}, \code{\link{classifyLDAsc}}.
}

\examples{
## Loading the dataset
data(gastro)

## Doing LDA classifier with 2 genes for the 6th gene group comparing
## the 2 categories from 'Type' sample label.
gastro.class = classifyLDA(gastro.summ, sLabelID="Type",
  gNameID="GeneName", nGenes=2, geneGrp=6)
gastro.class

## To do classifier with 3 genes for the 6th gene group comparing
## normal vs adenocarcinomas from 'Tissue' sample label
gastro.class = classifyLDA(gastro.summ, sLabelID="Tissue",
  gNameID="GeneName", nGenes=3, geneGrp=6,
  facToClass=list(Norm=c("Neso","Nest"), Ade=c("Aeso","Aest")))
}

\author{
  Elier B. Cristo, adapted by Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
