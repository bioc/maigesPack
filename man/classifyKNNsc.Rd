\name{classifyKNNsc}
\alias{classifyKNNsc}

\title{
  Function to do discrimination analysis, by the search and choose method
}

\description{
  Function to search by groups of few genes, also called cliques, that
  can discriminate (or classify) between two distinct biological
  sample types, using the k nearest neighbours method. This function
  uses the search and choose method.
}

\usage{
classifyKNNsc(obj=NULL, sLabelID="Classification", func="wilcox.test",
              facToClass=NULL, gNameID="GeneName", geneGrp=1, path=NULL,
              nGenes=3, cliques=100, kn=5)
}

\arguments{
  \item{obj}{object of class \code{\link{maiges}} to search the classifiers.}
  \item{sLabelID}{character string with the identification of the sample
    label to be used.}
  \item{func}{string specifying the function to be used to search by the
    initial one-dimensional classifiers, like 'wilcox.test' or 't.test'.}
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
  \item{cliques}{integer specifying the number of cliques or classifiers
    to be generated.}
  \item{kn}{number of neighbours for the \emph{knn} method.}
}

\value{
  The result of this function is an object of class \code{\link{maigesClass}}.
}

\details{
  This function implements the method known as Search and choose
  proposed by Cristo (2003). If you want to use an exhaustive search use
  the function \code{\link{classifyKNN}}.
  
  This function uses the function \code{\link[class]{knn.cv}} from
  package \emph{class} to construct k-nearest neighbour classifiers. It
  is possible to use the functions \code{\link{classifyLDAsc}} and
  \code{\link{classifySVMsc}} to search by classifiers using Fisher's
  linear discriminant analysis and support vector machinnes, respectively.
}

\references{
  Cristo, E.B. Metodos Estatisticos na Analise de Experimentos de
  Microarray. Masther's thesis, Instituto de Matematica
  e Estatistica - Universidade de Sao Paulo, 2003 (in portuguese).
}

\seealso{
  \code{\link[class]{knn.cv}}, \code{\link{classifyKNN}},
  \code{\link{classifyLDAsc}} and \code{\link{classifySVMsc}}.
}

\examples{
## Loading the dataset
data(gastro)

## Doing KNN classifier with 2 genes for the 6th gene group comparing
## the 2 categories from 'Type' sample label.
gastro.class = classifyKNNsc(gastro.summ, sLabelID="Type",
  gNameID="GeneName", nGenes=2, geneGrp=1, cliques=10)
gastro.class

## To do classifier with 3 genes for the 6th gene group comparing
## normal vs adenocarcinomas from 'Tissue' sample label
gastro.class = classifyKNNsc(gastro.summ, sLabelID="Tissue",
  gNameID="GeneName", nGenes=3, geneGrp=1, cliques=10,
  facToClass=list(Norm=c("Neso","Nest"), Ade=c("Aeso","Aest")))
}

\author{
  Elier B. Cristo, adapted by Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
