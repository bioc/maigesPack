\name{deGenes2by2Wilcox}
\alias{deGenes2by2Wilcox}

\title{
  Function to do differential expression analysis, comparing only two
  samples
}

\description{
  This function takes an object of class \code{\link{maiges}} and do
  differential expression analysis for the genes onto dataset, comparing
  only two samples, by the Wilcox (Mann-Withney) test.
}

\usage{
deGenes2by2Wilcox(data=NULL, sLabelID=names(data@Slabels)[1], sTypeComp=NULL,
                  doClust=TRUE, \dots)
}

\arguments{
  \item{data}{object of class \code{\link{maiges}}.}
  \item{sLabelID}{character string giving the sample label ID to be used.}
  \item{sTypeComp}{list with character vectors specifying the two
    sample types to be compared.}
  \item{doClust}{logical indicating if the object generated from this
    analysis will be used for cluster analysis. Defaults to TRUE.}
  \item{\dots}{additional parameters for functions
    \code{\link[stats]{t.test}}, \code{\link[stats]{wilcox.test}} or
    \code{\link{bootstrapT}}.}
}

\value{
  The result of this function is an object of class
  \code{\link{maigesDE}} if \code{doClust} if FALSE or of class
  \code{\link{maigesDEcluster}} if it is TRUE.
}

\details{
  This function calculate Wilcox statistics and p-values for the test
  comparing the equality of means using the function
  \code{\link[stats]{wilcox.test}}.

  There another function to do parametric t test ,
  \code{\link{deGenes2by2Ttest}}, and another option of non-parametric
  test doing bootstrap of t statistics, \code{\link{deGenes2by2BootT}}.
}

\seealso{
  \code{\link[stats]{wilcox.test}}, \code{\link{deGenes2by2Ttest}} and \code{\link{deGenes2by2BootT}}.
}

\examples{
## Loading the dataset
data(gastro)

## Doing wilcoxon test fot 'Type' sample label
gastro.wilcox = deGenes2by2Wilcox(gastro.summ, sLabelID="Type")
gastro.wilcox
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
