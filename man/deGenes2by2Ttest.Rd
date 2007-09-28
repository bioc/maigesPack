\name{deGenes2by2Ttest}
\alias{deGenes2by2Ttest}

\title{
  Function to do differential expression analysis, comparing only two
  samples
}

\description{
  This function takes an object of class \code{\link{maiges}} and do
  differential expression analysis for the genes onto dataset, comparing
  only two samples, by t test.
}

\usage{
deGenes2by2Ttest(data=NULL, sLabelID=names(data@Slabels)[1], sTypeComp=NULL,
                 doClust=TRUE, \dots)
}

\arguments{
  \item{data}{object of class \code{\link{maiges}}.}
  \item{sLabelID}{character string giving the sample label ID to be used.}
  \item{sTypeComp}{list with character vectors specifying the two
    sample types to be compared.}
  \item{doClust}{logical indicating if the object generated from this
    analysis will be used for cluster analysis. Defaults to TRUE.}
  \item{\dots}{additional parameters for function \code{\link[stats]{t.test}}.}
}

\value{
  The result of this function is an object of class
  \code{\link{maigesDE}} if \code{doClust} if FALSE or of class
  \code{\link{maigesDEcluster}} if it is TRUE.
}

\details{
  This function calculate t statistics and p-values for the test of
  difference in the means of the groups using the function
  \code{\link[stats]{t.test}}.

  There are other two functions, \code{\link{deGenes2by2Wilcox}} and
  \code{\link{deGenes2by2BootT}}, to do non-parametric tests by the
  Wilcox test an t statistic bootstrap, respectively.
}

\seealso{
  \code{\link[stats]{t.test}}, \code{\link{deGenes2by2Wilcox}} and \code{\link{deGenes2by2BootT}}.
}

\examples{
## Loading the dataset
data(gastro)

## Doing t test fot 'Type' sample label
gastro.ttest = deGenes2by2Ttest(gastro.summ, sLabelID="Type")
gastro.ttest
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
