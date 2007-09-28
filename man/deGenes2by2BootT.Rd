\name{deGenes2by2BootT}
\alias{deGenes2by2BootT}

\title{
  Function to do differential expression analysis, comparing only two
  samples
}

\description{
  This function takes an object of class \code{\link{maiges}} and do
  differential expression analysis for the genes onto dataset, comparing
  only two samples by a bootstrap of t statistics method.
}

\usage{
deGenes2by2BootT(data=NULL, sLabelID=names(data@Slabels)[1], sTypeComp=NULL,
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
  This function calculate t statistics and p-values by re-sampling of
  the data using the function \code{\link{bootstrapT}}.

  There is the option to do the t test directly, using the function
  \code{\link{deGenes2by2Ttest}}, or to do the non-parametric Wilcox
  test using the function \code{\link{deGenes2by2Wilcox}}.
}

\seealso{
  \code{\link{bootstrapT}}, \code{\link{deGenes2by2Ttest}} and \code{\link{deGenes2by2Wilcox}}.
}

\examples{
## Loading the dataset
data(gastro)

## Doing bootstrap from t statistic test fot 'Type' sample label, k=1000
## specifies one thousand bootstraps
gastro.boot = deGenes2by2BootT(gastro.summ, sLabelID="Type", k=1000)
gastro.boot
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
