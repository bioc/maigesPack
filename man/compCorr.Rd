\name{compCorr}
\alias{compCorr}

\title{
  Compute correlation differences and their p-values
}

\description{
  This function takes two correlation values (or matrices of correlation
  values) and calculate the differences between these values (term by
  term) and their respective p-values by a Fisher's Z transformation.
}

\usage{
compCorr(n1, r1, n2, r2)
}

\arguments{
  \item{n1}{numerical or matrix of sample sizes for group 1.}
  \item{r1}{numerical or matrix of correlation values for group 1.}
  \item{n2}{numerical or matrix of sample sizes for group 2.}
  \item{r2}{numerical or matrix of correlation values for group 2.}
}

\value{
  The result of this function is a list with two numerical items.
  \item{diff}{matrix (or a single number) of differences between
    correlation values from two groups}
  \item{pval}{matrix (or a single number) of p-values of differences
    between the correlation values}
}

\details{
  This function use Fisher's Z transformation from scripts adapted from
  the Internet: \cr
  \url{http://ftp.sas.com/techsup/download/stat/compcorr.html} \cr
  \url{http://www.fon.hum.uva.nl/Service/Statistics/Two_Correlations.html}
}

\examples{
compCorr(12, 0.9, 16, 0.73)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
