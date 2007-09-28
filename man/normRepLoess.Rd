\name{normRepLoess}
\alias{normRepLoess}

\title{
  Bootstrap of LOWESS normalisation
}

\description{
  This function normalises a microarray object re-doing the LOWESS
  fitting several times, selecting a pre-specified proportion of points
  aleatorily.
}

\usage{
normRepLoess(raw, span=0.4, propLoess=0.5, nRep=50, func="none",
             bkgSub="none", \dots)
}

\arguments{
  \item{raw}{an object of class \code{\link{maigesRaw}} to be normalised.}
  \item{span}{real number in (0,1) representing the proportion of points
    to use in the loess regression.}
  \item{propLoess}{real number in (0,1) representing the proportion of
    points (spots) to be used in each iteration of loess.}
  \item{nRep}{number of repetitions for loess procedure.}
  \item{func}{character string giving the function to estimate the final
    W value. You must use 'mean', 'median' or 'none' (default).}
  \item{bkgSub}{character with background subtraction method, using the
    function \code{\link[limma]{backgroundcorrect}} from \emph{limma}
    package.}
  \item{\dots}{additional parameters for function
    \code{\link[limma:loessfit]{loessFit}} from \emph{limma} package.}
}

\value{
  The result of this function is an object of class \code{\link{maiges}}.
}

\details{
  The LOWESS fitting for normalising microarray data is a computational
  intensive task, so pay attention to not specify a very large argument
  in \code{nRep}. If you do so, your process will take so much time to conclude.
}

\seealso{
  \code{\link[limma:loessfit]{loessFit}}.
}

\examples{
## Loading the dataset
data(gastro)

## Doing the repetition loess with default parameters. Be carefull, this
## is very time consuming
\dontrun{
gastro.norm = normRepLoess(gastro.raw2)
}

## Do the same normalization selecting 60\% dos spots with 10
## repetitions and estimating the W by the mean value.
\dontrun{
gastro.norm = normRepLoess(gastro.raw2, propLoess=0.6, nRep=10, func="mean")
}

}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
