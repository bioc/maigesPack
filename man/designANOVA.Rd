\name{designANOVA}
\alias{designANOVA}

\title{
  Function to construct design an contrasts matrices for ANOVA models
}

\description{
  This function takes an object of class \code{\link{maiges}} together
  with other arguments and construct the matrices of design and
  contrasts for adjusting ANOVA models. The design matrix are generated
  using the function \code{\link[stats]{model.matrix}}.
}

\usage{
designANOVA(data=NULL, factors=names(data@Slabels), model=NULL,
            contrasts=NULL, \dots)
}

\arguments{
  \item{data}{object of class \code{\link{maiges}}.}
  \item{factors}{vector of character strings specifying the sample
    labels IDs to be used as factors in the models.}
  \item{model}{a formula specifying the model to be fitted.}
  \item{contrasts}{character vector specifying the contrasts to be
    done. This is done by the function \code{\link[limma]{makeContrasts}}
    from package \emph{limma}. Pay attention that we use the
    treatment-control parametrisation.}
  \item{\dots}{additional parameters for function \code{\link[stats]{model.matrix}}.}
}

\value{
  The result of this function is an object of class \code{\link{maigesANOVA}}.
}

\seealso{
  \code{\link[stats]{model.matrix}}, \code{\link[limma]{makeContrasts}},
  \code{\link{deGenesANOVA}}.
}

\examples{
## Loading the dataset
data(gastro)

## Constructing a maigesANOVA object for the 'Tissue' sample label using
## default model (simple linear model with intercept) and contrasts (all
## parameters are equal between themselves)
gastro.ANOVA = designANOVA(gastro.summ, factors="Tissue")
gastro.ANOVA
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
