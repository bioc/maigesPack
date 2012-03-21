\name{plot}
\alias{plot.maigesRaw}
\alias{plot.maiges}
\alias{plot.maigesANOVA}
\alias{plot.maigesDE}
\alias{plot.maigesDEcluster}
\alias{plot.maigesClass}
\alias{plot.maigesRelNetB}
\alias{plot.maigesRelNetM}
\alias{plot.maigesActMod}
\alias{plot.maigesActNet}
\alias{plot}

\title{
  Method plot for objects defined in this package
}

\description{
  Generic function \code{\link[graphics]{plot}} to display scatter plots
  or other types of graphical representation for objects defined in this
  package.
}

\usage{
\method{plot}{maigesRaw}(x, bkgSub="subtract", z=NULL, legend.func=NULL,
    ylab="W", \dots)

\method{plot}{maiges}(x, z=NULL, legend.func=NULL, ylab="W", \dots)

\method{plot}{maigesANOVA}(x, z=NULL, legend.func=NULL, ylab="W", \dots)

\method{plot}{maigesDE}(x, adjP="none", idx=1, \dots)

\method{plot}{maigesDEcluster}(x, adjP="none", idx=1, \dots)

\method{plot}{maigesClass}(x, idx=1, \dots)

\method{plot}{maigesRelNetB}(x=NULL, cutPval=0.05, cutCor=NULL,
name=NULL, \dots)

\method{plot}{maigesRelNetM}(x=NULL, cutPval=0.05, names=NULL, \dots)

\method{plot}{maigesActMod}(x, type=c("S", "C")[2], keepEmpty=FALSE, \dots)

\method{plot}{maigesActNet}(x, type=c("score", "p-value")[1], \dots)
}

\arguments{
  \item{x}{an object of any class defined in this package, except \code{\link{maigesPreRaw}}.}
  \item{bkgSub}{string specifying the method for background
    subtraction. See function \code{\link[limma]{backgroundcorrect}} to
    find the available options.}
  \item{z}{accessor method for stratifying data, see \code{\link[marray]{maPlot}}.}
  \item{legend.func}{string specifying options to show legend in the figure.}
  \item{ylab}{character string specifying the label to y axis.}
  \item{adjP}{type of p-value adjustment, see function
    \code{\link[multtest]{mt.rawp2adjp}} in package multtest.}
  \item{idx}{index of the test statistic to be plotted in case of
    objects of classes \code{\link{maigesDE}} and
    \code{\link{maigesDEcluster}} or the index of the clique to be plotted
    in case of object with class \code{\link{maigesClass}}.}
  \item{cutPval}{real number in [0,1] specifying a cutoff p-value to
    show significant results from relevance network analysis. For class
    \code{\link{maigesRelNetB}}, if this parameter is specified the
    argument \code{cutCor} isn't used.}
  \item{cutCor}{real number in [0,1], specifying a coefficient
    correlation value cutoff (in absolute value) to show only absolute
    correlation values greater than this value. Pay attention, to use
    this cutoff it is necessary to specify \code{cutPval} as NULL.}
  \item{name}{character string giving a name for sample type tested to
    be plotted as a name in the method for class \code{\link{maigesRelNetB}}.}
  \item{names}{similar to the previous one, but it is a vector of length
    3.}
  \item{type}{string specifying the type of colour map to be plotted. For
    class \code{\link{maigesActMod}} it must be 'S' or 'C' for samples
    or biological conditions, respectively. For class
    \code{\link{maigesActNet}} it must be 'score' or 'p-value' for the
    statistics or p-values of the tests, respectively.}
  \item{keepEmpty}{logical, if true the results of all gene groups are
    displayed, else only the gene groups that present at least one
    significant result are displayed.}
  \item{\dots}{additional arguments for method
    \code{\link[marray]{maPlot}} or \code{\link[graphics]{plot}}}
}

\details{
  This method uses the function \code{\link[marray]{maPlot}} to display
  scatter plots ratio vs mean values for objects of class
  \code{\link{maiges}}, \code{\link{maigesRaw}} or
  \code{\link{maigesANOVA}}. For objects of class \code{\link{maigesDE}}
  or \code{\link{maigesDEcluster}}, this method display volcano
  plots. For objects of class \code{\link{maigesClass}} it do 2 or 3
  dimensions scatter plots that facilitate the visualisation of good
  classifying cliques of genes For objects of class
  \code{\link{maigesRelNetM}} the method displays 3 circular graphs
  representing the correlation values for the two groups tested and the
  p-values of the tests. For class \code{\link{maigesRelNetB}} it
  displays only one circular graph showing the correlation values for
  the type tested. In objects of class \code{\link{maigesActMod}} and
  \code{\link{maigesActNet}} the method do the same job as
  \code{\link{image}}.

  Pay attention that, even using the method \code{\link[marray]{maPlot}}
  from \emph{marray} package, we plot \emph{W} values against \emph{A}
  values instead of \emph{MA} plots.
  
}

\author{
  Gustavo H. Esteves <gesteves@vision.ime.usp.br>
}

\examples{
## Loading the dataset
data(gastro)

## Example with an object of class maigesRaw, without and with backgound
## subtraction, also we present a plot with normexp (from limma package)
## subtract algorithm.
plot(gastro.raw[,1], bkgSub="none")
plot(gastro.raw[,1], bkgSub="subtract")
plot(gastro.raw[,1], bkgSub="normexp")

## Example with an object of class maigesNorm.
plot(gastro.norm[,1])



## Example for objects of class maigesDE.

## Doing bootstrap from t statistic test fot 'Type' sample label, k=1000
## specifies one thousand bootstraps
gastro.ttest = deGenes2by2Ttest(gastro.summ, sLabelID="Type")

plot(gastro.ttest) ## Volcano plot


## Example for object of class maigesClass.

## Doing LDA classifier with 3 genes for the 6th gene group comparing
## the 2 categories from 'Type' sample label.
gastro.class = classifyLDA(gastro.summ, sLabelID="Type",
  gNameID="GeneName", nGenes=3, geneGrp=6)

\dontrun{
plot(gastro.class) ## plot the 1st classifier
plot(gastro.class, idx=7) ## plot the 7th classifier
}


## Example for object of class maigesActNet

## Doing functional classification of gene groups for 'Tissue' sample label
gastro.mod = activeMod(gastro.summ, sLabelID="Tissue", cutExp=1,
  cutPhiper=0.05)

plot(gastro.mod, "S", margins=c(15,3)) ## Plot for individual samples
plot(gastro.mod, "C", margins=c(21,5)) ## Plot for unique biological conditions



## Example for object of class maigesRelNetB

## Constructing the relevance network (Butte's method) for sample
## 'Tissue' equal to 'Neso' for the 1st gene group
gastro.net = relNetworkB(gastro.summ, sLabelID="Tissue", 
  samples="Neso", geneGrp=1, type="Rpearson")

plot(gastro.net, cutPval=0.05)




## Example for object of class maigesRelNetM

## Constructing the relevance network for sample
## 'Tissue' comparing 'Neso' and 'Aeso' for the 1st gene group
gastro.net = relNetworkM(gastro.summ, sLabelID="Tissue", 
  samples = list(Neso="Neso", Aeso="Aeso"), geneGrp=11,
  type="Rpearson")

plot(gastro.net, cutPval=0.05)
plot(gastro.net, cutPval=0.01)



## Example for objects of class maigesActNet

## Doing functional classification of gene networks for sample Label
## given by 'Tissue'
gastro.net = activeNet(gastro.summ, sLabelID="Tissue")

plot(gastro.net, type="score", margins=c(21,5))
plot(gastro.net, type="p-value", margins=c(21,5))
}

\seealso{
  \code{\link[multtest]{mt.rawp2adjp}},
  \code{\link[limma]{backgroundcorrect}}, \code{\link[marray]{maPlot}} in
  the package marray, \code{\link[graphics]{plot}} in the base package.
}

\keyword{array}
