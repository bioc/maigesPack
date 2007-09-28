\name{gastro}
\alias{gastro}
\alias{gastro.raw}
\alias{gastro.raw2}
\alias{gastro.norm}
\alias{gastro.summ}

\title{
  Gastro-esophagic dataset
}

\description{
  This dataset is composed of five objects of classes
  \code{\link{maigesPreRaw}}, \code{\link{maigesRaw}} and
  \code{\link{maiges}} containing a piece of the dataset that was
  analysed by the MAIGES (\url{http://www.maiges.org/}) and published at
  \emph{Cancer Research} (see the reference below). In this dataset we
  evaluated the profiles of expression for 72 gastro-esophagical
  samples. Were used cDNA microarrays with slides containing 4800 spots
  representing approximately 4400 unique genes. The original dataset was
  composed from 86 samples with dye swap (what makes 172 chips). The
  dataset available here contain a subgroup of 20 samples (40 chips) and
  500 spots (486 unique genes).

  The objects available are:
  - \code{gastro}: \code{\link{maigesPreRaw}} class containing the raw dataset.
  - \code{gastro.raw}: \code{\link{maigesRaw}} class containing the
  dataset ready for filtration and normalization.
  - \code{gastro.raw2}: same to the above, after filtration.
  - \code{gastro.norm}: \code{\link{maiges}} containing the normalized
  data.
  - \code{gastro.summ}: same to the above, after summary of replicates.
}

\usage{
data(gastro)
}

\references{
  Gomes, L.I.; Esteves, G.H.; Carvalho, A.F.; Cristo, E.B.; Hirata Jr.,
  R.; Martins, W.K.; Marques, S.M.; Camargo, L.P.; Brentani, H.;
  Pelosof, A.; Zitron, C.; Sallum, R.A.; Montagnini, A.; Soares, F.A.;
  Neves, E.J. & Reis, L.F. Expression Profile of Malignant and
  Nonmalignant Lesions of Esophagus and Stomach: Differential Activity
  of Functional Modules Related to Inflammation and Lipid Metabolism,
  \bold{Cancer Research}, 65, 7127-7136, 2005
  (\url{http://cancerres.aacrjournals.org/cgi/content/abstract/65/16/7127})
}

\keyword{datasets}
