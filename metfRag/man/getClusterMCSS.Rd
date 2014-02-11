\name{getClusterMCSS}
\alias{getClusterMCSS}
\title{Chemical Clustering}
\usage{
  getClusterMCSS(hcluster, mols, h = 0.2)
}
\arguments{
  \item{hcluster}{Previously created hclust for
  \code{mols}}

  \item{mols}{The a list of rCDK \code{mols}}

  \item{h}{dendrogram cutoff}
}
\description{
  Calculate the Maximum Common Substructure
}
\examples{
library(rcdk)
       smiles <- c('CCC', 'CCN', 'CCN(C)(C)',
                   'c1ccccc1Cc1ccccc1',
                   'C1CCC1CC(CN(C)(C))CC(=O)CC')
       mols <- parse.smiles(smiles)
       cluster <- hclust.mols(mols)
       clusterreps <- getClusterMCSS(cluster, mols)
}

