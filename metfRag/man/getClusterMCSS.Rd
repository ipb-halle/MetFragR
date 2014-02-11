\name{getClusterMCSS}
\alias{getClusterMCSS}
\alias{getMCSS}
\title{Chemical Clustering}
\usage{
  getMCSS(mols)
}
\arguments{
  \item{hcluster}{Previously created hclust for
  \code{mols}}

  \item{mols}{The a list of rCDK \code{mols}}

  \item{h}{dendrogram cutoff #'}
}
\description{
  Calculate the Maximum Common Substructure
}
\examples{
smiles <- c('CCC', 'CCN', 'CCN(C)(C)',
                   'c1ccccc1Cc1ccccc1',
                   'C1CCC1CC(CN(C)(C))CC(=O)CC')
       mols <- parse.smiles(smiles)
       cluster <- metfRag:::hclust.mols(mols)
       clusterreps <- getMCSS(mols, cluster)
}

