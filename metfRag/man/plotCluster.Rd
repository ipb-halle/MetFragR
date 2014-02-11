\name{plotCluster}
\alias{plotCluster}
\title{Chemical Clustering}
\usage{
  plotCluster(mols, scoreprop="Score", idprop="DatabaseID")
}
\arguments{
  \item{mols}{The a list of rCDK \code{mols}}

  \item{scoreprop}{The name of the property of the
  molecules where the score is kept}

  \item{idprop}{The name of the property of the molecules
  where the database ID is kept}
}
\description{
  Cluster a set of compounds and plot the results
}
\examples{
smiles <- c('CCC', 'CCN', 'CCN(C)(C)',
                   'c1ccccc1Cc1ccccc1',
                   'C1CCC1CC(CN(C)(C))CC(=O)CC')
       mols <- parse.smiles(smiles)
       dummy <- mapply(set.property, mols, "Score", c(1,2,3,4,5))
       dummy <- mapply(set.property, mols, "DatabaseID", c("C1", "C2", "C3", "C4", "C5"))
       plotCluster(mols)
}

