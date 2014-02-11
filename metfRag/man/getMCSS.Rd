\name{getMCSS}
\alias{getMCSS}
\title{Chemical Clustering}
\usage{
  getMCSS(mols)
}
\arguments{
  \item{mols}{The a list of rCDK \code{mols}}
}
\description{
  Calculate the Maximum Common Substructure
}
\examples{
smiles <- c('CCC', 'CCN', 'CCN(C)(C)',
                   'c1ccccc1Cc1ccccc1',
                   'C1CCC1CC(CN(C)(C))CC(=O)CC')
       mols <- parse.smiles(smiles)
       dummy <- getMCSS(mols)
}

