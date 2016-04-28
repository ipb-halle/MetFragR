\name{hclust.mols}
\alias{hclust.mols}
\title{Chemical Clustering of MetFrag and MetFusion results SDF files}
\usage{
  hclust.mols(mols = NULL, smiles = NULL, filename = NULL,
    scoreprop = "Score", idprop = "Identifier")
}
\arguments{
  \item{mols}{The a list of rCDK \code{mols}.}

  \item{smiles}{A character vector of SMILES codes.}

  \item{filename}{A filename to an SD File. If boths
  \code{mols} and \code{filename} are given, filename is
  ignored. Files are read by \code{cdk::load.molecules()}}

  \item{scoreprop}{The name of the property of the
  molecules where the score is kept}

  \item{idprop}{The name of the property of the molecules
  where the database ID is kept}
}
\description{
  Cluster a set of compounds ranked with MetFrag and
  MetFusion
}

