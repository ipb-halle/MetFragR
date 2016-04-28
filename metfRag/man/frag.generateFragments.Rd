\name{frag.generateFragments}
\alias{frag.generateFragments}
\title{
  Generate fragments of an IAtomContainer object
}
\description{
  \code{frag.generateFragments} generates fragments upto a specified tree depth.
}
\usage{
frag.generateFragments(molecule, treeDepth = 2)
}
\arguments{
  \item{molecule}{
    molecule as IAtomContainer to be fragmented
  }
  \item{treeDepth}{
    the tree depth of the fragmentation tree
  }
}
\value{
  Returns a list of IAtomContainer objects representing the generated fragments. The number of fragments can be quite large when using this function and complex molecule
  structures. 
}
\references{
Ruttkies C., Schymanski E.L. et al, MetFrag relaunched: incorporating strategies beyond in silico fragmentation. Journal of Cheminformatics, 2016, 8:3.
}
\author{
  Christoph Ruttkies
}
\examples{
#get molecule from smiles
smiles <- "CN(C)CC(C1=C=C(C=C1)OC)C2(CCCCC2)O"
molecule<-parse.smiles(smiles)[[1]]

#calculate the fragments
fragments<-frag.generateFragments(molecule)

}
\keyword{in silico fragmentation, MetFrag}
