\name{frag.generateFragments}
\alias{frag.generateFragments}
%- Also NEED an '\alias' for EACH other topic documented here
\title{
  Generate fragments of an IAtomContainer object
}
\description{
  \code{frag.generateFragments} generates fragments upto a specified tree depth.
}
\usage{
frag.generateFragments(mol, treeDepth = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mol}{
    molecule as IAtomContainer to be fragmented
  }
  \item{treeDepth}{
    the tree depth of the fragmentation tree
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a list of IAtomContainer objects representing the generated fragments. The number of fragments can be quite large when using this function and complex molecule
  structures. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Christoph Ruttkies
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#get molecule from smiles
smiles <- "CN(C)CC(C1=C=C(C=C1)OC)C2(CCCCC2)O"
mol<-parse.smiles(smiles)[[1]]

#calculate the fragments
fragments<-frag.generateFragments(mol)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fragmentation}
