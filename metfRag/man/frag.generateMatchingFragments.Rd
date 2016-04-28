\name{frag.generateMatchingFragments}
\alias{frag.generateMatchingFragments}
%- Also NEED an '\alias' for EACH other topic documented here
\title{
  Generate fragments of an IAtomContainer object matching the MS/MS peak list
}
\description{
  \code{frag.generateMatchingFragments} generates fragments upto a specified tree depth matching the m/z values of a defined MS/MS peak list.
}
\usage{
frag.generateMatchingFragments(molecule, mzs, neutralMonoisotopicMass, mzabs = 0.01, mzppm = 10.0, posCharge = TRUE, ionMode = 1, treeDepth = 2)
}
\arguments{
  \item{molecule}{molecule as IAtomContainer to be fragmented}
  \item{mzs}{m/z peak list to which generated fragments are matched}
  \item{neutralMonoisotopicMass}{neutral monoisotopic mass of the precursor molecule}
  \item{mzabs}{absolute mass deviation used for the peak to fragment matching}
  \item{mzppm}{relative mass deviation used for the peak to fragment matching}
  \item{posCharge}{defines whether ionisation used for MS/MS data acquisition was positive or negative}
  \item{ionMode}{adduct type of the precursor ion (1 -> [M+H], -1 -> [M-H]; 0 -> [M])}
  \item{treeDepth}{the tree depth of the fragmentation tree}
}
\value{
  Returns a list of IAtomContainer objects representing the generated fragments matching to the m/z peak list.
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

#define m/z peak list
mzs<-c(119.051, 123.044, 147.044, 153.019, 179.036, 189.058, 273.076, 274.083)

#calculate the matching fragments
matching.fragments<-frag.generateMatchingFragments(molecule, mzs, 272.06847, mzabs = 0.01, mzppm = 10.0, posCharge = TRUE, ionMode = 1, treeDepth = 2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{in silico fragmentation, MetFrag}
