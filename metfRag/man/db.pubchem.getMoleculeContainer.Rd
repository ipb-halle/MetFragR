\name{db.pubchem.getMoleculeContainer}
\alias{db.pubchem.getMoleculeContainer}
\title{
  Get AtomContainer by PubChem-IDs
}
\description{
  This function receives data from PubChem by PubChem Entrys and 
  return IAtomContainer.
}
\usage{
  db.pubchem.getMoleculeContainer(ids)
}
\arguments{
  \item{ids}{
    A vector with PubChem Entrys.
  }  
}
\value{
  The function returns a list of IAtomContainer.
}
\author{
  Adrian Helmchen
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
params <- list(mass=174.05, range=0.001);
mol.ids <- db.pubchem.getId(params);
db.pubchem.getMoleculeContainer(mol.ids);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
