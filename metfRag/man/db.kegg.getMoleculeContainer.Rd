\name{db.kegg.getMoleculeContainer}
\alias{db.kegg.getMoleculeContainer}
\title{
  Get AtomContainer by KEGG-IDs
}
\description{
  This function receives data from KEGG by KEGG Entrys and 
  return IAtomContainer.
}
\usage{
  db.kegg.getMoleculeContainer(ids)
}
\arguments{
  \item{ids}{
    A vector with KEGG Entrys.
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
mol.ids <- db.kegg.getId(params);
db.kegg.getMoleculeContainer(mol.ids);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

