\name{db.kegg.getId}
\alias{db.kegg.getId}
\title{
  Get Compound ID from PubChem
}
\description{
  This function returns a PubChem Compound ID by exact mass, 
  monoisotopic mass, formula or name.
}
\usage{
  db.kegg.getId(seek)
}
\arguments{
  \item{seek}{
    The search value for a specified database. Both, the exact mass 
    and the monoisotopic mass should be numeric and the other types 
    are characters. 
  }  
}
\value{
  PubChem IDs are returned.
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
params <- list(mass=174.05, range=0);
db.pubchem.getId(params);

params <- list(formula='H2O');
db.pubchem.getId(params);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

