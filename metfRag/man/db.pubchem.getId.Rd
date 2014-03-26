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
  db.kegg.getId(type, value)
}
\arguments{
  \item{type}{
    A character string describes the available search methods.
    For PubChem applies 'mimass' for monoisotopic mass, 'exmass' for
    exact mass, or 'formula'.
  }
  \item{value}{
    The search value for a specified database. Both, the exact mass and the monoisotopic 
    mass should be numeric and the other types are characters. 
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
type <- 'mimass';
val  <- 174.05;
db.pubchem.getId(type, val);

type <- 'exmass';
val  <- 174.05;
db.pubchem.getId(type, val);

type <- 'formula';
val  <- 'H2O';
db.pubchem.getId(type, val);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

