\name{db.kegg.getId}
\alias{db.kegg.getId}
\title{
  Get Compound ID from KEGG
}
\description{
  This function returns a KEGG Compound ID by exact mass or formula.
}
\usage{
  db.kegg.getId(type, value)
}
\arguments{
  \item{type}{
    A character string describes the available search methods.
    For Kegg applies 'exact_mass' or 'formula'.
  }
  \item{value}{
    The search value for a specified database. The mass should be numeric and
    the other types are characters. 
  }  
}
\value{
  KEGG Entrys are returned.
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
type <- 'exact_mass';
val  <- 174.05;
db.kegg.getId(type, val);

type <- 'formula';
val  <- 'H2O';
db.kegg.getId(type, val);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
