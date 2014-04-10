\name{db.kegg.getId}
\alias{db.kegg.getId}
\title{
  Get Compound ID from KEGG
}
\description{
  This function returns a KEGG Compound ID by exact mass or formula.
}
\usage{
  db.kegg.getId(seek)
}
\arguments{
  \item{seek}{
    A character string describes the available search methods.
    Using KEGG the types are 'exact_mass', 'entry' and 'formula'.
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
params <- list(mass=174.05,range=0);
db.kegg.getId(params);

params <- list(formula='H2O')
db.kegg.getId(params);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
