\name{db.kegg.getMoleculeContainer}
\alias{db.kegg.getMoleculeContainer}
\title{
  Saves Mol-File by KEGG C Numbers
}
\description{
  This function receives data from KEGG by KEGG Entrys and 
  either write the result in a mol-file or return
  an IAtomContainer.
}
\usage{
  db.kegg.getMoleculeContainer(ids)
}
\arguments{
  \item{ids}{
    A matrix with one column contains KEGG Entrys.
  }  
}
\value{
  The function returns either a boolean whether a file was
  successfully written or a list of IAtomContainer.
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
ids <- c("C00493","C04236","C16588");
mat <- matrix(c(ids),3);
db.kegg.getMoleculeContainer(mat);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

