\name{db.pubchem.MoleculeToFile}
\alias{db.pubchem.MoleculeToFile}
\title{
  Saves SDF-File by PubChem Ids
}
\description{
  This function receives data from PubChem by PubChem IDs and save the result
  in a sdf-file.
}
\usage{
  db.pubchem.MoleculeToFile(file, ids)
}
\arguments{
  \item{file}{
    A character string describes a path to sdf-file.
  }
  \item{ids}{
    A matrix with one column contains PubChem IDs.
  }  
}
\value{
  The function returns a boolean whether a file was
  successfully written.
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
ids <- c("73264919","73264709","73264708");
mat <- matrix(c(ids),3);
file <- paste(getwd(),"molecule.sdf",sep="/");
db.pubchem.MoleculeToFile(file, mat);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
