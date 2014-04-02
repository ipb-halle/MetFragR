\name{db.kegg.MoleculeToFile}
\alias{db.kegg.MoleculeToFile}
\title{
  Saves Mol-File by KEGG C Numbers
}
\description{
  This function receives data from KEGG by KEGG Entrys and 
  write the result in a mol-file.
}
\usage{
  db.kegg.MoleculeToFile(file, ids)
}
\arguments{
  \item{file}{
    A character string describes a path to mol-file.
  }
  \item{ids}{
    A matrix with one column contains KEGG Entrys.
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
ids <- c("cpd:C00493","cpd:C04236","cpd:C16588");
mat <- matrix(c(ids),3);
file <- paste(getwd(),"molecule.mol",sep="/");
db.kegg.MoleculeToFile(file, mat);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

