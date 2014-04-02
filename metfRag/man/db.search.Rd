\name{db.search}
\alias{db.search}
\title{
  Find a molecule in a compound database
}
\description{
  This function saves information about a compound in a file in SDK or MOL 
  format. Giving mass, formula or a database id to find a molecule.  
  The search and information was received by PubChem or KEGG.
}
\usage{
  db.search(file, db, type, value)
}
\arguments{
  \item{file}{
    Specify a path and file as character string to save the chemical informations. 
    Using the mol-extension for KEGG and the sdk-extension for PubChem.
  }
  \item{db}{
    This parameter specify the database which was used. This could be
    'kegg' or 'pubchem'
  }
  \item{type}{
    A character string describes the available search methods.
    For Kegg applies 'exact_mass', 'entry', 'formula' and 
    for PubChem applies 'formula', 'mimass', 'exmass'.
  }
  \item{value}{
    The search value for a specified database. mass should be numeric and
    the other types are character.
  }  
}
\value{
  A file with informations about a specified molecule
  is written. If a file was successfully written the function
  returns TRUE else the function returns FALSE.
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
file <- paste(getwd(),"molecule.mol",sep="/");
db   <- 'kegg';
type <- 'exact_mass';
val  <- 174.05;
db.search (file, db, type, val);

file <- paste(getwd(),"molecule.sdf",sep="/");
db   <- 'pubchem';
type <- 'formula';
val  <- 'H2O';
db.search (file, db, type, val);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
