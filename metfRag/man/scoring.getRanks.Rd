\name{scoring.getRanks}
\alias{scoring.getRanks}
%- Also NEED an '\alias' for EACH other topic documented here
\title{
  Rank Molecules by IAtomContainer properties
}
\description{
  \code{scoring.getRanks} ranks molecules in IAtomContainer format
  by specific properties or an own function.
}
\usage{
scoring.getRanks(mols, sorting, condition, split=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mols}{
    molecules as IAtomContainer in \code{list} format 
  }
  \item{sorting}{
    \code{list} of properties in molecules or own functions. 
    Hence, the order of the properties represents the order of sorting
    the container.
  }
  \item{condition}{
    a property \code{list} which describes the ranked molecule. The parameter name
    is the property name and the parameter value is the (splitted) property value
  }
  \item{split}{
    a \code{list} with \code{sep} which describes the splitting parameter
    and \code{pos} which represents the position number of splitted part
    of molecule.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returned a list with optimitistic rank, pessimistic rank and relative rank
  inclusive the better, worst, equal and all candidates.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Adrian Helmchen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#get molecules
params <- list(mass=174.05, range=0.001)
conts <- db.pubchem.getId(params)
mols.pc <- db.pubchem.getMoleculeContainer(conts)

#score molecules
mzs <-c(119.051,123.044,147.044,153.019,179.036,189.058,273.076,274.083)
ints <-c(467.616,370.662,6078.145,10000.0,141.192,176.358,10000.000,318.003)
mols.pc <- score.molecules.from.container(mols.pc, mzs, ints, 272.06847)

#rank given molecule
params <- list('Score', list('Score', function(x) x * 5), 'de.ipbhalle.metfrag.score.PeakScore')
condition <- list(IUPAC_INCHIKEY="YOGIRYKSHRZUBS")
scoring.getRanks(mols.pc, params, condition, list(sep="-", pos=1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ranking}
