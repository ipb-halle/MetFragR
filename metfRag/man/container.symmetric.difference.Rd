\name{container.symmetric.difference}
\alias{container.symmetric.difference}
\title{
Symmetric difference of two sets of IAtomContainer
}
\description{
This function build the symmetric difference of two sets of IAtomContainer.
}
\usage{
container.symmetric.difference(set.a, set.b, link)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{set.a}{a \code{list} of IAtomContainern}
  \item{set.b}{a \code{list} of IAtomContainern}
  \item{link}{a \code{string} which is a container property. If your leave this
  blank \code{comm.lib.showLinkOptions} are called and show you a list
  which lists the container properties}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A \code{list} with combined sets of IAtomContainer
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Adrian Helmchen
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{help}}
%}
\examples{
params <- list(mass=174.05, range=0.001)
mol.ids <- db.kegg.getId(params)
mol.cont <- db.kegg.getMoleculeContainer(mol.ids)
container.symmetric.difference(mol.cont, mol.cont, "FORMULA");
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{symmetric difference}