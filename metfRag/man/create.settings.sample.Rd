\name{create.settings.sample}
\alias{create.settings.sample}
\title{Create a settings object used by run.metfrag function}
\usage{
	create.settings.sample()
}
\arguments{
  \item{settingsObject}{A \code{list} of parameter settings}
}
\description{
  Returns a list with defined settings.
}
\value{
Returns a list with defined settings.
}
\references{
Ruttkies C., Schymanski E.L. et al, MetFrag relaunched: incorporating strategies beyond in silico fragmentation. Journal of Cheminformatics, 2016, 8:3.
}
\author{
  Christoph Ruttkies (\email{cruttkie@ipb-halle.de})
}
\examples{
#create settings sample object
settingsObject<-create.settings.sample()
#run MetFrag
scored.candidates<-run.metfrag(settingsObject)
#scored.candidates is a data.frame with scores and candidate properties
}