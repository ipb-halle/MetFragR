\name{run.metfrag}
\alias{run.metfrag}
\title{Run MetFrag by specifying all parameters given with the settings object as list}
\usage{
	run.metfrag(settingsObject)
}
\arguments{
  \item{settingsObject}{A \code{list} of parameter settings}
}
\description{
The function uses defined settings like database, fragmentation and scoring
parameters to initialise and run the MetFrag process. The FragmenterScore
calculation is based on the match of in silico generated fragments of the
candidate molecules to the given tandem mass spectrum. A data.frame of candidates
molecules is returned sorted by their final score.
}
\value{
Returns a data.frame containing ranked candidates represented by different properties and identifiers.
}
\references{
Ruttkies C., Schymanski E.L. et al, MetFrag relaunched: incorporating strategies beyond in silico fragmentation. Journal of Cheminformatics, 2016, 8:3.
}
\author{
  Christoph Ruttkies (\email{cruttkie@ipb-halle.de})
}
\examples{
#define settings
settingsObject<-list()
settingsObject[["DatabaseSearchRelativeMassDeviation"]]<-5.0
settingsObject[["FragmentPeakMatchAbsoluteMassDeviation"]]<-0.001
settingsObject[["FragmentPeakMatchRelativeMassDeviation"]]<-5.0
settingsObject[["MetFragDatabaseType"]]<-"PubChem"
settingsObject[["NeutralPrecursorMass"]]<-253.966126
settingsObject[["PeakList"]]<-matrix(c(
90.97445, 681,
106.94476, 274,
110.02750, 110,
115.98965, 95,
117.98540, 384,
124.93547, 613,
124.99015, 146,
125.99793, 207,
133.95592, 777,
143.98846, 478,
144.99625, 352,
146.00410, 999,
151.94641, 962,
160.96668, 387,
163.00682, 782,
172.99055, 17,
178.95724, 678,
178.97725, 391,
180.97293, 999,
196.96778, 720,
208.96780, 999,
236.96245, 999,
254.97312, 999), ncol=2, byrow=TRUE)
settingsObject[["NeutralPrecursorMolecularFormula"]]<-"C7H5Cl2FN2O3"
settingsObject[["PrecursorCompoundIDs"]]<-c("50465", "57010914", "56974741", "88419651", "23354334")
#run MetFrag
scored.candidates<-run.metfrag(settingsObject)
#scored.candidates is a data.frame with scores and candidate properties
}
\keyword{fragmentation, MetFrag}
