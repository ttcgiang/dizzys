\name{persNewInfec}

\alias{persNewInfec}

\title{persNewInfec in a Metapopulation}

\description{Generic function for calculing number of subpopulations in a metapopulation not extincts at time t.}

\usage{persNewInfec(object)}

\arguments{
  \item{object}{a seir object.}
}

\details{The function 'persNewInfec' gives us the number of subpopulations, in a metapopualtion, in which the disease is not extinct at time t by using the parameter 'object' given.}

\value{'persNewInfec' returns a seir object with the value of the slot 'persNewInfec'. The slot 'persNewInfec' is a data frame. This data frame has five columns, 'time' means being the time at which there are some extinct subpopulations, at time t, we have, 'nbVilles' is the number of not extinct subpopulations before the time t, 'ndie' means being the number of died subpopulations, 'dieVille' gives us which subpopulation is died, 'remain' is the number of not extinct subpopulations after the time t.}

\references{Matt J. Keeling and Pejman Rohani (2008) Modeling Infectious Diseases IN HUMANS AND ANIMALS. Princeton University Press.}

\author{TRAN Thi Cam Giang}

\examples{
	obj1<-obj1<-globSEIRSimulNewInfec(nbVilles=5,N=1e5,nbCONTACT0=100,duration=365*30)
	objper<-persNewInfec(obj1)
	objper@persistence
}

\keyword{ metapopulation }
\keyword{ persNewInfec }% __ONLY ONE__ keyword per line
