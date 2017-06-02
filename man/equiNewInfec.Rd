\name{equiNewInfec}

\alias{equiNewInfec}

\title{Finding Endemic Equilibrium of a seasonally-forced SEIR/SIR model}

\description{'equiNewInfec' returns the values of the state variables of an SEIR/SIR model at endemic equilibrium.}

\usage{equiNewInfec(duration=100*365,unitTIME=1,N=10e6,mu=1/(70*365),
    nbCONTACT0=300,nbCONTACT1=.1,probINFECTER=0.1,sigma=1/7,gamma=1/7,phiPHASE=c(0),periDISE=365)}

\arguments{
  \item{duration}{time values over which to perform the numerical integration.}
  \item{unitTIME}{time unit of simulation.}
  \item{N}{number of population.}
  \item{mu}{per capita birth and death rates (per day).}
  \item{nbCONTACT0}{mean value of the number of contact per day}
  \item{nbCONTACT1}{amplitude of the number of contact per day}
  \item{probINFECTER}{the probability that a susceptible individual native from i being in contact with another infected individual native from k gets infected}
  \item{sigma}{transition rate from exposed (E) to infected (I), per day (inverse of the average duration of the latency period.}
  \item{gamma}{recovery rate (per day).}
  \item{phiPHASE}{phase of the contact rate (in radians).}
  \item{periDISE}{period of the contact rate (in days).}
}
\details{The host population is supposed to be at demographic equilibrium with the births balancing the deaths (births and death occur with the same rate 'mu') and a total population size constant and equal to 1. The contact rate is forced by a sinusoid of period 'periDISE' and phase 'phiPHASE'. The SEIR epidemiological model is defined by the following set of differential equations:

	beta = -nbCONTACT0*log(1-probINFECTER);

	dS/dt =  mu*N - nbCONTACT0 * (1 + nbCONTACT1 * cos(2 * pi * t / periDISE + phiPHASE))*log[1-(I*probINFECTER/N)] * S * I - mu * S

	dE/dt = nbCONTACT0 * (1 + nbCONTACT1 * cos(2 * pi * t / periDISE + phiPHASE))*log[1-(I*probINFECTER/N)] * S * I - mu * E - sigma * E

	dI/dt = sigma *  E - mu * I - gamma * I

	dR/dt = gamma * I - mu * R
	

'nbCONTACT0' is the mean value of the contact number and 'nbCONTACT1' is the amplitude (in percentage of the mean) around the mean. The endemic equilibrium value of this sinusoidally-forced epidemiological system is found in two steps. First the endemic equilibium point (S*,E*,I*,R*) of the unforced system ('nbCONTACT1 = 0') is found from the following analytical equations:

	S* = (gamma+mu)*(sigma+mu)/(beta * sigma)

	E* = mu*((1/(sigma+mu)) - ((gamma+mu)/(beta * sigma)))

	I* = mu*((beta*sigma - (gamma+mu)*(sigma+mu))/(beta *(gamma+mu)*(sigma+mu)))

	R* = 1 - S* - E* - I*

This endemic equilibium point is used as a starting point for numerical integration of the system of differential equations. This is done by the 'detSEIRNewInfec' function that uses the 'ode' function of the 'deSolve' package. The last values of the of the state variable are returned. These correspond to a given point on the equilibrium limit cycle (assuming that this equilibrium limit cycle is reached) and the value of this point will depend on the phase 'phiPHASE'.}

\value{'equiNewInfec' returns a named vector of 3 numerical proportions containing the values of the state variables S, E and I on one point of the equilibrium limit cycle (provided that this limit cycle is reached). In case this limit cycle is reached, the value of the S, E and I state variables will depend on the input phase value 'phi'.}

\author{TRAN Thi Cam Giang}

\references{Anderson RM & May RM (1991) Infectious Diseases of Humans - Dynamics and Control. Oxford University Press.}

\seealso{'seir' in package 'dizzyz' and 'ode' in package 'deSolve'.}

\examples{
## The point on the limit cycle depends on the input phase value 'phi':
	res<-equiNewInfec(duration=100*365,unitTIME=1,N=10e6,mu=1/(70*365),
    nbCONTACT0=300,nbCONTACT1=.1,probINFECTER=0.1,sigma=1/7,gamma=1/7,phiPHASE=c(0),periDISE=365)	
	print(res)
}

\keyword{SEIR/SIR model}
\keyword{limit cycle equilibrium}% __ONLY ONE__ keyword per line
