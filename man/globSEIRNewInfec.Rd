\name{globSEIRNewInfec}

\alias{globSEIRNewInfec}

\title{Creat a seir Object}

\description{This is important function, it creates a seir object for 'seirNewInfec' class.}

\usage{
globSEIRNewInfec(typeSIMU="stoch",duration=5*365,method="direct",unitTIME=1,mu=1/(70*365),
			nbCONTACT0=100,nbCONTACT1=0.0,probINFECTER=0.01,probVISITER=0.01,			
			sigma=1/8,gamma=1/5,
			periDISE=365,phiPHASE=c(0),nbVilles=1,seed=as.numeric(Sys.time()),typeRNG="good",
			S=NULL,E=NULL,I=NULL,R=NULL,N=1e5)
}

\arguments{
  \item{typeSIMU}{type of simulation, 'deterministic' or 'stochastic'.}
  \item{duration}{time values over which to perform the numerical integration.}
  \item{S}{initial value of the state variable 'S', this is the number of susceptible individuals at the time t0.}
  \item{E}{initial value of the state variable 'E', this is the number of exposed individuals at the time t0.}
  \item{I}{initial value of the state variable 'I', this is the number of infected individuals at the time t0.}
  \item{R}{initial value of the state variable 'R', this is the number of recovered individuals in immunity at the time t0.}
  \item{N}{number of population.}
  \item{mu}{per capita birth and death rates (per day).}
  \item{nbCONTACT0}{mean value of the number of contact per day}
  \item{nbCONTACT1}{amplitude of the number of contact per day}
  \item{probINFECTER}{the probability that a susceptible individual native from i being in contact with another infected individual native from k gets infected}
  \item{probVISITER}{the probability that an individual from subpopulation i visits subpopulation j}
  \item{sigma}{transition rate from exposed (E) to infected (I), per day (inverse of the average duration of the latency period).}
  \item{gamma}{recovery rate (per day).}
  \item{periDISE}{period of the contact rate (in days).}
  \item{phiPHASE}{phase of the contact rate (in radians).}
  \item{nbVilles}{number of subpopulations in the metapopulation.}
  \item{seed}{random seed.}
  \item{typeRNG}{random number generator, 'good' or 'fast'. If we chose typeRNG="good", it means that we use the random number generator in C++. In the other side, we chose typeRNG="fast", it means that we use the random number generator of the professor Yann Chevaleyre at university of Paris 13. With the random generator's Yann Chevaleyre, we find that it is faster than the random generator in C++, in contrast, its accuracy is less than that of the random generator's Yann Chevaleyre.}
  \item{method}{simulation algorithm,'direct' or 'adaptivetau'. If we chose method="direct", it means that we implement the direct algorithm of Gillespie in 1977. In contrast, if we chose method="adaptivetau", it means that we implement the adaptive tau-leaping approximation for simulating the trajectory of a continuous-time Markov process.}
  \item{unitTIME}{unit of time for simulation.}  
  \item{\dots}{further parameters.}
}

\details{To implement this function, we have three step. The first step, in R,  is to get the initial values of state variables and the values of parameter, then, sends those values to the C++ functions. The second step is to do simulation in the C++ implementation. The third step is to send the result from the C++ function to R. We will get the simulation resultat in R.

Moreover, in the set of arguments, there are some arguments we can give each subpopulation each value's.

  \bold{S, E, I, R, N}: the initial values of state variables for 'nbVilles' subpopulations. Multiple values can be specified so that each subpopulation can be given its own state variables. If there are fewer values than subpopulations they are recycled in the standard fashion. Subpopulations will all be simulated in the first value specified.

  \bold{mu, nbCONTACT0, nbCONTACT1,probINFECTER,probVISITER, sigma, gamma, phiPHASE}: the parameters of simulation for 'nbVilles' subpopulations. Multiple values can be specified so that each subpopulation can be given its own parameters. If there are fewer values than subpopulations they are recycled in the standard fashion. Subpopulations will all be simulated in the first value specified.
}

\value{Result returned is a seir object of the 'seirNewInfec' class with all values for all slots.}

\references{Norman Matloff (2009), The Art of R Programming.}

\author{TRAN Thi Cam Giang}

\seealso{
	'globSEIRNewInfec', 'detSEIRNewInfec', 'stoSEIRNewInfec' fuctions in 'dizzysNewInfec' package.
}

\examples{
	obj<-globSEIRNewInfec()
	plot(obj)
	obj<-globSEIRNewInfec(typeSIMU="sto",nbVilles=3,N=c(1e7,1e6))
	plot(obj)
	obj<-globSEIRNewInfec(typeSIMU="deter",nbVilles=3,N=c(1e7,1e6))
	plot(obj)	
}

\keyword{ deterministic model }
\keyword{ stochastic model }% __ONLY ONE__ keyword per line
