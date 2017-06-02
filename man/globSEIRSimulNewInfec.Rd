\name{globSEIRSimulNewInfec}

\alias{globSEIRSimulNewInfec}

\title{Redoing or Continuing a Simulation.}

\description{Redoing or continuing a simulation by using values of parameters given.}

\usage{
globSEIRSimulNewInfec(object,typeSIMU="stoch",continue=F,
	duration=5*365,method="direct",unitTIME=1,mu=1/(70*365),
	nbCONTACT0=200,nbCONTACT1=.1,probINFECTER=0.01,probVISITER=0.1,
	sigma=1/8,gamma=1/5,
	periDISE=365,phiPHASE=c(0),nbVilles=1,seed=as.numeric(Sys.time()),
	typeRNG="good", append=TRUE, t0=NULL,
	S=NULL,E=NULL,I=NULL,R=NULL,N=1e7)
}

\arguments{
  \item{object}{a seir object. If missing this object, we create a new seir object by using the values of state variables and of parameters given.}
  \item{typeSIMU}{type of simulation, 'deterministic' or 'stochastic'.}
  \item{continue}{logical (T/F). If continue=F, it means to redo silulation. Else, if continue=T, it means to continue to do silulation.}
  \item{duration}{time values over which to perform the numerical integration.}
  \item{method}{simulation algorithm,'direct' or 'adaptivetau'. If we chose method="direct", it means that we implement the direct algorithm of Gillespie in 1977. In contrast, if we chose method="adaptivetau", it means that we implement the adaptive tau-leaping approximation for simulating the trajectory of a continuous-time Markov process.}
  \item{unitTIME}{unit of time.}
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
  \item{seed}{random seed for random number generator.}
  \item{typeRNG}{random number generator, 'good' or 'fast'. If we chose typeRNG="good", it means that we use the random number generator in C++. In the other side, we choose typeRNG="fast", it means that we use the random number generator of the professor Yann Chevaleyre at university of Paris 13. With the random generator's Yann Chevaleyre, we find that it is faster than the random generator in C++, in contrast, its accuracy is less than that of the random generator's Yann Chevaleyre.}
  \item{S}{initial value of the state variable 'S', this is the number of susceptible individuals at the time t0.}
  \item{E}{initial value of the state variable 'E', this is the number of exposed individuals at the time t0.}
  \item{I}{initial value of the state variable 'I', this is the number of infected individuals at the time t0.}
  \item{R}{initial value of the state variable 'R', this is the number of recovered individuals in immunity at the time t0.}
  \item{N}{number of population.}
  \item{phiPHASE}{phase of the contact rate (in radians).}
  \item{append}{logical(T/F). It is available when continue=T.
	If continue=T and append=T, it means that we want our new object contains the old data + the new.
        If continue=T and append=F, we only get the new.}
  \item{t0}{numeric. It means that if we want the start time after what has already been simulated ("t0 = NULL") or if we prefer that time restarts any value (for example "t0 = 1").}
}

\details{For this fucntion, to get the values of arguments, first, for the arguments missed, we get its values from the seir object, in contrast, we get its given values. Moreover, in the set of arguments, there are some arguments we can give each subpopulation each value's.
  \bold{S, E, I, R, N}: the initial values of state variables for 'nbVilles' subpopulations. Multiple values can be specified so that each subpopulation can be given its own state variables. If there are fewer values than subpopulations they are recycled in the standard fashion. Subpopulations will all be simulated in the first value specified.
  \bold{mu, nbCONTACT0, nbCONTACT1, probVISITER, probINFECTER, sigma, gamma, phiPHASE}: the parameters of simulation for 'nbVilles' subpopulations. Multiple values can be specified so that each subpopulation can be given its own parameters. If there are fewer values than subpopulations they are recycled in the standard fashion. Subpopulations will all be simulated in the first value specified.
}

\value{
	Result returned is an object of the 'seir' class.
}

\author{TRAN Thi Cam Giang}

\seealso{
	'simul' function in the 'dizzys' package.
}

\examples{
	#STO, STO
	sto<-globSEIRNewInfec(N=1e6,typeSIMU="stoch",duration=5*365,nbVilles=2)
	plot(globSEIRSimulNewInfec(sto,typeSIMU="stoch",continue=TRUE,duration=5*365,nbCONTACT1=0,phiPHASE=c(pi/2,0)),col=c(1,2))
	#DET, DET
	det<-globSEIRNewInfec(N=10e4,typeSIMU="deter",duration=50*365)
	plot(globSEIRSimulNewInfec(det,typeSIMU="deter",continue=TRUE,duration=5*365,nbCONTACT1=0.1,phiPHASE=pi))
}

\keyword{ seir model }
\keyword{ R package }% __ONLY ONE__ keyword per line
