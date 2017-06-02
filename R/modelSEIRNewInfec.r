#
# Defining the seir class in R
# this defines necessary parameters and variables 
#
setClass("seirNewInfec",
	# names of parameters and variables 
   representation(sigma="numeric", # average latent period
		gamma="numeric", # average infectious period
		mu="numeric", # birth and death rate	
		seed="numeric", # number 'seed'		
		S="numeric",E="numeric",I="numeric",R="numeric",N="numeric",	# value of variables 
  	               nbCONTACT0="numeric", 
		nbCONTACT1="numeric", 
		phiPHASE="numeric", 
		probVISITER="numeric",
		probINFECTER="numeric", 
		duration="numeric",# simulation time
		nbVilles="numeric", #number of population in a metapopulation
		unitTIME="numeric", #unit of time
		periDISE="numeric", # disease period 
		typeSIMU="character", # type of simulation, 'deterministic' or 'stochastic'
		method="character",  # method of simulation, 'direct' or 'adaptivetau'
		typeRNG="character", # type of the random number generator 0: generator of C++ (good), 1: generator of Yann (fast)
		statSTATE="logical", #stationary distribution
		localExtPOP="list",
		disRecolPOP="list",
		pop="list", # containing (seir) data.frames of all populations over time	
		persistence="data.frame" # global persistence in a metapopulation			
		),

	# initial values of parameters and variables
  prototype(sigma=1/7,gamma=1/7,mu=1/(70*365),seed=23456,S=NULL,E=NULL,I=NULL,R=NULL,N=1e7,
		nbCONTACT0=300, nbCONTACT1=0.1,phiPHASE=c(0),
		probVISITER=0.01, probINFECTER=0.01, duration=1000, nbVilles=2, unitTIME=1,periDISE=365,
		typeRNG="good",typeSIMU="stoch",method="direct",statSTATE=FALSE,pop=NULL)
)

###################################################
# additional fucntion for other main fucntions
#
####
# checking a number is integer or no
#
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}
####
# sequence of time due to simulation time
#
timeseqNewInfec <- function(duration=10*365, unitTIME=1){
	return(seq(0,duration,le=duration/unitTIME))
}
#####
##pause()
pause <- function() readline("press Enter to continue!")
#####
# function 'pause'
#
pauseNewInfec <- function() readline("press Enter to continue!")
#####
# resize a vector with length given
# exemple:
# V=c(1,10),n=6, so Output:  1  1  1 10 10 10
# V=c(1,10),n=5, so Output:  1 10  1 10  1
# V=c(1,10),n=3, so Output:  1 10  1 
resizeVectorNewInfec <- function(V=c(1,10),n=6)
{
	if(is.null(V)) return(V)
	else{
		legV <- length(V)
		size <- ceiling(n/legV)
		if(n > legV){
			res<-c()
			if((n %% 2)==0){			
				for(i in 1:legV) res<-cbind(res,rep(V[i],size))
				res <- res[1:n]		
			}
			else{
				res<-rep(V,n)[1:n]		
			}
			return(res)
		}
		else
			if(n<legV) return(V[1:n])
		else return(V)
	}
}
###################################################
#
# Initializing the values of al variables for one METAPOPULATION
# according to the input or the equilibrium values after 100 years we do a deterministic simulation
#
initVarSEIRNNewInfec <- function(nbVilles=1,S=NULL,E=NULL,I=NULL,R=NULL,N=NULL,
				mu=1/(70*365),sigma=1/8,gamma=1/5,periDISE=365,phiPHASE=c(0),
				nbCONTACT0=300, nbCONTACT1=0.1, probINFECTER=0.01){
	#resize vectors of variables	
	S<-resizeVectorNewInfec(V=S,n=nbVilles)
	E<-resizeVectorNewInfec(V=E,n=nbVilles)
	I<-resizeVectorNewInfec(V=I,n=nbVilles)
	R<-resizeVectorNewInfec(V=R,n=nbVilles)
	N<-resizeVectorNewInfec(V=N,n=nbVilles)
	#parameter
	mu<-resizeVectorNewInfec(V=mu,n=nbVilles)
	nbCONTACT0<-resizeVectorNewInfec(V=nbCONTACT0,n=nbVilles)
	nbCONTACT1<-resizeVectorNewInfec(V=nbCONTACT1,n=nbVilles)
	sigma<-resizeVectorNewInfec(V=sigma,n=nbVilles)
	gamma<-resizeVectorNewInfec(V=gamma,n=nbVilles)
	phiPHASE<-resizeVectorNewInfec(V=phiPHASE,n=nbVilles)
	
	
	vecVar<-c(S=S,E=E,I=I,R=R,N=N) 
	nbVar <- length(vecVar)/nbVilles
	
	# checking values of variables
	if(nbVar==0) stop("There are not any initial variables!")
	else if(length(subset(vecVar,vecVar<0))>0) stop("Error, there are negative values!")
	else if(nbVar==5){# case 1: there are enough values of variables
		#validN = TRUE
		for(i in 1:nbVilles){
			if(N[i]!=(S[i]+E[i]+I[i]+R[i])){
			     print(paste("Values of initial variables should verify N = S + E + I + R! in city=",i))
			     N[i]<- S[i]+E[i]+I[i]+R[i]
			} 
		}
		return(list(S=S,E=E,I=I,R=R,N=N))	
	}
	else if(nbVar==4){# case 2: there are 4 values of variables
		if(is.integer0(which(vecVar==S))){ # missing  S
			#print("missing  S")
			S<-N-E-I-R	
			return(list(S=S,E=E,I=I,R=R,N=N))	
		}
		else if(is.integer0(which(vecVar==E))){# missing  E
			#print("missing  E")
			E<-N-S-I-R
			return(list(S=S,E=E,I=I,R=R,N=N))	
		}
		else if(is.integer0(which(vecVar==I))){# missing  I
			#print("missing  I")
			I<-N-S-E-R
			return(list(S=S,E=E,I=I,R=R,N=N))	
		}
		else if(is.integer0(which(vecVar==R))){# missing  R
			#print("missing  R")
			R<-N-S-I-E	
			return(list(S=S,E=E,I=I,R=R,N=N))	
		}
		else if(is.integer0(which(vecVar==N))){# missing  N
			#print("missing  N")
			N<-S+E+I+R
			return(list(S=S,E=E,I=I,R=R,N=N))	
		}
	}
	else if(nbVar<4){# case 3: there are less than 4 values of variables
		if((nbVar==1)&(!is.null(N))){ # there is N (population size for population)
		# here, we do a deterministic simulation in 100 years,
		# then get initial valeus for variables at the last time 
			S<-c(); E<-c(); I<-c(); R<-c()
			for(i in 1:nbVilles){				
				valEqual <- equiNewInfec(duration=100*365,N=N[i],mu=mu[i],nbCONTACT0=nbCONTACT0[i],probINFECTER=probINFECTER,sigma=sigma[i],gamma=gamma[i],phiPHASE=phiPHASE[i],periDISE=periDISE)

				Sdet <- valEqual[1]; S = c(S,round(Sdet));
				Edet <- valEqual[2]; E = c(E,round(Edet));
				Idet <- valEqual[3]; I = c(I,round(Idet));
				R <- c(R,N[i] - round(Sdet) - round(Edet) - round(Idet)); 

			}				
			#print(vecVar)
			return(list(S=S,E=E,I=I,R=R,N=N))		
		}
		else # there is no N, there is error
			stop("Number of initial variables should be more than 4!")
	}

}
###################
# Basic function
# function does simulation 'determionistic' or 'stochastic'
# according to the parameter 'typeSIMU'
# 
globSEIRNewInfec<-function(typeSIMU="stoch",duration=5*365,method="direct",unitTIME=1,mu=1/(70*365),
			nbCONTACT0=100,nbCONTACT1=0.0,probINFECTER=0.01,probVISITER=0.01,			
			sigma=1/8,gamma=1/5,
			periDISE=365,phiPHASE=c(0),nbVilles=1,seed=as.numeric(Sys.time()),typeRNG="good",
			S=NULL,E=NULL,I=NULL,R=NULL,N=1e5,...){
	# stochastic simulation
	if(!is.integer0(grep(typeSIMU,"stochastic")))
		stoSEIRNewInfec(method=method,duration=duration,unitTIME=unitTIME,mu=mu,
		nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1, probINFECTER=probINFECTER,probVISITER=probVISITER,
		sigma=sigma,gamma=gamma,periDISE=periDISE,phiPHASE=phiPHASE,nbVilles=nbVilles,
		seed=seed,typeRNG=typeRNG,S=S,E=E,I=I,R=R,N=N)
	# deterministic simulation
	else if(!is.integer0(grep(typeSIMU,"deterministic")))
		detSEIRNewInfec(nbVilles=nbVilles,duration=duration,unitTIME=unitTIME,S=S,E=E,I=I,R=R,N=N,periDISE=periDISE,
			mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,probINFECTER=probINFECTER,sigma=sigma,gamma=gamma,phiPHASE=phiPHASE)
	# error
	else
		stop("invalid value of typeSIMU!")
}
#######
#
# small function does deterministic simulation, that is called in the main function 'seir' above
# seir.det() produit un object de seir/sir en déterminist
detSEIRNewInfec<-function(nbVilles=1,duration=10*365,unitTIME=1,S=NULL,E=NULL,I=NULL,R=NULL,N=1e7,periDISE=365,
			mu=1/(70*365),nbCONTACT0=300,nbCONTACT1=.1,probINFECTER=0.01,sigma=1/8,gamma=1/5,phiPHASE=c(0),...){
	# checking the parametres
	sigma<-resizeVectorNewInfec(sigma,nbVilles); #sigma
	for(i in 1:length(sigma))
		if(sigma[i]<0) stop("Error, invalid value of sigma, sigma should be more than zero.")

	gamma<-resizeVectorNewInfec(gamma,nbVilles); #gamma
	for(i in 1:length(gamma))
		if(gamma[i]<0) stop("Error, invalid value of gamma, gamma should be more than zero.")

	mu<-resizeVectorNewInfec(mu,nbVilles); #mu
	for(i in 1:length(mu))
		if(mu[i]<0) stop("Error, invalid value of mu, mu should be more than zero.")

	nbCONTACT0<-resizeVectorNewInfec(nbCONTACT0,nbVilles); #nbCONTACT0
	for(i in 1:length(nbCONTACT0))
		if(nbCONTACT0[i]<0) stop("Error, invalid value of nbCONTACT0, nbCONTACT0 should be more than zero.")

	nbCONTACT1<-resizeVectorNewInfec(nbCONTACT1,nbVilles); #nbCONTACT1
	for(i in 1:length(nbCONTACT1))
		if(nbCONTACT1[i]<0) stop("Error, invalid value of nbCONTACT1, nbCONTACT1 should be more than zero.")

	# phases	
	phiPHASE<-resizeVectorNewInfec(phiPHASE,nbVilles); #phiPHASE
	# initial values for all variables
	initVar<-initVarSEIRNNewInfec(nbVilles=nbVilles,S=S,E=E,I=I,R=R,N=N,mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,probINFECTER=probINFECTER,
				sigma=sigma,gamma=gamma,periDISE=periDISE,phiPHASE=phiPHASE)
	#S
	S<-initVar$S
	#E
	E<-initVar$E
	#I
	I<-initVar$I
	#R
	R<-initVar$R
	#N
	N<-initVar$N

	times <- timeseqNewInfec(duration,unitTIME)	
	# Solving the system of differential equations 
	# (with deSolve package and code in integration between C++ and R)
	populations<-vector("list",nbVilles)
	for(i in 1:nbVilles){# simulation in the independant form
		#parameters
		pars<-c(nbCONTACT0=nbCONTACT0[i],nbCONTACT1=nbCONTACT1[i],probINFECTER=probINFECTER,periDISE=periDISE,phiPHASE=phiPHASE[i],mu=mu[i],sigma=sigma[i],gamma=gamma[i],N=N[i])

#		print(pars)
		#variables
		yini<-c(t=times[1],S=S[i],E=E[i],I=I[i])
		# calling the function 'ode' in the package 'deSolve' 
		require("deSolve")
		# and at the same time, the function 'derivs' in the code file C++
		pops.det<-as.data.frame(ode(func = "derivs", y = yini, parms = pars, times = times, dllname = "dizzysNewInfec", initfunc = "initmod"))
		valueN<-rep(N[i],nrow(pops.det))
		# Gathering the result
		valueR<-valueN-pops.det[,3]-pops.det[,4]-pops.det[,5]
		pops.det<-data.frame(pops.det[,1], pops.det[,3],pops.det[,4],pops.det[,5],valueR,valueN)
		# adding the names of variables into the result
		names(pops.det)<-c("time","S","E","P","R","N")		
		names(populations)<-paste("pop",i)
		populations[[i]]<-pops.det
	}
	#return a seir object
	return(new("seirNewInfec",nbVilles=nbVilles,duration=duration,unitTIME=unitTIME,
				N=N,S=S,E=E,I=I,R=R,periDISE=periDISE,mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,probINFECTER=probINFECTER,
					sigma=sigma,gamma=gamma,phiPHASE=phiPHASE,pop=populations,typeSIMU="deterministic"))
}
##############
#R on the server can't work with C++11 (this is a new library in C++ from 2011)
#This function uses the uniform distribution of R to generate random numbers. (trainsition matrix)
#After, we find the stationary distribution
#

#transition matrix
transMATRIX<-function(nbVilles,lower, upper)
{
	lower<-0.0
	upper<-1/nbVilles
	for (i in 1:nbVilles) {
		vecrow<-seq(nbVilles,0.0)
		n<-(nbVilles-1)
		vecGENEtp<-runif(n,min=lower,max=upper)
		valVIL=1.0-sum(vecGENEtp)
		if(i==1) vecGENE<-c(valVIL,vecGENEtp)
		else 
		if(i<nbVilles)
			vecGENE<-c(vecGENEtp[1:(i-1)],valVIL,vecGENEtp[i:n])
		else
		if(i==nbVilles)
			vecGENE<-c(vecGENEtp,valVIL)
	
		if(i==1){
#			vecGENE<-as.vector(vecGENE)
			resMATRIX<-as.data.frame(t(vecGENE))
		}
		else{
			resMATRIX<-rbind(resMATRIX,vecGENE)
		}
    }

    return (as.matrix(resMATRIX));
}
# stationary distribution
# multifly 100times
statDIS<-function(transMATRIX=matrix()){
	x<-transMATRIX
	y<-x
	for(i in 1:100)
		y<-y%*%x
	return(y)
}
###########
#
# function does one stochastic simulation,
# result is a stochastic seir object
# parameter: 'typeRNG' has two names, "good" or "fast", 'good' is the random number generator of C++, 'fast' is this of Yann
# parameter: 'method' also has two names, "direct" or "adaptivetau"
# all parameters left are vectors,
# sigma may be 'infinity' if we want SIR simulation
stoSEIRNewInfec<-function(sigma=1/8,gamma=1/5,mu=1/(70*365),seed=as.numeric(Sys.time()),S=NULL,E=NULL,I=NULL,R=NULL,N=1e5,
		nbCONTACT0=300, nbCONTACT1=0.1,phiPHASE=c(0),
		probVISITER=0.01, probINFECTER=0.1, duration=1000, nbVilles=2, unitTIME=1,periDISE=365,
		typeRNG="good",typeSIMU="stoch",method="direct",statSTATE=FALSE,...){
	#finding the stationary distribution
	initS<-c(); initE<-c(); initI<-c(); initR<-c(); initN<-c()	
	if(statSTATE){
		#statinary matrix
		matrix01<-transMATRIX(nbVilles,0,1)
		statMatriw<-statDIS(matrix01)
		statline<-statMatriw[1,]
		# find the initial values of variables for a single city
		# initial value of variables for all populations	
		initVar<-initVarSEIRNNewInfec(nbVilles=1,S=S,E=E,I=I,R=R,N=N,mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,sigma=sigma,gamma=gamma,periDISE=periDISE,phiPHASE=phiPHASE)
		initS <- initVar$S * statline
		initE <- initVar$E * statline
		initI <- initVar$I * statline
		initR <- initVar$R * statline
		initN <- initVar$N * statline
	}
	else{
	# initial value of variables for all populations	
	initVar<-initVarSEIRNNewInfec(nbVilles=nbVilles,S=S,E=E,I=I,R=R,N=N,mu=mu,nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,sigma=sigma,gamma=gamma,periDISE=periDISE,phiPHASE=phiPHASE)
	initS <- initVar$S
	initE <- initVar$E
	initI <- initVar$I
	initR <- initVar$R
	initN <- initVar$N		
	}

	# stochastic simulation of n populations
	# for the DIRECT algorithm
	if(!is.integer0(grep(method,"direct"))){
		#
		listDisRecolPOP <- vector("list",nbVilles)
		listLocalEXTPOP <- vector("list",nbVilles)
		#
		if(!is.integer0(grep(typeRNG,"good"))) rtypeRNG=0
		else
		  if(!is.integer0(grep(typeRNG,"fast"))) rtypeRNG=1
		else
			stop("invalid value of typeRNG, it should be 'good' or 'fast'")
		
		#gathering the result
		# here, we call the function 'getValeurPOPS' from C++
		resMetaPOP <- .Call("getStoSEIRNewInfec",sigma,gamma,mu,seed,initS,initE,initI,initR,
				nbCONTACT0,nbCONTACT1,phiPHASE,probVISITER,probINFECTER,
				duration,nbVilles,unitTIME,rtypeRNG,periDISE)

		populations <- vector("list",nbVilles)
		names(populations)<-paste("pop",c(1:nbVilles))
		nbElement<-nbVilles+2

		for(i in 1:nbElement){
			if(i==(nbVilles+1)) {				
				resLocalExt<- resMetaPOP[[i]]				
			}
			else if(i==(nbVilles+2)) {
				resRecol<-resMetaPOP[[i]]
			}
			else {
				populations[[i]] <- data.frame(resMetaPOP[[i]])
				names(populations[[i]]) <- c("time","S","E","P","R","N","I")
			}
		}
		#calculating the difference between the moment of the beginning of the local extinction and of the end of the local extinction$
		#
		listEveLocalExtPOP <- vector("list",nbVilles)
		listEveRecolPOP <- vector("list",nbVilles)
		#
		for( ivil in 1:nbVilles){
		#print(paste("vil==",ivil))
		#print(resLocalExt[[ivil]])
		#print(resRecol[[ivil]])
			nbLocalExt<-length(resLocalExt[[ivil]])
			nbRecol<-length(resRecol[[ivil]])
			if(nbLocalExt==1){# in the case: there is one single extinction event 
				listEveLocalExtPOP[[ivil]]<-resLocalExt[[ivil]][1]-resRecol[[ivil]][1]
				listEveRecolPOP[[ivil]]<-NULL
			}
			else if((nbLocalExt==nbRecol) & (nbLocalExt>1)){
		#there are many local extinction event and one global extinction
		listEveLocalExtPOP[[ivil]]<-unlist(resLocalExt[[ivil]]-resRecol[[ivil]])
		listEveRecolPOP[[ivil]]<-unlist(resRecol[[ivil]][2:nbLocalExt]-resLocalExt[[ivil]][1:(nbLocalExt-1)])
			}
			else if((nbLocalExt < nbRecol) & (nbLocalExt>1)){
			#there are many local extinction event, but no global extinction
		listEveLocalExtPOP[[ivil]]<-unlist(resLocalExt[[ivil]][1:nbLocalExt]-resRecol[[ivil]][1:nbLocalExt])
		listEveRecolPOP[[ivil]]<-unlist(resRecol[[ivil]][2:nbRecol]-resLocalExt[[ivil]][1:nbLocalExt])
			}
			else {
			# there are no local extinction and no global extinction
		listEveLocalExtPOP[[ivil]]<-NULL
		listEveRecolPOP[[ivil]]<-NULL
			}
		
		
		}	
	}
	else
	# for the ADAPTIVETAU algorithm
	if(!is.integer0(grep(method,"adaptivetau"))){
			
		# stochastic simulation for SEIR model
		#calling the function 'sses.AdaptiveTau' from C++
		populations <- sses.adaptivetau(sigma=sigma,gamma=gamma,mu=mu,S=S,E=E,I=I,R=R,N=N,
		nbCONTACT0=nbCONTACT0, nbCONTACT1=nbCONTACT1,phiPHASE=phiPHASE,
		probVISITER=probVISITER, probINFECTER=probINFECTER, duration=duration, 			nbVilles=nbVilles,periDISE=periDISE,
		typeSIMU=typeSIMU,method=method,statSTATE=statSTATE)
		#
		listDisRecolPOP <- vector("list",nbVilles)
		listLocalEXTPOP <- vector("list",nbVilles)
		#
		for( ivil in 1:nbVilles){
			listLocalEXTPOP[[ivil]]<-NULL
			listDisRecolPOP[[ivil]]<-NULL
		}

	}
	else
	stop("Object method should be 'direct' or 'adaptivetau'!")	
	
	
	#return an object_stoch
	return(new("seirNewInfec",sigma=sigma,gamma=gamma,mu=mu,
		seed=seed,S=initS,E=initE,I=initI,R=initR,N=initN,
		nbCONTACT0=nbCONTACT0, nbCONTACT1=nbCONTACT1,
		phiPHASE=phiPHASE,
		probVISITER=probVISITER, 
		probINFECTER=probINFECTER,
		duration=duration, nbVilles=nbVilles, 			    unitTIME=unitTIME,periDISE=periDISE,		typeRNG=typeRNG,typeSIMU="stochastic",method=method,statSTATE=statSTATE,localExtPOP=listEveLocalExtPOP,disRecolPOP=listEveRecolPOP,
		pop=populations))
#the end
}

###############
################################################
#
# Here, we are interested in the global disease persistence in a metapopulation
# We will arrange persistence time of each population in an increase.
# We save the name of population and its disease persistence time
# 
# Input: a seir object
#
setGeneric("persNewInfec", function(object,...) standardGeneric("persNewInfec"))
setMethod("persNewInfec", "seirNewInfec",
function(object,...){
	# number of population
	nbVilles <- object@nbVilles
	# population in the metapopulation
	populations <- object@pop
	#time of simulation
	dura<-object@duration
	# calculating the persistence time of each population
	persis <- sort(sapply(populations,function(x) tail(subset(x,x$I>0),1)[1,1]))
	nbresVil<-nbVilles
	vecresVil<-c()
	nbVildie<-c()
	reVecTime<-c()
	for(i in 1:nbVilles)
	{
		 if(persis[i]<dura){
			reVecTime<-c(reVecTime,persis[i])
			vecresVil<-c(vecresVil,nbresVil)
			nbVildie<-c(nbVildie,1)
			nbresVil<-nbresVil-1
		}
	}
	if(is.vector(reVecTime)){
		pers <- cbind(time=reVecTime)
		allVilles <- dimnames(pers)[[1]]
		allVilles <- as.numeric(lapply(strsplit(allVilles," "), function(x) x[2]))
		pers <- cbind(time=reVecTime,resVil=vecresVil,ndie=nbVildie,ville=allVilles)
		object@persistence <-as.data.frame(pers)
		
	}
	else
		object@persistence <-data.frame(time=c(),resVil=c(),ndie=c(),ville=c())
	# returning the result
	return(object)
}
)
##################
#
# plotting trajectory of global persistence time in a metapopulation
# There are two types we can plot
# type I: Kaplan–Meier curve of the disease persistence time
# type II: finding the entier extinction position of each population 
#
setGeneric("plot.persNewInfec", function(object,...) standardGeneric("plot.persNewInfec"))
setMethod("plot.persNewInfec", "seirNewInfec",
function(object,x="time",y="P",type="s",col="red",xlim=c(),ylim=c(),curvetype="KM",vilabline=c(),add=F,
		xlab="time",ylab="#populations non extinct",unitTIME=1,...){
	# gathering the value of persistence in a metapopulation
	pers <- object@persistence
	#time according to unit of time
	pers$time<-pers$time/unitTIME
	leng=nrow(pers)
	if(leng>0){# if there is extinction
		#number of population
		nbVilles <- object@nbVilles
		#time of simulation
		duration <- object@duration
		# TYPE I: Kaplan–Meier curve
		if(!is.integer0(grep(curvetype,"KM"))){
			if(add==F) # creating a new plot
			plot(pers$time,pers$resVil,type=type,col=col,xlab=xlab,ylab=ylab,...)
			else	# adding a new Kaplan–Meier curve on a plot given
				lines(pers$time,pers$resVil,type=type,col=col,...)
		}
		else # TYPE I: population curve
		if(!is.integer0(grep(curvetype,"population"))){
			# plotting all fluctuations of the number of infected for all populations
			plot.seirNewInfec(object,y=y,type=type,col=col,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
			pop<-c(1:nbVilles)
			allVillesEX <- pers$ville
			nbVilExt<-length(allVillesEX)
			nbvilabline <- length(vilabline)	
			if(is.null(pop)){# 'pop', population that we want plot
				for(i in 1:nbvilabline)
					for(j in 1:nbVilExt) 
						if(vilabline[i]==allVillesEX[j])
							abline(v=pers$time[j])	
			}
			else{
				for(i in 1:nbvilabline){# finding position of extinction
					flagEX_i <- FALSE
					if(!is.integer0(which(pop==vilabline[i]))){
						for(j in 1:nbVilExt)
							if(vilabline[i]==allVillesEX[j]){
								abline(v=pers$time[j])	
								flagEX_i=TRUE
							}							
					}
					if(!flagEX_i) print(paste("city ",vilabline[i],"has no local extinction"))
				}
			}
		}
		else
			stop("Error, we have two curve types, curvetype should be 'KM' or 'population' !")
	}
	else  #if there is no extinction
		stop("Error, there is no extinction!")
}
)
###################################################################################
# Estimating the survival parameters.
# calculing Survival Probability
# using the package 'survival'
#
####
# Kaplan-Meier estimate with 95% confidence bounds
setGeneric("surv.probNewInfec", function(object,...) standardGeneric("surv.probNewInfec"))
setMethod("surv.probNewInfec", "seirNewInfec",
   function(object,...){	
	pers<-object@persistence
 	my.fit <- survfit(Surv(pers[,1],pers[,3])~1)
	return(my.fit)
  }
)
# plotting  estimated Kaplan-Meier curve with 95% confidence bounds
setGeneric("plot.surv.probNewInfec", function(object,...) standardGeneric("plot.surv.probNewInfec"))
setMethod("plot.surv.probNewInfec", "seirNewInfec",
   function(object,...){
	my.fit <- surv.prob(object)
	survival:::plot.survfit(my.fit, main="Kaplan-Meier estimate with 95% confidence bounds",xlab="time", ylab="survival probability")
}
)
# estimating the global persistence rate
# by basing on Parametric Regression Models
# 
# INPUT: a seir object that has persistence
# OUTPUT: estimated persistence rate
pers.rate.objNewInfec<-function(object){
	perobj<-persistence(object)@persistence
	leng<-nrow(perobj)
	# time and ndie
	time <- perobj$time	
	ndie<- perobj$ndie
	para.obj<-survreg(Surv(time,ndie)~1, dist="exp")
	return(para.obj)
}
#####
# Confidence interval of estimated persistence rate
# INPUT:
#	pers.rate.obj : is a object of the estimated persistence rate
#	level: level of confidence

# OUTPUT:
#	confidence interval of the estimated persistence rate		
confint.pers.rateNewInfec<-function(pers.rate.obj, level=0.95){
	if(is.null(pers.rate.obj)){
		print(pers.rate.obj)
		stop("See the persistence, no local extinction in the metapopulation!")		
	}
	else{
		return(confint(pers.rate.obj,level=level))

	}
}

#############################################################################################

