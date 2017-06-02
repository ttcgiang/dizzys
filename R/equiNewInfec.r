#
# Created by TRAN Thi Cam Giang 
# Function permit us to transmit the values of parameters, following the formule of YANN \beta = -K *log(1-c)
# Then, call the function C++ in R, to caculating the equilibrium variables of the SEIR model
#
equiNewInfec<- function(duration=100*365,unitTIME=1,N=10e6,mu=1/(70*365),
    nbCONTACT0=300,nbCONTACT1=.1,probINFECTER=0.1,sigma=1/7,gamma=1/7,phiPHASE=c(0),periDISE=365)
{

	#Call the function C++ in R and at the same time, transmit the values of parameters
	# return S*, E*, I*
	sei_eq <- .Call("getEquiNewInfec",as.numeric(mu),as.numeric(nbCONTACT0),as.numeric(probINFECTER),as.numeric(sigma),as.numeric(gamma))	
	#S*, E*, I* (*N)
	sei_eq <- sei_eq*N
	#print(sei_eq)
	#Doing one deterministic simulation in 100 years
	obj_det <- detSEIRNewInfec(duration=duration,unitTIME=unitTIME,N=N,S=sei_eq[1],E=sei_eq[2],I=sei_eq[3],periDISE=periDISE,mu=mu,
			nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,probINFECTER=probINFECTER,sigma=sigma,gamma=gamma,phiPHASE=phiPHASE[1])
	
	#getting the values of variables of the population 1
	output <- obj_det@pop[[1]]
        #returning a vector of the equilibrium values of S, E, I
	return(unlist(tail(output,1)[-1]))
}

#################The end #######################

