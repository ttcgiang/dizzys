#####################################################################################
# This is a method in the package "dizzys" that is called 'pop'
# This method permits to create a new data.frame (like a new population)
# by using a seirNewInfec object, position of populations chosen and a function to calculate

# Initializing the method
pop <- function (object,...) 
  UseMethod("pop")

# Defining the method
# object: one seirNewInfec object
# subset: position of populations chosen
# fct: name of function that we want to calculate
pop.seirNewInfec <- function(object,subset=c(),fct="sum",...){
	#get the value of POPULATION in the seirNewInfec object
	pops.stoch <- object@pop
	#number of POPULATION
	nbVilles <- length(pops.stoch)
	if(nbVilles == 1) return(pops.stoch) # if there is one population
	else
	if(nbVilles >1){ #if there is many populations
		if(is.null(subset)){#if no population chosen
			pops.sto.subset <- pops.stoch
		}
		else{#if have populations chosen
			pops.sto.subset <-  vector("list")
			namepop<-c()
			nbvalidvil <- 0; nsubset <-length(subset)
			#checking city at subset[i] existing or non in popsall list
			for(i in 1:nsubset){			
				ixvil <- which(c(1:nbVilles)==subset[i])
				#print(ixvil)			
				if(!is.integer0(ixvil)){
					namepop <- cbind(namepop,paste("pop",ixvil))
					nbvalidvil <- nbvalidvil+1
					pops.sto.subset[[nbvalidvil]] <- pops.stoch[[ixvil]]				
				}
				else
				print("Invalid city index in subset!")			
			}
			names(pops.sto.subset)<-namepop				
		}
		#working with fct
		if(is.null(fct)){#returning all 
			return(pops.sto.subset)
		}
		else{#returning a seul dataframe according to the parametre "fct"
			#rejecting the parametre "time" in list		
			tmp <- pops.sto.subset
			# number of populations is 
			nbpop <- length(tmp)
			# the length is 
			l <- nrow(tmp[[1]])
			# the number of variables is 
			nbvar <- ncol(tmp[[1]])
			# converting a list to array:
			tmp <- array(unlist(tmp),c(l,nbvar,nbpop))
			# apply any function "fct" for all populations:
			#print(class(apply(tmp,c(1,2),fct)))
			return(apply(tmp,c(1,2),fct))				
		}
	}
	else
		stop("invalid number of cities!")	
}
# printing data of an object of seir stochastic
########################
# Defining all operations for seir object
#adding
'+.seirNewInfec'<-function(a=object, b=object){
	cat("addition of two seirNewInfec objects.\n")
}

#subtract
'-.seirNewInfec'<-function(a=object, b=object){
	cat("subtraction of two seirNewInfec objects.\n")
}

#multiply
'*.seirNewInfec'<-function(a=object, b=object){
	cat("multiplication of two seirNewInfec objects.\n")
}

#divide
'/.seirNewInfec'<-function(a=object, b=object){
	cat("division of two seirNewInfec objects.\n")
}
#exponentiation (right to left)
'^.seirNewInfec'<-function(a=object, b=object){
	cat("exponentiation of two seirNewInfec objects.\n")
}
#exponentiation (right to left)
'**.seirNewInfec'<-function(a=object, b=object){
	cat("exponentiation of two seirNewInfec objects.\n")
}
#division
'%%.seirNewInfec'<-function(a=object, b=object){
	cat("modulus of two seirNewInfec objects.\n")
}
#division
'%/%.seirNewInfec'<-function(a=object, b=object){
	cat("integer division of two seirNewInfec objects.\n")
}
#comparing
'<.seirNewInfec'<-function(a=object, b=object){
	cat("'less than' operator of two seirNewInfec objects.\n")
}
#
'<=.seirNewInfec'<-function(a=object, b=object){
	cat("'less than or equal to' operator of two seirNewInfec objects.\n")
}
#
'>.seirNewInfec'<-function(a=object, b=object){
	cat("'more than or equal to' operator of two seirNewInfec objects.\n")
}
#
'>=.seirNewInfec'<-function(a=object, b=object){
	cat("'more than or equal to' operator of two seirNewInfec objects.\n")
}
#
'==.seirNewInfec'<-function(a=object, b=object){
	cat("'exactly equal to' operator of two seirNewInfec objects.\n")
}
#
'!=.seirNewInfec'<-function(a=object, b=object){
	cat("'not equal to' operator of two seirNewInfec objects.\n")
}
#
'!.seirNewInfec'<-function(a=object){
	if(is.null(a)) return(TRUE)
	else return(FALSE)
}
#logic and or
'|.seirNewInfec'<-function(a=object, b=object){
	if(is.null(a)&&is.null(b)) return(FALSE)
	else return(TRUE)
}
#
'&.seirNewInfec'<-function(a=object, b=object){
	if(!is.null(a)&& !is.null(b)) return(TRUE)
	else return(FALSE)
}

########################
# additional function that helps an other main function
# extracting the identical values in the vector 'param'
# and the same time, the time is increased.
# exemple: 
# Input: param=c(0.2,0.2,0.1,0.1,0.3,0.3),duration=c(5,10,8)
# Output: 
#  new_param= c(0.2, 0.1, 0.3)
#  new_duration= c(5, 15, 23)
#
paramdurNewInfec<-function(ville=1,nbVilles=2,param=c(0.2,0.2,0.1,0.1,0.3,0.3),duration=c(5,10,8)){
	nbfois<-length(param)/nbVilles
	pamvil<-c()
	for(i in 0:(nbfois-1)) pamvil<-c(pamvil,param[ville+i*nbVilles]) 
	
	newpamvil<-c()
	newdurvil<-c()
	if(length(pamvil)>1){
		for(i in 1:(length(pamvil)-1)){
			if(pamvil[i]!=pamvil[i+1]){
				newpamvil<-c(newpamvil,pamvil[i])
				newdurvil<-c(newdurvil,sum(duration[1:i]))
			}
			#else newdurvil<-c(newdurvil,sum(duration[i],duration[i+1]))

			if(i==(length(pamvil)-1)&&(pamvil[i]!=pamvil[i+1])){
				newpamvil<-c(newpamvil,pamvil[i+1])
				newdurvil<-c(newdurvil,sum(duration[1:(i+1)]))
			}
			if(i==(length(pamvil)-1)&&(pamvil[i]==pamvil[i+1])){
				newpamvil<-c(newpamvil,pamvil[i])
				 newdurvil<-c(newdurvil,sum(duration[1:(i+1)]))
			}	
		}
	
	}
	else{
		newpamvil<-pamvil
		newdurvil<-duration
	}
	return(list(newpamvil,newdurvil))
}
###########
# Defining the function 'summary' for a seirNewInfec object.
# This is a generic function used to produce result summaries of the results of various model fitting functions.
# Here, the function 'summary' of a seirNewInfec object will show all informations of parameters and initial values of variables (such as: name, signification, value of each parameter and each variable
#
summary.seirNewInfec<-function(object,...)
{
	#population over time
	pops.stoch <- object@pop
	#number of populations
	nbVilles <- object@nbVilles
	#type of simulation (deterministic or stochastic)
	cat(paste("Type:",object@typeSIMU),fill=60)
	#simulation algorithm,'direct' or 'adaptivetau'	
	cat(paste("Method:",object@method),fill=60)
	# random number generator
	cat(paste("RNG: built-in C++"),fill=60)
	# time of simulation
	cat(paste("Duration:",sum(object@duration),"days"),fill=60)
	#unit of time ('day' or 'year')
	unitTIME <- object@unitTIME
	if(unitTIME<=1) cat(paste("Time unit:",unitTIME,"day"),fill=60)
	else cat(paste("Time unit:",unitTIME,"days"),fill=60)
	cat(paste("Seed:",round(object@seed,digits=5)),fill=60)
	cat(paste("Number of populations:",object@nbVilles),fill=60)
	unityear<-365
	
	# parameters \nbCONTACT0, \nbCONTACT1, \mu, \gamma, \sigma for all populations
	for(i in 1: nbVilles){
		cat(paste("$-----pop",i),fill=60)
		cat("Coefficients:\n")
		#population size 
		cat(paste(format(paste("Population size N of population ",i,":"),width=50),format(object@N[i],width=10)),fill=60)
		#\nbCONTACT0 : mean value of the number of contact per unit time
		nbCONTACT0dura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@nbCONTACT0,object@duration)		
		vilnbCONTACT0<-nbCONTACT0dura[[1]]
		vilduration<-round(nbCONTACT0dura[[2]]/unityear,digits=2)
		nbsimulnbCONTACT0<-length(vilnbCONTACT0)
		if(nbsimulnbCONTACT0==1){
			cat(paste(format("Mean contact number nbCONTACT0:",width=50),format(round(vilnbCONTACT0,digits=5),width=10),
					"/day/ind",sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulnbCONTACT0){
			if(j==1) {
				    cat(paste(format("Mean contact rate nbCONTACT0:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(paste(round(vilnbCONTACT0[j],digits=5),"/day/ind",sep=""),width=22)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(paste(round(vilnbCONTACT0[j],digits=5),"/day/ind",sep=""),width=22)),fill=60)
			}
			}
		}	
		#\nbCONTACT1 : amplitude of the number of contact per unit time \nbCONTACT1
		nbCONTACT1dura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@nbCONTACT1,object@duration)		
		vilnbCONTACT1<-nbCONTACT1dura[[1]]
		vilduration<-round(nbCONTACT1dura[[2]]/unityear,digits=2)
		nbsimulnbCONTACT1<-length(vilnbCONTACT1)
		if(nbsimulnbCONTACT1==1){
			cat(paste(format("Amplitude of contact number nbCONTACT1:",width=50),format(round(vilnbCONTACT1,digits=5),width=10),
				sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulnbCONTACT1){
			if(j==1) {
				    cat(paste(format("Amplitude of contact rate nbCONTACT1:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(round(vilnbCONTACT1[j],digits=5),width=7)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(round(vilnbCONTACT1[j],digits=5),width=7)),fill=60)
			}
			}
		}
		#\probINFECTER : is the probability that a susceptible individual native from i being in contact with another infected individual native from k gets infected. 
		probINFECTERdura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@probINFECTER,object@duration)		
		vilprobINFECTER<-probINFECTERdura[[1]]
		vilduration<-round(probINFECTERdura[[2]]/unityear,digits=2)
		nbsimulprobINFECTER<-length(vilprobINFECTER)
		if(nbsimulprobINFECTER==1){
			cat(paste(format("Amplitude of contact number probINFECTER:",width=50),format(round(vilprobINFECTER,digits=5),width=10),
				sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulprobINFECTER){
			if(j==1) {
				    cat(paste(format("Amplitude of contact rate probINFECTER:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(round(vilprobINFECTER[j],digits=5),width=7)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(round(vilprobINFECTER[j],digits=5),width=7)),fill=60)
			}
			}
		}
		#probVISITER[i,j] : the probability that an individual from subpopulation i visits subpopulation j
		probVISITERdura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@probVISITER,object@duration)		
		vilprobVISITER<-probVISITERdura[[1]]
		vilduration<-round(probVISITERdura[[2]]/unityear,digits=2)
		nbsimulprobVISITER<-length(vilprobVISITER)
		if(nbsimulprobVISITER==1){
			cat(paste(format("Amplitude of contact number probVISITER:",width=50),format(round(vilprobVISITER,digits=5),width=10),
				sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulprobVISITER){
			if(j==1) {
				    cat(paste(format("Amplitude of contact rate probVISITER:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(round(vilprobVISITER[j],digits=5),width=7)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(round(vilprobVISITER[j],digits=5),width=7)),fill=60)
			}
			}
		}
	
		# \sigma : average latent period
		sigmadura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@sigma,object@duration)		
		vilsigma<-sigmadura[[1]]
		vilduration<-round(sigmadura[[2]]/unityear,digits=2)
		nbsimulsigma<-length(vilsigma)
		if(nbsimulsigma==1){
			cat(paste(format("Latency rate sigma:",width=50),format(round(vilsigma,digits=5),width=10),
					"/day\n(mean duration of latency period: ",1/vilsigma," days)",sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulsigma){
			if(j==1) {
				    cat(paste(format("Latency rate sigma:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(paste(round(vilsigma[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(paste(round(vilsigma[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			}
		}	
		#\gamma : average infection period
		gammadura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@gamma,object@duration)		
		vilgamma<-gammadura[[1]]
		vilduration<-round(gammadura[[2]]/unityear,digits=2)
		nbsimulgamma<-length(vilgamma)
		if(nbsimulgamma==1){
			cat(paste(format("Recovery rate gamma:",width=50),format(round(vilgamma,digits=5),width=10),
					"/day\n(Mean duration of infectious period: ",1/vilgamma," days)",sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulgamma){
			if(j==1) {
				    cat(paste(format("Recovery rate gamma:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(paste(round(vilgamma[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(paste(round(vilgamma[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			}
		}	
		#\mu : birth rate and death rate
		mudura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@mu,object@duration)		
		vilmu<-mudura[[1]]
		vilduration<-round(mudura[[2]]/unityear,digits=2)
		nbsimulmu<-length(vilmu)
		if(nbsimulmu==1){
			cat(paste(format("Birth and death rate mu:",width=50),format(round(vilmu,digits=7),width=10),
					"/day\n(Mean life expectancy: ",1/(vilmu*365)," years)",sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulmu){
			if(j==1) {
				    cat(paste(format("Birth and death rate mu:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(paste(round(vilmu[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(paste(round(vilmu[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			}
		}	
		#synchrony \phiPHASE, phase of forcing
		phiPHASEdura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@phiPHASE,object@duration)		
		vilphiPHASE<-phiPHASEdura[[1]]
		vilduration<-round(phiPHASEdura[[2]]/unityear,digits=2)
		nbsimulphiPHASE<-length(vilphiPHASE)
		if(nbsimulphiPHASE==1){
			cat(paste(format("Synchrony phiPHASE:",width=50),format(round(vilphiPHASE,digits=5),width=10),
					sep=""),fill=60)
		}
		else{
			for(j in 1:nbsimulphiPHASE){
			if(j==1) {
				    cat(paste(format("Synchrony phiPHASE:",width=35),
				    format(paste("[0,",vilduration[j],"] years",sep=""),width=15),
				    format(paste(round(vilphiPHASE[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			else {   
				 cat(paste(format("                          ",width=35),
			    format(paste("(",vilduration[j-1],",",vilduration[j],"] years",sep=""),width=15),
			    format(paste(round(vilphiPHASE[j],digits=5),"/day",sep=""),width=22)),fill=60)
			}
			}
		}	
		# State variables:
		cat("\nState variables:",fill=60)			
		popi <- pops.stoch[[i]]
		nbcol <- ncol(popi)
		popi <-popi[,c(2:nbcol)]
		print(apply(popi,2,summary))
	}	
}
#####################################################################################
# printing an object of seirNewInfec stochastic
print.seirNewInfec <- function(object,...) {
#printing the values of the populations as "time	S	E	P	R	N	I"
	object
   }

#####################################################################################
# printing coef of an object of seirNewInfec.stoch
coef.seirNewInfec<-function(object,...){
	#cat("\nPrinting the values of parameters of an object in the seirNewInfec stochastic class\n")
	cat("@initial values of variables:",fill=60)
	varis <- c(N = as.numeric(object@N),
			S = as.numeric(object@S),
			E = as.numeric(object@E),
			I = as.numeric(object@I),
			R = as.numeric(object@R)
			)
	print(varis)
	cat("\n@values of parameters :",fill=60)
	parameters <- c(duration = object@duration,
			unitTIME = object@unitTIME,
			nbVilles = object@nbVilles,			
			duration = object@duration,
			mu = object@mu,
			nbCONTACT0 = object@nbCONTACT0,
			nbCONTACT1 = object@nbCONTACT1,
			probVISITER = object@probVISITER,
			probINFECTER = object@probINFECTER,
			sigma = object@sigma,
			gamma = object@gamma,
			phiPHASE = object@phiPHASE,
			seed = object@seed)
	print(parameters)
#	cat(paste("rng:",object@rng))
  }
########################################
#
# function 'str' here is inherited from the function 'str' in the package 'utils'
# the function 'str' here permits to compactly display the structure of a seirNewInfec object. It is similaire to the fucntion 'summary'.
# displaying the beginning values of variables for all populations
# displaying the values of parameters
#
str.seirNewInfec<-function(object,...){
	cat("@populations:",fill=60)
	nbVilles <- object@nbVilles
	for(i in 1:nbVilles){	
		# beginning values of variables for all populations
		cat(paste("..$pop",i,": data frame"),fill=60)	
		cat("..$coef:",fill=60)	
		N<-object@N[i]
		S<-object@S[i]
		E<-object@E[i]
		I<-object@I[i]
		R<-object@R[i]
		# time of simulation
		duration<-sum(object@duration)
		#unit of time
		unitTIME<-object@unitTIME
		#nbCONTACT0		
		nbCONTACT0dura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@nbCONTACT0,object@duration)		
		vilnbCONTACT0<-nbCONTACT0dura[[1]]
		vilduration<-nbCONTACT0dura[[2]]
		nbsimulnbCONTACT0<-length(vilnbCONTACT0)
		if(nbsimulnbCONTACT0==1){
			nbCONTACT0<-round(vilnbCONTACT0,digits=5)
		}
		else{
			nbCONTACT0<-c()
			for(j in 1:nbsimulnbCONTACT0){
				nbCONTACT0<-c(nbCONTACT0,round(vilnbCONTACT0[j],digits=5))
			}
		}	
		#nbCONTACT1
		nbCONTACT1dura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@nbCONTACT1,object@duration)		
		vilnbCONTACT1<-nbCONTACT1dura[[1]]
		vilduration<-nbCONTACT1dura[[2]]
		nbsimulnbCONTACT1<-length(vilnbCONTACT1)
		if(nbsimulnbCONTACT1==1){
			nbCONTACT1<-round(vilnbCONTACT1,digits=5)
		}
		else{
			nbCONTACT1<-c()
			for(j in 1:nbsimulnbCONTACT1){
				nbCONTACT1<-c(nbCONTACT1,round(vilnbCONTACT1[j],digits=5))
			}
		}
		#probVISITER, 
		probVISITERdura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@probVISITER,object@duration)		
		vilprobVISITER<-probVISITERdura[[1]]
		vilduration<-probVISITERdura[[2]]
		nbsimulprobVISITER<-length(vilprobVISITER)
		if(nbsimulprobVISITER==1){
			probVISITER<-round(vilprobVISITER,digits=5)
		}
		else{
			probVISITER<-c()
			for(j in 1:nbsimulprobVISITER){
				probVISITER<-c(probVISITER,round(vilprobVISITER[j],digits=5))
			}
		}
		#probINFECTER
		probINFECTERdura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@probINFECTER,object@duration)		
		vilprobINFECTER<-probINFECTERdura[[1]]
		vilduration<-probINFECTERdura[[2]]
		nbsimulprobINFECTER<-length(vilprobINFECTER)
		if(nbsimulprobINFECTER==1){
			probINFECTER<-round(vilprobINFECTER,digits=5)
		}
		else{
			probINFECTER<-c()
			for(j in 1:nbsimulprobINFECTER){
				probINFECTER<-c(probINFECTER,round(vilprobINFECTER[j],digits=5))
			}
		}

		# \sigma
		sigmadura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@sigma,object@duration)		
		vilsigma<-sigmadura[[1]]
		vilduration<-sigmadura[[2]]
		nbsimulsigma<-length(vilsigma)
		if(nbsimulsigma==1){
			sigma<-round(vilsigma,digits=5)
		}
		else{
			sigma<-c()
			for(j in 1:nbsimulsigma){
				sigma<-c(sigma,round(vilsigma[j],digits=5))
			}
		}
		# \gamma				
		gammadura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@gamma,object@duration)		
		vilgamma<-gammadura[[1]]
		vilduration<-gammadura[[2]]
		nbsimulgamma<-length(vilgamma)
		if(nbsimulgamma==1){
			gamma<-round(vilgamma,digits=5)
		}
		else{
			gamma<-c()
			for(j in 1:nbsimulgamma){
				gamma<-c(gamma,round(vilgamma[j],digits=5))
			}
		}
		# \mu			
		mudura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@mu,object@duration)		
		vilmu<-mudura[[1]]
		vilduration<-mudura[[2]]
		nbsimulmu<-length(vilmu)
		if(nbsimulmu==1){
			mu<-round(vilmu,digits=5)
		}
		else{
			mu<-c()
			for(j in 1:nbsimulmu){
				mu<-c(mu,round(vilmu[j],digits=5))
			}
		}	
		#\phiPHASE			
		phiPHASEdura<-paramdurNewInfec(ville=i,nbVilles=nbVilles,object@phiPHASE,object@duration)		
		vilphiPHASE<-phiPHASEdura[[1]]
		vilduration<-phiPHASEdura[[2]]
		nbsimulphiPHASE<-length(vilphiPHASE)
		if(nbsimulphiPHASE==1){
			phiPHASE<-round(vilphiPHASE,digits=5)
		}
		else{
			phiPHASE<-c()
			for(j in 1:nbsimulphiPHASE){
				phiPHASE<-c(phiPHASE,round(vilphiPHASE[j],digits=5))
			}
		}	
		######
		parami<-c(N=N,S=S,E=E,I=I,R=R,duration=duration,unitTIME=unitTIME,
		nbCONTACT0=nbCONTACT0,nbCONTACT1=nbCONTACT1,probINFECTER=probINFECTER,probVISITER=probVISITER,
		sigma=sigma,gamma=gamma,mu=mu,phiPHASE=phiPHASE)
		names(parami)[1:5]<-c("N","S","E","I","R")	
		print(parami)
		if(i!=nbVilles) cat("\n",fill=60)
	}	
	#coef(object)
}
########################################

