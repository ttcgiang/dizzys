/* 
1. Function allows to calculate the equilibrium values of variables S, E, I, R 
based on the values of parameters.

2. This function permits only to call C code  from R.

3. If building library outside of R package (i.e. for debugging): R CMD SHLIB equilibrium.cpp
*/ 

#include <stdexcept>
#include "equiNewInfec.h"


// From R, we transmit directly the values of parameters to C++. After to make calculations in C++
// we re-transmit result to R.
//
// m0: Average death rate and also Average birth rate
// beta0: Contact rate
// sigma0: Average latent period
// gamma0: Average infectious period
//   
// Result: equilibrium values of S, E, I
SEXP getEquiNewInfec(SEXP mu0,SEXP nbCONTACT0,SEXP probINFECTER,SEXP sigma0,SEXP gamma0){
	double mu, beta, gamma, sigma;
	double cnbCONTACT0, cprobINFECTER;
	//convert the R  value type to the C++ value type
	mu = NUMERIC_VALUE(mu0); 
	cnbCONTACT0 = NUMERIC_VALUE(nbCONTACT0);
	cprobINFECTER = NUMERIC_VALUE(probINFECTER);
	sigma = NUMERIC_VALUE(sigma0);
	gamma = NUMERIC_VALUE(gamma0);
 	
	//parameter holding result 
	SEXP rsei_eq;
   	double *csei_eq;
   	int len = 3;

   	// Allocating storage space:
   	PROTECT(rsei_eq = NEW_NUMERIC(len));
	csei_eq = NUMERIC_POINTER(rsei_eq);
	
	//calculating beta=-Klog(1-c)
	beta = -cnbCONTACT0*log(1-cprobINFECTER);
	double sir_R0 = beta/(gamma+mu); // calculating the value of R0
	if(sigma==INFINITY){// for SIR model, 
		if(sir_R0>1.0){
			csei_eq[0] = 1/sir_R0;// equilibrium value of S
			csei_eq[1] = 0.0; // E is zero
			csei_eq[2] =  mu*(sir_R0 - 1)/beta; //equilibrium value of I
		}
		else{
			 error("The equilibrium value R0 is less than 1, R0 = beta/(gamma+mu)");
		}
			
	}
	else{// for SEIR model,
		csei_eq[0] = (gamma+mu)*(sigma+mu)/(beta * sigma);// equilibrium value of S
		csei_eq[1] = mu*((1/(sigma+mu)) - ((gamma+mu)/(beta * sigma)));// equilibrium value of E
		csei_eq[2] = mu*((beta*sigma - (gamma+mu)*(sigma+mu))/(beta *(gamma+mu)*(sigma+mu)));// equilibrium value of I
	}   	

   UNPROTECT(1);
   return rsei_eq;
}
//the end
