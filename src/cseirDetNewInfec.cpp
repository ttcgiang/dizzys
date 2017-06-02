/*
 *  Created by TRAN Thi Cam Giang on 03/2013.
 *  Goal: calculate the deterministic equations of the SEIR/SIR model.
 *  By associating C++ with R.
 *  we dip R in C++
 */

//Library of R
#include <R.h>
#include <Rdefines.h>
#include "cseirDetNewInfec.h"

extern "C" {

//Define the station of each parameter in the vector 'parms'
static double newparms[9];
#define nbCONTACT0 newparms[0]
#define nbCONTACT1 newparms[1]
#define probINFECTER newparms[2]
#define periDISE newparms[3]
#define phiPHASE newparms[4]
#define mu newparms[5]
#define sigma newparms[6]
#define gamma newparms[7]
#define N newparms[8]


/* initialize */
void initmod(void (* odeparms)(int *, double *))
{
    int nbparames=9;
	odeparms(&nbparames, newparms);
}

/* Derivatives*/
/* There are four ordinal different equations*/
#define time y[0]
#define S y[1]
#define E y[2]
#define I y[3]
#define dtime  ydot[0]
#define dS ydot[1]
#define dE ydot[2]
#define dI ydot[3]

void derivs(int *neq, double *t, double *y, double *ydot,double *yout, int *ip)
{
    //Initialize the values of parameters
	double pi=3.141593;
	dtime = 1;

    // the sinusoidal function of the contact parameter \beta
    double nbCONTACT = nbCONTACT0*(1+nbCONTACT1*cos(2*pi*time/periDISE+phiPHASE));
    double beta= -nbCONTACT * log(1-probINFECTER);

    // the ODE of \dS
	dS =  N*mu - beta*S*I/N - mu*S;	

    if(sigma==INFINITY){// for the SIR model
        // The ODE of \dE
		dE = 0.0;
        // The ODE of \dI
		dI = beta*S*I/N - gamma*I - mu*I;
	}
    else{// for the SEIR model
        // The ODE of \dE
		dE = beta*S*I/N  - sigma*E - mu*E;
        // The ODE of \dI
        dI = sigma*E - gamma*I - mu*I;
	}
}
}
//The end

