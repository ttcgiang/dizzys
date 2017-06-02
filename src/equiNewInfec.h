//bibliothèque sous C++
#include <stdio.h>
#include <stdlib.h>
//bibliothèque sous R
#include <R.h>
#include <Rdefines.h>
  
using namespace std;
extern "C" {

/*
Goal:
	Function allows to calculate the equilibrium values of variables S, E, I, R based on the values of parameters.
*/
SEXP getEquiNewInfec(SEXP mu0,SEXP nbCONTACT0,SEXP probINFECTER,SEXP sigma0,SEXP gamma0);

/***********************************************************************************************/
}


