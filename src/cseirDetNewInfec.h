#include <R.h>
#include <Rdefines.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
extern "C" {
/* initializer */
void initmod(void (* odeparms)(int *, double *));

/* Derivatives */
void derivs(int *neq, double *t, double *y, double *ydot,double *yout, int *ip);
}

/***********************************************************************************************/



