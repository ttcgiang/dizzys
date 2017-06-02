/*
 *
 *  Created by yann chevaleyre on 24/11/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "randYANN.h"
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <assert.h>

using namespace std;

/* Below: functions to generate random variables with a normal distribution */

/* File:      gasdev.C
 Date:      12/31/01  (written fall 2001)
 Author:    C. Paciorek
 Purpose:   generate a standard normally-distributed (i.e. Gaussian) random variable (scalar)
 Compile:   g++ -c -Wall -ansi
 Usage:     N/A 
 Arguments: N/A
 Requires:  see just below
 Reference: Numerical Recipes in C - code taken directly from the source
 */


double normal_distribution() {  // generates 2 random deviates at a time; idum is seed
	static int iset=0;
    static int seed = 0;
    int  *idum = &seed;			/*idum is the random seed*/
	static float gset; 
	float fac,rsq,v1,v2; 

	if (*idum < 0) iset=0; 
	if (iset == 0) {  //We don't have an extra deviate handy, so 
		do { v1=2.0*ran1(idum,GENERATE,(char*)"a")-1.0; // pick two uniform numbers in the square extending from -1 to +1 in each direction, 
			v2=2.0*ran1(idum,GENERATE,(char*)"a")-1.0; 
			rsq=v1*v1+v2*v2; // see if they are in the unit circle,
		} 
		while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again. 
		fac=sqrt(-2.0*log(rsq)/rsq); /* Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time. */
		gset=v1*fac; 
		iset=1; // Set ag. 
		return (double) v2*fac; } 
	else { // We have an extra deviate handy, 
		iset=0; // so unset the ag, 
		return (double) gset; // and return it. 
	} 
}


/* File:      ran1.C
 Date:      01/04/02  (written fall 2001)
 Author:    C. Paciorek 
 Purpose:   random uniform generator
 Compile:   g++ -c -Wall -ansi
 Usage:     N/A
 Arguments: N/A
 Requires:  see just below
 Reference: Numerical Recipes in C  - code taken directly from the source and modified by C. Paciorek to read and write the shuffle table so that one can restart a simulation with the same random sequence as one stopped with
 */

#define IA 16807 
#define IM 2147483647 
#define AM (1.0/IM) 
#define IQ 127773 
#define IR 2836 
#define NTAB 32 
#define NDIV (1+(IM-1)/NTAB) 
#define EPS 1.2e-7 
#define RNMX (1.0-EPS) 

float ran1(int *idum,int mode,char* line){ /* "Minimal" random number generator of Park and Miller with Bays-Durham shu e and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter idum between successive deviates in a sequence. RNMX should approximate the largest oating value that is less than 1. */
	/* if mode==READ or mode==WRITE, the function reads or writes the shuffle table to line */
	int j; int k; 
	static int iy=0; 
	static int iv[NTAB]; 
	float temp; 
	
	if(mode==GENERATE){
		if (*idum <= 0 || !iy) { //Initialize. 
			if (-(*idum) < 1) *idum=1; //Be sure to prevent idum = 0. 
			else *idum = -(*idum); 
			for (j=NTAB+7;j>=0;j--) { // Load the shu e table (after 8 warm-ups). 
				k=(*idum)/IQ; 
				*idum=IA*(*idum-k*IQ)-IR*k; 
				if (*idum < 0) *idum += IM; 
				if (j < NTAB) iv[j] = *idum; } 
			iy=iv[0]; } 
		k=(*idum)/IQ; //Start here when not initializing. 
		*idum=IA*(*idum-k*IQ)-IR*k; //Compute idum=(IA*idum) % IM without over- ows by Schrage's method. 
		if (*idum < 0) *idum += IM; 
		j=iy/NDIV; //Will be in the range 0..NTAB-1. 
		iy=iv[j]; //Output previously stored value and re ll the shu e table. 
		iv[j] = *idum; 
		if ((temp=AM*iy) > RNMX) return RNMX; //Because users don't expect endpoint values. 
		else return temp; 
	} 
	return(0);
}

/* --- Unit Test ---*/
/* testing of random normal distribution generator - displays a normal distribution on screen */
void unit_test_normal_distribution()
{
	int tab[64],maxi=0,i,j;
	
	for (i=0; i < 64 ; i++)
		tab[i] = 0;
	
	for (i=0; i < 100000; i++) {
		int j = normal_distribution()*8+32;
		if (j<0)	j=0;
		if (j>=64)	j=63;
		tab[j]++;
		if (tab[j] > maxi)
			maxi = tab[j];
	}
	
	for (i=30; i >= 0; --i) {
		std::cout << "|";
		for (j=0; j < 64; j++) 
			if (tab[j] > (int)(i*maxi/30.0))
				std::cout << "#";
			else
				std::cout << " ";
		
		std::cout << "|\n";
	}
}

/********* The end ******/


