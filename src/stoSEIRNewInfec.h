#ifndef MODELNewInfec_H
#define MODELnewInfec_H
//bibliothèque de C++
#include <ctime>
#include <stdio.h>
#include <vector>
#include <map>
#include <stdlib.h>
#include <cmath>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <stdexcept>
#include "index.h"
#include "randYANN.h"
#include "initialFUNC.h"

using namespace std;
// here, we do functions of YANN and GIANG

struct resMETAPOP{
    vector<vector<vector<unsigned long> > > tableResVilles;
    vector<vector<double > > tableExtLocal; //vector total for the data seir of each city
    vector<vector<double > > tableRecolLocal; //vector total for the data seir of each city
};

/**************************************************/
   /*
  TROUVER double** probVisit
Object:
    Fonction crée un table en deux dimensions qui contient les taux de couplage entre ville i et j,
    qui montre la probabilite pour que un ind x de la ville i visite la ville j.

Entrée:
    rho0: taux de couplage initiale

*/
  double** calculerProbVISITER(int nbVilles, double prbVisit0);

   /******************************************/
      /* TROUVER
       * Fonction  calculer les valeurs de \xi
        \xi qui montre la probabilite
        pour que un indi y rencontre un indi x dans ville j qui vient de ville i
      */
  void calculerProbVENIRk(int nbVilles, unsigned long **y,
                              double **arr_probVISIT, double**arr_probVENIRk);
 /******************************************************************/

/*
 Goal: calculating the contact rate for each city following the fomula YANN et GIANG

 Input:
    nbVilles : number of cities in the metapopulation
    prbVisit : \rho_{i,j}, the probability that an individual from subpopulation i
  visits subpopulation j.
    mtVarb : matrix of all variables S, E, I, R at the time t.
    nbCont : \kappa_{j} is the average number of contacts per unit of time a susceptible will have when visiting city j.
    prbInf : c_{i,k} is the probability that a susceptible individual native from i being in contact with another infected individual native from k
  gets infected.
    prbFromk: \xi_{jk} refers to the probability that an individual y meeting x in C_{j} comes from C_{k}.

Output:
    a vector of the contact rates for all cities in the metapopulation.
   */

void calLamdaNewInfec(int nbVilles, double** prbVisit, unsigned long** mtVarb,
                          double** nbCont, double** prbInf, double** prbFromk, double *&lamda);

 /******************************************************************/
 /*
  * GOAL : calculating the initial values of all variables for all cities.
  * INPUT:
  *     nbVilles: number of cities in the metapopulation
  *     arr_statDIS : a stationary distribution with nbVilles elements
  *     S0 : initial number of susceptibles in the metapopulation
  *     E0 : intitial number of exposed in the metapopualtion
  *     I0 : initial number of infected in the metapopulation
  *     R0 : initial number of recovered in the metapopulation
  * OUTPUT:
  *     a matrix in two dimentions for all initial values of all variavles for all cities in the metapopulation
  *
  */

void initializeStatSEIRNewInfec(int nbVilles, double* statMatrix, unsigned long S0, unsigned long I0, unsigned long E0, unsigned long R0,
                        unsigned long **&y);
 /******************************************************************/
   /*
    * GOAL : calculating the number of contact (with the proposition: The coefficient \kappa should also depend on i,
    *  because an individual native from city i meets more people in his own city than abroad
    *  (\kappa_{i,i}>\kappa_{i,j}).
    * INPUT :
    *   nbVilles : number of cities
    *   nbCONTACT : the initial value of the number of CONTACT
    *   power : the number of time here \kappa_{i,i} = power * \kappa_{i,j}
    * */
  double** calculerNbCONTACT(int nbVilles, double nbCONTACT0, double multiplicative);
   /******************************************************************/
  /*
  Calcule le nomber de contact following the forcing
   avec K_i(t) = K0 *(1+ K1*cos(2*pi*t + phi_i).

   */
   void calculateNbCONTACTForcing(int nbVilles, double** nbCONTACT0, double** nbCONTACT1,
                         double periCONTACT, double t, double*phi, double **arr_NbCONTACT);
      /******************************************************************/
    /*
    Objectif:
        Fonction est utilisée dans le cas, les valeurs des parameters, et des variables sont différents pour chaque ville.
        Simuler le modèle SEIR avec plusieurs populations et beta sinusoidal.

    Entrée:
        nbVilles: nombre de sous-populations.
            gamma: vecteur de taux de guérison (/jour).
            sigma: vecteur de inverse de la durée moyenne d'infection (/jour).
            beta0: vecteur de valeur moyenne du taux de contact (/ind/jour).
            beta1: vecteur de phase du forçage
            tmax: durée maximale de la simulation (jours).
            unitTemps: pas de temps auquel les valeurs des variables sont enregistrées 				(jours).
            S0: valeur initiale de susceptibles dans une sous-population.
            E0: valeur initiale d'exposé dans une sous-population.
            I0: valeur initiale d'infectieux dans une sous-population.
            R0: valeur initiale de guéris dans une sous-population.
            N0: valeur initiale de la taille d'une sous-population.
            rmu: taux de mortalité égal au taux de natalité (/jour).
            epsilon: [0;1] taux d'infection de l'extrieur.
            rho: [0;1] taux de couplage entre les sous-populations
            graine: graine du générateur de nombres aléatoires.
            phiMIN: phase minimum dans la métapopulation.
            phiMAX: phase maximun dans la métapopulation.

    Sortie:
        Un vecteur en trois dimensions:
            la dimension 1 pour temps
            la dimension 2 pour valeurs des variables
            la dimension finalles pour le nombre de villes
    */
     resMETAPOP simulerSEIRNewInfec(
                                  // for model SEIR
                                  int nbVilles,double sigma,double gamma,double mu, unsigned long **valSEIR0,
                                  double beta0,double beta1,double *phi,
                                  // for simulation
                                  double seed,double unitTIME, double tmax, double typeRNG,double periDISE,
                                  // for metapopulation
                                  double** arr_probVISITER,double** arr_nbCONTACT0,double** arr_nbCONTACT1, double** arr_probINFECTER, double** arr_probVENIRk);
          /******************************************************************/

#endif // MODELNewInfec_H
