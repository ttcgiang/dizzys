/*
The goal of the file is:
1. d'attribuer des noms à des index des parameters, des variables. 
Cela nous faciliter de programmer. Au lieu de retenir les index des parameters et des variables, on retient ses noms.
2. de définir des valeurs pour des variables nécessaires.
*/

#ifndef INDEX_H
#define INDEX_H

// extern "C": parce que R ne comprend que la langage C, malgré que notre code est sous C++ 
//

// DESCRIRE LES INDEX DES PARAMETRES
/**********************************************************************************************/
// noms pour les index des variables
const int it=0,iS=1,iE=2,iI=3,iR=4,iN=5,iP=6;

// noms pour les index des parametres
const int itmax = 0, inbVilles = 1,
          isigma = 2, igamma = 3,imu = 4,iepsilon = 5,
          irho = 6,iunitTIME = 7, iseed = 8,
          iS0 = 9, iI0 = 10, iE0 = 11, iR0 = 12, iN0 = 13,
          ibeta0 = 14, ibeta1 = 15, iphiMIN = 16, iphiMAX =17, inbSimu =18,
          iPeriDISE = 19, itypeRNG=20,
          //
         //for metepopulation of YANN et GIANG
         iprobVISITER=21, // prob visiting other city
         inbCONTACT0=22, //nbCONTACT0 is the mean transmission rate
         inbCONTACT1=23, //and nbCONTACT1 controls the amplitude of seasonality
         inbMulCONTACT=24,// multiplier, Kii=A*Kij
         iprobINFECTER=25; // prob of transmission of disease
         //iprobVENIRk=24;//prob of comming from city k

// définir le numéro pi
#define Pi 3.141592653589

// noms pour des générateurs de nombres aléatoires
const int GOOD=0; // utiluser le générateur de C++
const int FAST=1; // utiliser le générateur de Yann
const int GENERATE=0; // parameter dans la fonction de Yann


/**********************************************************************************************/
#endif // INDEX_H
