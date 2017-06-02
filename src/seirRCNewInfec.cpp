//#include <Rcpp.h>
//using namespace Rcpp;
//#include "varstationaire.h"
#include "stoSEIRNewInfec.h"
#include "initialFUNC.h"

//Base library of R
#include <R.h>
#include <Rdefines.h>

using namespace std;

//arr_rho is matrix or a number
//epsilon is a matrix or a number

extern "C" {

  // Gathering all values of parameters and variables from R to C++
  SEXP getStoSEIRNewInfec(SEXP sigma, SEXP gamma,SEXP mu,SEXP seed,
                          SEXP initS, SEXP initE, SEXP initI, SEXP initR,
                          SEXP nbCONTACT0, SEXP nbCONTACT1,
                          SEXP phiPHASE, SEXP probVISITER, SEXP probINFECTER, SEXP tmax,SEXP nbVilles,
                          SEXP unitTIME, SEXP typeRNG, SEXP periDISE)
  {
    try{

      //Number of cities in the metapopulation
      if(!isReal(nbVilles) || length(nbVilles)!=1)
        error("invalid number of cities.");
      int cnbVilles = INTEGER_VALUE(nbVilles);

      //-sigma 0.125   //sigma
      if (!isReal(sigma) || length(sigma)!=1) {
        error("invalid parameter sigma");
      }
      double csigma = NUMERIC_VALUE(sigma);

      //-gamma 0.2  //gamma
      if (!isReal(gamma) || length(gamma)!=1) {
        error("invalid parameter gamma");
      }
      double cgamma = NUMERIC_VALUE(gamma);

      //-mu 0.000039   //mu
      if (!isReal(mu) || length(mu) !=1) {
        error("invalid parameter mu");
      }
      double cmu = NUMERIC_VALUE(mu);

      //-seed 10    //number 'seed'
      if (!isReal(seed)  ||  length(seed) != 1) {
        error("invalid value of seed,it should be a real number");
      }
      double cseed=NUMERIC_VALUE(seed);

      //-initS, vector includes all initial values of variable S
      if (!isVector(initS)  ||  !isReal(initS)) {
        error("invalid vector of initial values of S (initS)");
      }
      double *rinitS = REAL(initS);
      double* cinitS; create1D(cnbVilles,0.0,cinitS);
      for(int i=0; i<cnbVilles; i++) cinitS[i] = rinitS[i];

      //-initE, vector includes all initial values of variable E
      if (!isVector(initE)  ||  !isReal(initE)) {
        error("invalid vector of initial values of E (initE)");
      }
      double *rinitE = REAL(initE);
      double* cinitE; create1D(cnbVilles,0.0,cinitE);
      for(int i=0; i<cnbVilles; i++) cinitE[i] = rinitE[i];

      //-initI, vector includes all initial values of variable I
      if (!isVector(initI)  ||  !isReal(initI)) {
        error("invalid vector of initial values of I (initI)");
      }
      double *rinitI = REAL(initI);
      double* cinitI; create1D(cnbVilles,0.0,cinitI);
      for(int i=0; i<cnbVilles; i++) cinitI[i] = rinitI[i];

      //-initR, vector includes all initial values of variable R
      if (!isVector(initR)  ||  !isReal(initR)) {
        error("invalid vector of initial values of R (initR)");
      }
      double *rinitR = REAL(initR);
      double* cinitR; create1D(cnbVilles,0.0,cinitR);
      for(int i=0; i<cnbVilles; i++) cinitR[i] = rinitR[i];

      //-nbCONTACT0 300
      if (!isReal(nbCONTACT0)  ||  length(nbCONTACT0) != 1) {
        error("invalid value of nbCONTACT0,it should be a real number");
      }
      double cnbCONTACT0=NUMERIC_VALUE(nbCONTACT0);

      //-nbCONTACT1 0.1
      if (!isReal(nbCONTACT1)  ||  length(nbCONTACT1) != 1) {
        error("invalid value of nbCONTACT1,it should be a real number");
      }
      double cnbCONTACT1=NUMERIC_VALUE(nbCONTACT1);

      //-phiPHASE  is a vector of n elements
      if (!isVector(phiPHASE)  ||  !isReal(phiPHASE)) {
        error("invalid vector of phase (phiPHASE)");
      }
      double *rphiPHASE = REAL(phiPHASE);
      double* cphiPHASE; create1D(cnbVilles,0.0,cphiPHASE);
      for(int i=0; i<cnbVilles; i++) cphiPHASE[i] = rphiPHASE[i];

      //-probVISITER 0.01
      if (!isReal(probVISITER)  ||  length(probVISITER) != 1) {
        error("invalid value of probVISITER,it should be a real number");
      }
      double cprobVISITER=NUMERIC_VALUE(probVISITER);

      //-probINFECTER 0.01
      if (!isReal(probINFECTER)  ||  length(probINFECTER) != 1) {
        error("invalid value of probINFECTER,it should be a real number");
      }
      double cprobINFECTER=NUMERIC_VALUE(probINFECTER);


      // time unit to gather the values of variables SEIR
      if (!isReal(unitTIME)  ||  length(unitTIME) != 1) {
        error("invalid value of unitTIME,it should be a real number");
      }
      double cunitTIME=NUMERIC_VALUE(unitTIME);

      //simulation time tmax
      if (!isReal(tmax)  ||  length(tmax) != 1) {
        error("invalid value of simulation time,it should be a real number");
      }
      double ctmax=NUMERIC_VALUE(tmax);

      //type of the random number generator 1: Yann, 0: C++
      if (!isReal(typeRNG)  ||  length(typeRNG) != 1) {
        error("invalid value of type generating a random,it should be 0 or 1");
      }
      double ctypeRNG=NUMERIC_VALUE(typeRNG);

      // Disease period
      if (!isReal(periDISE)  ||  length(periDISE) != 1) {
        error("invalid value of period, it should be a real number");
      }
      double cPeriDISE=NUMERIC_VALUE(periDISE);

      // for metapopulation
      double ** arr_probVISITER = calculerProbVISITER(cnbVilles,cprobVISITER);

      // xay dung ma tran INFECTER
      double **arr_probINFECTER; create2D(cnbVilles,cnbVilles,cprobINFECTER,arr_probINFECTER);

      // finding the stationary distribution
      //vector stationary line (statline)
      //if (!isVector(statline)  ||  !isReal(statline)) {
      //  error("invalid vector of stationary line (statline)");
      //}

      // initialiser les valeurs des variables
      int nvar= 6;
      unsigned long **valSEIR0;
      valSEIR0 = new unsigned long*[cnbVilles];
      for (int i = 0; i <  cnbVilles; i++)
      {
        valSEIR0[i] = new unsigned long[nvar];
        valSEIR0[i][0]= 0;
        valSEIR0[i][iS]= cinitS[i];
        valSEIR0[i][iE]= cinitE[i];
        valSEIR0[i][iI]= cinitI[i];
        valSEIR0[i][iR]= cinitR[i];
        valSEIR0[i][iN]= cinitS[i]+cinitE[i]+cinitI[i]+cinitR[i];
      }
      /*
      cout<<"cnbVilles"<<cnbVilles<<endl;
      cout<<"valSEIR0"<<endl;
      for(int i=0; i<cnbVilles; i++){
      for (int j=0; j<nvar; j++)
      cout<<"    " << valSEIR0[i][j];
      cout<<endl;
      }
      */


      // xay dung ma tran VENIR de la ville k
      double **arr_probVENIRk;  create2D(cnbVilles,cnbVilles,0, arr_probVENIRk);
      calculerProbVENIRk(cnbVilles,valSEIR0,arr_probVISITER,arr_probVENIRk);

      // xay dung ma tran NUMBER OF CONTACT in forcing
      // double nbCONTACT0 = valParSIM[inbCONTACT0];
      //double nbCONTACT1 = valParSIM[inbCONTACT1];

      double **arrNbCONTACT0 = calculerNbCONTACT(cnbVilles,cnbCONTACT0,1);
      double **arrNbCONTACT1 = calculerNbCONTACT(cnbVilles,cnbCONTACT1,1.0);

      //SIMULATION
      // Doing simulation by using the function 'simulerSEIR'
      resMETAPOP resTOTAL =  simulerSEIRNewInfec(cnbVilles,csigma,cgamma,cmu,valSEIR0,
                                                 0.0,0.0,cphiPHASE,cseed,cunitTIME,ctmax,ctypeRNG,cPeriDISE,
                                                 arr_probVISITER,arrNbCONTACT0,arrNbCONTACT1,arr_probINFECTER,arr_probVENIRk);

      vector<vector<vector<unsigned long> > > resSEIRMeta = resTOTAL.tableResVilles;
      //RESULTAT
      // Gathering the result, then send this from C++ to R
      int size2D = resSEIRMeta[0].size();
      SEXP listMETAPOP, list_names_villes;
      char names_villes[255];
      int nbElements=cnbVilles+2;
      // objects in our list:
      PROTECT(list_names_villes = allocVector(STRSXP,nbElements));
      for(int i=0; i<nbElements;i++){
        sprintf(names_villes,"pop%d",i);
        if(i==cnbVilles)  SET_STRING_ELT(list_names_villes,cnbVilles,mkChar("localExtPOP"));
        else if(i==(cnbVilles+1))  SET_STRING_ELT(list_names_villes,cnbVilles+1,mkChar("RecolPOP"));
        else SET_STRING_ELT(list_names_villes,i,mkChar(names_villes));
      }

      PROTECT(listMETAPOP = allocVector(VECSXP, nbElements));
      SEXP myT,myS,myE,myI,myR,myN,myP, listSEIR, list_names;
      double *valeurT, *valeurS, *valeurE,*valeurI,*valeurR,*valeurN,*valeurP;
      const char *names[7] = {"time(day)", "S", "E", "I", "R", "N","P"};

      for(int idxVille=0; idxVille<cnbVilles;idxVille++){
        // creating an double vector:
        PROTECT(myT = NEW_NUMERIC(size2D));
        valeurT = NUMERIC_POINTER(myT);

        PROTECT(myS = NEW_NUMERIC(size2D));
        valeurS = NUMERIC_POINTER(myS);

        PROTECT(myE = NEW_NUMERIC(size2D));
        valeurE = NUMERIC_POINTER(myE);

        PROTECT(myI = NEW_NUMERIC(size2D));
        valeurI = NUMERIC_POINTER(myI);

        PROTECT(myR = NEW_NUMERIC(size2D));
        valeurR = NUMERIC_POINTER(myR);

        PROTECT(myN = NEW_NUMERIC(size2D));
        valeurN = NUMERIC_POINTER(myN);

        PROTECT(myP = NEW_NUMERIC(size2D));
        valeurP = NUMERIC_POINTER(myP);

        // Creating a character string vector of the "names" attribute of the objects in out list:
        PROTECT(list_names = allocVector(STRSXP,7));
        for(int i = 0; i < 7; i++)
          SET_STRING_ELT(list_names,i,mkChar(names[i]));

        // Creating a list with 2 vector elements:
        PROTECT(listSEIR = allocVector(VECSXP, 7));

        for(int j=0; j<size2D; j++){
          valeurT[j]= resSEIRMeta[idxVille][j][0];
          valeurS[j]= resSEIRMeta[idxVille][j][1];
          valeurE[j]= resSEIRMeta[idxVille][j][2];
          valeurI[j]= resSEIRMeta[idxVille][j][3];
          valeurR[j]= resSEIRMeta[idxVille][j][4];
          valeurN[j]= resSEIRMeta[idxVille][j][5];
          valeurP[j]= resSEIRMeta[idxVille][j][6];
        }

        // attaching myint vector to list:
        SET_VECTOR_ELT(listSEIR, 0, myT);
        // attaching mydouble vector to list:
        SET_VECTOR_ELT(listSEIR, 1, myS);
        SET_VECTOR_ELT(listSEIR, 2, myE);
        SET_VECTOR_ELT(listSEIR, 3, myI);
        SET_VECTOR_ELT(listSEIR, 4, myR);
        SET_VECTOR_ELT(listSEIR, 5, myN);
        SET_VECTOR_ELT(listSEIR, 6, myP);

        // and attaching the vector names:
        setAttrib(listSEIR, R_NamesSymbol, list_names);
        SET_VECTOR_ELT(listMETAPOP, idxVille, listSEIR);
        UNPROTECT(9);
      }

      //localEXT
      vector<vector<double > > resExtLocal = resTOTAL.tableExtLocal; //vector total for the data seir of each city
      SEXP names_loclExt;
      char tpnames[255];
      // objects in our list:
      PROTECT(names_loclExt = allocVector(STRSXP,cnbVilles));
      for(int i=0; i<cnbVilles;i++){
        sprintf(tpnames,"localExt%d",i);
        SET_STRING_ELT(names_loclExt,i,mkChar(tpnames));
      }

      double *valeurExt;
      SEXP listLocalExt,myLocalExt;
      PROTECT(listLocalExt = allocVector(VECSXP, cnbVilles));

      for(int idxVil=0; idxVil<cnbVilles; idxVil++){
        int sizeExtVil = resExtLocal[idxVil].size();
        PROTECT(myLocalExt = NEW_NUMERIC(sizeExtVil));
        valeurExt = NUMERIC_POINTER(myLocalExt);
        for(int j=0; j<sizeExtVil; j++){
          valeurExt[j]=resExtLocal[idxVil][j];
        }
        // attaching myint vector to list:
        setAttrib(listLocalExt, R_NamesSymbol, names_loclExt);
        SET_VECTOR_ELT(listLocalExt, idxVil, myLocalExt);
        UNPROTECT(1);
      }
      SET_VECTOR_ELT(listMETAPOP,cnbVilles, listLocalExt);
      UNPROTECT(2);

      // Recol
      vector<vector<double > > resRecolLocal = resTOTAL.tableRecolLocal;
      SEXP namesRecol;
      // objects in our list:
      PROTECT(namesRecol = allocVector(STRSXP,cnbVilles));
      for(int i=0; i<cnbVilles;i++){
        sprintf(tpnames,"recol%d",i);
        SET_STRING_ELT(namesRecol,i,mkChar(tpnames));
      }
      double *valeurRecl;
      SEXP listRecol,myRecol;
      PROTECT(listRecol = allocVector(VECSXP, cnbVilles));
      for(int idxVil=0; idxVil<cnbVilles; idxVil++){
        int sizeExtVil = resRecolLocal[idxVil].size();
        PROTECT(myRecol = NEW_NUMERIC(sizeExtVil));
        valeurRecl = NUMERIC_POINTER(myRecol);
        for(int j=0; j<sizeExtVil; j++){
          valeurRecl[j]=resRecolLocal[idxVil][j];
        }
        // attaching myint vector to list:
        setAttrib(listRecol, R_NamesSymbol, namesRecol);
        SET_VECTOR_ELT(listRecol, idxVil, myRecol);
        UNPROTECT(1);
      }
      SET_VECTOR_ELT(listMETAPOP,cnbVilles+1, listRecol);
      UNPROTECT(2);

      setAttrib(listMETAPOP,R_NamesSymbol,list_names_villes);
      UNPROTECT(2);

      // freeing the memory of pointer
      /*

      for (int i = 0; i <  cnbVilles; i++)
      {
      delete []valSEIR0[i];
      }
      delete []valSEIR0;
      */

      return listMETAPOP;

    }catch (exception &e) {
      error(e.what());
      return(0);
    }
  }
}

/***************The end*********************/
