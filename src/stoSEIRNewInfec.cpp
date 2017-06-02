#include "stoSEIRNewInfec.h"

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

void calLamdaNewInfec(int nbVilles, double** prbVisit,unsigned long** mtVarb,
                          double** nbCont, double** prbInf, double** prbFromk, double *&lamda){
   // double* resContRate; create1D(nbVilles,0.0, resContRate);
    for(int i=0; i<nbVilles; i++){
        lamda[i] = 0.0;
        for(int j=0; j<nbVilles; j++){
            double sumPrbContInf=0.0;
            for(int k=0; k<nbVilles; k++){
               // cout<<"k=="<<k<<" mtVarb[k][iI]="<<mtVarb[k][iI]<<" mtVarb[k][iN]= "<<mtVarb[k][iN]<<endl;
                //cout<<"prbInf[i][k]="<<prbInf[i][k]<<endl;
                  //cout<<"prbFromk[j][k]="<<prbFromk[j][k]<<endl;

                sumPrbContInf += (mtVarb[k][iI]*prbInf[i][k]*prbFromk[j][k])/mtVarb[k][iN];
                 // cout<<"sumPrbContInf="<<sumPrbContInf<<endl;

            }
            double tplog= log(1-sumPrbContInf);
           // cout<<"tplog="<<tplog<<endl;
            //cout<<"nbContij="<<nbCont[i][j]<<endl;
            lamda[i] += - prbVisit[i][j]*nbCont[i][j]*tplog;
        }
    }
//for(int i=0; i<nbVilles; i++){cout<<"LAMAMAMAM="<<lamda[i]<<endl;}
    //return (resContRate);
}
/***************************************************************************/
   /*
  TROUVER double** probVisit
Object:
    Fonction crée un table en deux dimensions qui contient les taux de couplage entre ville i et j,
    qui montre la probabilite pour que un ind x de la ville i visite la ville j.

Entrée:
    rho0: taux de couplage initiale

*/
  double** calculerProbVISITER(int nbVilles, double prbVisit0){
      double **arr_prbVisit;
      arr_prbVisit = new double*[nbVilles];
      for (int i = 0; i<nbVilles; i++)
      {
          arr_prbVisit[i] = new double[nbVilles];
          for(int j = 0; j< nbVilles ; j++) {
              double totaltp= (double)((nbVilles-1)*prbVisit0);
              //\rho_ii
              if(j==i) arr_prbVisit[i][j] = 1.0 -totaltp;
              else
              arr_prbVisit[i][j] = prbVisit0;
               //   cout<<"arr_prbVisit"<<i<<j<<" =   "<< arr_prbVisit[i][j];
          }
       //  cout<<endl;
      }
      return arr_prbVisit;
  }

 /******************************************************************/
  /* TROUVER
   * Fonction  calculer les valeurs de \xi
    \xi qui montre la probabilite
    pour que un indi y rencontre un indi x dans ville j qui vient de ville i
  */
  void calculerProbVENIRk(int nbVilles, unsigned long **y,
                          double **arr_probVISIT, double**arr_probVENIRk){
      for (int i = 0; i<nbVilles; i++)
      {
         //calculer le demominateur
          double denom=0.0;
          for(int k=0; k<nbVilles; k++){
              //cout<<"y[k][iN] =="<<y[k][iN]<<endl;
              //cout<<"arr_probVISIT =="<<arr_probVISIT[k][i]<<endl;
              denom = denom + y[k][iN]*arr_probVISIT[k][i];
          }
         // cout<<"i="<<i<<"  denom == "<<denom<<endl;
          //
          for(int j = 0; j< nbVilles ; j++) {
              arr_probVENIRk[i][j]= y[j][iN]*arr_probVISIT[j][i]/denom;
          }
      }
  }
  /**************************************************************************/
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

    void initializeStatSEIRNewInfec(int nbVilles, double* statMatrix, unsigned long S0, unsigned long E0,unsigned long I0, unsigned long R0,
                            unsigned long **&y){
        for(int i=0; i<nbVilles; i++){
            y[i][iS]= statMatrix[i]*S0;
            y[i][iE]= statMatrix[i]*E0;
            y[i][iI]= statMatrix[i]*I0;
            y[i][iR]= statMatrix[i]*R0;
            y[i][iN]= y[i][iS]+y[i][iE]+y[i][iI]+y[i][iR];
        }
    }

  /******************************************************************/
    /*
     * GOAL : calculating the number of contact (with the proposition: The coefficient \kappa should also depend on i,
     *  because an individual native from city i meets more people in his own city than abroad
     *  (\kappa_{i,i}>\kappa_{i,j}).
     * INPUT :
     *   nbVilles : number of cities
     *   nbCONTACT : the initial value of the number of CONTACT
     *   multiplicative : the number of time here \kappa_{i,i} = multiplicative * \kappa_{i,j}
     * */
   double** calculerNbCONTACT(int nbVilles, double nbCONTACT0, double multiplicative){
       double **arr_NbCONTACT;
       arr_NbCONTACT = new double*[nbVilles];
       for (int i = 0; i<nbVilles; i++)
       {
           arr_NbCONTACT[i] = new double[nbVilles];
           for(int j = 0; j< nbVilles ; j++) {
                if(j!=i) arr_NbCONTACT[i][j]=nbCONTACT0;
                else
                    arr_NbCONTACT[i][j]=multiplicative*nbCONTACT0;
           }
       }
       return arr_NbCONTACT;

   }
   /******************************************************************/
   /*
   Calcule le nomber de contact following the forcing
    avec K_i(t) = K0 *(1+ K1*cos(2*pi*t + phi_i).

    */
     void calculateNbCONTACTForcing(int nbVilles, double** nbCONTACT0, double** nbCONTACT1,
                          double periCONTACT,double t, double*phi, double** arr_NbCONTACT){
         double coef=2*Pi*t/periCONTACT;

         //double** arr_NbCONTACT;
         //arr_NbCONTACT = new double*[nbVilles];

         for (int i = 0; i<nbVilles; i++)
         {
            // arr_NbCONTACT[i] = new double[nbVilles];
             for(int j = 0; j< nbVilles ; j++) {
/*
                 cout<<"t="<<t<<endl;
                 cout<<"i="<<i<<"  j= "<<j<<endl;
                 cout<<"nbCONTACT0="<<nbCONTACT0[i][j]<<endl;
                 cout<<"nbCONTACT1="<<nbCONTACT1[i][j]<<endl;
                 cout<<"coef="<<coef<<endl;
                 cout<<"phi="<<phi[i]<<endl;
*/
                 arr_NbCONTACT[i][j]= nbCONTACT0[i][j]*(1 + nbCONTACT1[i][j] *cos(coef + phi[j]));
             //   cout<<"arr_NbCONTACT="<<arr_NbCONTACT[i][j]<<endl;
             }

         }

         // return(arr_NbCONTACT);
     }
     /***************************************************/
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
                                 double** arr_probVISITER,double** arr_nbCONTACT0,double** arr_nbCONTACT1, double** arr_probINFECTER, double** arr_probVENIRk){

       /*
       cout<<"trong code C++"<<endl;
	 cout<<"GGLDLDLKDKJDDJDJKJDKDJKJKDJKJK"<<endl;

       for(int i=0; i<nbVilles; i++){
             cout<<"phi=    "<<phi[i]<<endl;
         }


       for(int i=0; i<nbVilles; i++){
           for (int j=0; j<6; j++)
               cout<<"    " << valSEIR0[i][j];
           cout<<endl;
       }
*/
 
       // for model SEIR

       // for simulation

       // Parameters for pour SIMULATION
       // Number of varibales = 5 (S, E, I, R, N);
       // Number of events = 9
         int nevent = 9; //nbVarSEIRN = 6;
       //unsigned long S0 = S, E0 =E, I0 =I, R0 = R;
       double p1 = 0.0, p2 = 0.0, p3 = 0.0;
       double t = 0.0, tstep = 0.0;
       double fsumVilles = 0.0;
       int event = 0, nextville = 0;
       double pointTemps = unitTIME;
       //Initializing necessaire vectors
       double *lamda = new double[nbVilles];
       int *matrCumulI = new int[nbVilles];
       for(int i=0; i<nbVilles; i++) matrCumulI[i]=0.0;
       double** arr_nbCONTACTFORCING; arr_nbCONTACTFORCING = new double*[nbVilles];
       for (int i = 0; i <  nbVilles; i++) {
           arr_nbCONTACTFORCING[i] = new double[nbVilles];
           for(int j = 0; j< nbVilles ; j++) arr_nbCONTACTFORCING[i][j] = 0.0;
       }
       calculateNbCONTACTForcing(nbVilles,arr_nbCONTACT0,arr_nbCONTACT1,periDISE,t,phi,arr_nbCONTACTFORCING);
       //double *beta_sinusoidal;
       //double** arr_xi;create2D(nbVilles,nbVilles,0.0,arr_xi);

       // Saving the values of all variables and all propensity functions at moment t
       unsigned long **y;
       double **f;
       f = new double*[nbVilles];
       for (int i = 0; i <  nbVilles; i++)
       {
           f[i] = new double[nevent];
           for(int j= 0; j< nevent ; j++) f[i][j] = 0.0;
       }
       //Initialising the values of variables
       // for metapopulation
       y = valSEIR0;


       //Creating a table that saves the values of variables when simulation runs
       vector< unsigned long> tpmSEIRNP;
       vector<vector<vector<unsigned long> > > tableResVilles(nbVilles); //vector total for the data seir of each city
       for(int i=0; i<nbVilles; i++){
           tpmSEIRNP.push_back(0.0);tpmSEIRNP.push_back(y[i][iS]); tpmSEIRNP.push_back(y[i][iE]); tpmSEIRNP.push_back(y[i][iI]);tpmSEIRNP.push_back(y[i][iR]);tpmSEIRNP.push_back(y[i][iN]);tpmSEIRNP.push_back(0);
           tableResVilles[i].push_back(tpmSEIRNP);
           tpmSEIRNP.clear();
       }

       vector<vector<double > > tableExtLocal(nbVilles); //vector total for the data seir of each city
       vector<vector<double > > tableRecolLocal(nbVilles); //vector total for the data seir of each city
       vector<bool> vecFlagExtLocal (nbVilles,true);

       double xpro = 0.0, fville = 0.0, fprop =0.0;
       //PROGRAM starts.......................................................

       //Calculating the time where the program runs according to the time of CPU
       clock_t start, end;
       start = clock();
       static int seed0 = seed;
       int  *idum = &seed0;

       //number 'seed'
       srand(seed);
       //simulation
       while(t<=tmax){
             // recalculating the matrix "COME from the city k"

            calculerProbVENIRk(nbVilles,y,arr_probVISITER,arr_probVENIRk);

            //recalculing the nbCONATCT forcing
            calculateNbCONTACTForcing(nbVilles,arr_nbCONTACT0,arr_nbCONTACT1,periDISE,t,phi,arr_nbCONTACTFORCING);

/*
           char line[255];
           FILE * pFile;
           pFile = fopen ("beta.txt","a+");

           sprintf(line, "%.5f \t\t %.5f \t\t %.5f \n",t,arr_nbCONTACTFORCING[0][0],arr_nbCONTACTFORCING[1][1]);
           fputs (line,pFile);

           fclose (pFile);
   */
           //Calculating Value of \lamda at time t
           //calculerLamda(nbVilles,arr_xi,beta_sinusoidal,arr_rho,y,lamda);
           calLamdaNewInfec(nbVilles,arr_probVISITER,y,arr_nbCONTACTFORCING,arr_probINFECTER,arr_probVENIRk,lamda);

           //Calculating Values of the propensity functions
           calculerF(nbVilles,mu,sigma,gamma,lamda,y,f);

           //Generating the uniform random numbers
           if(typeRNG==FAST) {// using the method of New Infec
               p1 = ran1(idum,GENERATE, (char*)"a");//probability for time
               p2 = ran1(idum,GENERATE, (char*)"a");//probability for time
               p3 = ran1(idum,GENERATE, (char*)"a");//probability to choose event
           }
           else
               if(typeRNG==GOOD){// using the method of C++
                   p1 = rand()/(double)RAND_MAX;//probability for time
                   p2 = rand()/(double)RAND_MAX;//probability for time
                   p3 = rand()/(double)RAND_MAX;//probability to choose event
               }

           //Calculating the sum of all propensity probabilities of all cities
           fsumVilles = 0.0;
           for(int i = 0; i<nbVilles; i++){
               for(int j = 0; j<nevent; j++){
                   fsumVilles += f[i][j];
               }
           }

           //Calculating the time interval \tstep and updating time
           if(fsumVilles > 0.0) {
               tstep = (double)(-log(p1)/fsumVilles);
           }
           else {
               cout<<"no event"<<endl;
               cout<<"parce que: fsumVilles = "<<fsumVilles<<endl;
               break;
           }

           //Choosing randomly a city where event can occur
           nextville = choisirVille(nbVilles,nevent,p2,f);

           // Choosing randomly one event in the city chosen
           if(nextville != -1){
               event = choisirEvenement(nevent,nextville,p3,f);
           }
           else{
               printf ("(nextville == -1)fsumVilles = %.20f \n", fsumVilles);
               printf ("p2 = %.5f \n", p2);
               printf ("xpro = %.20f \n", xpro);
               printf ("fville = %.20f \n", fville);
                for(int i=0; i<nbVilles; i++){
                    cout<<"ville = "<<i<<endl;
                    for(int j=0; j<nevent; j++)
                        cout<<"f"<<j<<" = "<<f[i][j]<<"  ";
                    cout<<endl;
                }
           }

           //Doing the event \m
           if(event != -1){
               fairEvenementM(event,sigma,nextville,y,matrCumulI);
               // Check the local extinction and the stop of the local extinction
		//
	//for(int vil=0; vil<nbVilles; vil++){
	//cout<<"vil=="<<vil<<" S= "<<y[vil][iS]<<"  E="<<y[vil][iE]<<"  I="<<y[vil][iI] <<"    R="<<y[vil][iR] <<"   N="<<y[vil][iN]<<endl;}

                for(int i= 0; i< nbVilles; i++){
                    if((y[i][iE]==0) && (y[i][iI]==0) && (vecFlagExtLocal[i]==false)){
	//cout<<"khi co extinciton : <<"<<endl;
	//cout<<" S= "<<y[i][iS]<<"  E="<<y[i][iE]<<"  I="<<y[i][iI] <<"    R="<<y[i][iR] <<"   N="<<y[i][iN]<<endl;
                        tableExtLocal[i].push_back(t);
                        vecFlagExtLocal[i]=true;
                    }
                    if(((y[i][iE]!=0) || (y[i][iI]!=0)) && (vecFlagExtLocal[i]==true)){
	//cout<<"khi co STOP extinciton : <<"<<endl;
	//cout<<" S= "<<y[i][iS]<<"  E="<<y[i][iE]<<"  I="<<y[i][iI] <<"    R="<<y[i][iR] <<"   N="<<y[i][iN]<<endl;
                        tableRecolLocal[i].push_back(t);
                        vecFlagExtLocal[i]=false;
                    }
                }

               //Saving the values of variables at time t in a 3D table.
               if(t> pointTemps){
                  for(int i= 0; i< nbVilles; i++){
                       tpmSEIRNP.clear();
                       tpmSEIRNP.push_back(pointTemps);tpmSEIRNP.push_back(y[i][iS]); tpmSEIRNP.push_back(y[i][iE]); tpmSEIRNP.push_back(y[i][iI]);tpmSEIRNP.push_back(y[i][iR]);tpmSEIRNP.push_back(y[i][iN]);tpmSEIRNP.push_back(matrCumulI[i]);
                       tableResVilles[i].push_back(tpmSEIRNP);
                       matrCumulI[i] = 0.0;
                   }
                   pointTemps +=unitTIME;
               }
           }
           else{
               printf ("(event == -1)fville = %.20f \n", fsumVilles);
               printf ("p2 = %.5f \n", p3);
               printf ("xpro = %.20f \n", xpro);
               printf ("fprop = %.20f \n", fprop);
                for(int i=0; i<nbVilles; i++){
                    cout<<"ville = "<<i<<endl;
                    for(int j=0; j<nevent; j++)
                        cout<<"f"<<j<<" = "<<f[i][j]<<"  ";
                    cout<<endl;
                }
           }
           // update time t
             t = t + tstep;
          //continue
       }

       end = clock();
       // CPU time is used:
       double totalDure = (double)(end - start)/CLOCKS_PER_SEC;
       //Freeing the memory of the pointers
       delete []lamda;// 1D pointer
       delete []matrCumulI;// 1D pointer
       // Freeing the 2D pointers
       for (int i = 0; i <  nbVilles; i++)
       {
           delete []y[i];
           delete []f[i];
           delete []arr_nbCONTACTFORCING[i];
       }
       delete []y; delete[]f; delete arr_nbCONTACTFORCING;
       // result
       resMETAPOP listMetaPOP;
       listMetaPOP.tableResVilles = tableResVilles;
       listMetaPOP.tableExtLocal = tableExtLocal;
       listMetaPOP.tableRecolLocal = tableRecolLocal;

       return(listMetaPOP);
   }

     /***** The end of the main program******/

     /*****************************************************************************/
