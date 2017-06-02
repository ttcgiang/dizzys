/* 
$Id: seir_adaptivetau.cpp$
Created by TRAN Thi Cam Giang

--------------------------------------------------------------------------
Here, we model the SIR and SEIR models. This is a C++ implementation of the "adaptive tau-leaping"
 algorithm described by Cao Y, Gillespie DT, Petzold LR. The Journal of Chemical Physics (2007).
This programme interchanges C++ and R.
If building library outside of R package (i.e. for debugging):
        R CMD SHLIB adaptivetau.cpp
 --------------------------------------------------------------------------
*/
//libraries of R
#include <iostream>
#include <vector>
#include <limits>
#include <sstream>
#include <stdexcept>
#include "stoSEIRNewInfec.h"
#include "index.h"

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>
#include "Rwrappers.h"


using namespace std;


enum EStepType {
    eExact = 0,
    eExplicit,
    eImplicit
};

const bool debug = false;

//use below rather than R's "error" directly (which will not free memory, etc.)
#ifdef throwError
#undef throwError
#endif
#define throwError(e) { ostringstream s; s << e; throw runtime_error(s.str()); }

// Call R function without responding to user interrupts (or other, I
// suppose).  This lets use catch them and free memory in the heap.
// Unfortunately, R_interrupts_suspended is not directly part of the R
// API (although it is used by the BEGIN/END_SUSPEND_INTERRUPTS macros
// in R_ext/GraphicsDevice.h)
extern "C" { LibExtern Rboolean R_interrupts_suspended;}

SEXP evalWithoutInterrupts(SEXP expr) {
    //R_interrupts_suspended = TRUE;
    SEXP res = eval(expr, NULL);
   // R_interrupts_suspended = FALSE;
    return res;
}

/************************************************************/
// CODE of TRAN THI CAM GIANG that are added in the c++ codes of the package "adaptive tau"
extern "C" { 
/****/
/*
    Goal: calculating the propensity functions for all subpopulation.
        Used for SIR model
    */  
double* seirratefunc(double* x,int nbVilles,double** arr_probVISITER,double** arr_probVENIRk,double** arr_probINFECTER,
		double** arr_nbCONTACT0,double** arr_nbCONTACT1,
        double periDISE, double m_T,double *phi,double** arr_nbCONTACTFORCING,double mu,double sigma, double gamma){
    // 4 variables (S, I, R, N) and 6 events
    int nbVarSEIRN=6, nevent=8;
    //Initializing
    double *lamda = new double[nbVilles];
    for(int i = 0; i< nbVilles ; i++) lamda[i]=0.0;
    //Converting 1D vector to 2D vector **y
    unsigned long **y;
    y = new unsigned long  *[nbVilles];
    for (int i = 0; i<nbVilles; i++)
    {
        y[i] = new unsigned long [nbVarSEIRN];
        y[i][iS] = (unsigned long)x[0*nbVilles+i];//iS
        y[i][iE] = (unsigned long)x[1*nbVilles+i];//iE
        y[i][iI] = (unsigned long)x[2*nbVilles+i];//iI
        y[i][iR] = (unsigned long)x[3*nbVilles+i];//iR
        y[i][iN] = y[i][iS]+y[i][iE]+y[i][iI]+y[i][iR];//iN

    }
    // recalculating the matrix "COME from the city k
    calculerProbVENIRk(nbVilles,y,arr_probVISITER,arr_probVENIRk);
    //recalculing the nbCONATCT forcing
    /*
    //calculateNbCONTACTForcing(nbVilles,arr_nbCONTACT0,arr_nbCONTACT1,periDISE,m_T,phi,arr_nbCONTACTFORCING);
    cout<<"DDDD  "<<m_T<<"  "<<arr_nbCONTACTFORCING[0][0]<<endl;
    //Calculating Value of \lamda at time t
    for (int i = 0; i<nbVilles; i++){
        for(int j=0;j<nbVarSEIRN;j++)
            cout<<"yyyy"<<i<<j<<"  =="<<y[i][j]<<endl;
    }
    */
    
    calLamdaNewInfec(nbVilles,arr_probVISITER,y,arr_nbCONTACTFORCING,arr_probINFECTER,arr_probVENIRk,lamda);
    //Calculating the values of all propensity function
    double **f;
    f = new double*[nbVilles];
    for (int i = 0; i <  nbVilles; i++)
    {
        f[i] = new double[nevent];
        for(int j = 0; j< nevent ; j++) f[i][j] = 0.0;
    }

    for(int i=0; i< nbVilles;i++){
        //S born
        f[i][0] = mu*y[i][iN];
        //S die
        f[i][1] = mu*y[i][iS];
        //E die
        f[i][2] = mu*y[i][iE];
        //I die
        f[i][3] = mu*y[i][iI];
        //R die
        f[i][4] = mu*y[i][iR];
        //S-> E
        f[i][5] = lamda[i]*y[i][iS];
        //E-> I
        f[i][6] = sigma*y[i][iS];
        //I->R
        f[i][7] = gamma*y[i][iI];
    }
    // vector of result
    vector<double> vecRes;
    for(int j = 0; j< nevent ; j++)
    {
        for (int i = 0; i <  nbVilles; i++)
            vecRes.push_back(f[i][j]);
    }
    int nbvecRes=vecRes.size();
    double *resRates = new double[nbvecRes];
    for(int i=0; i<vecRes.size();i++){
        resRates[i] = vecRes[i];
    }
    // Freeing the memory of all dynnamic pointers
    delete []lamda;
    for (int i = 0; i<nbVilles; i++)
    {
        delete []y[i];
        delete []f[i];
    }
    delete []y; delete[]f;
    return(resRates);
}
/****/
    
  /*
    Goal: calculating the propensity functions for all subpopulation.
        Used for SIR model
    */
double* sirratefunc(double* x,int nbVilles,double** arr_probVISITER,double** arr_probVENIRk,double** arr_probINFECTER,
		double** arr_nbCONTACT0,double** arr_nbCONTACT1,
        double periDISE, double m_T,double *phi,double** arr_nbCONTACTFORCING, double mu,double gamma){
    // 4 variables (S, I, R, N) and 6 events
    int nbVarSIRN=5, nevent=6;
    //Initializing
    double *lamda = new double[nbVilles];
    //Converting 1D vector to 2D vector **y
    unsigned long **y;
    y = new unsigned long  *[nbVilles];
    for (int i = 0; i<nbVilles; i++)
    {
        y[i] = new unsigned long [nbVarSIRN];
        y[i][iS] = (unsigned long)x[0*nbVilles+i];//iS
        y[i][iI] = (unsigned long)x[1*nbVilles+i];//iI
        y[i][iR] = (unsigned long)x[2*nbVilles+i];//iR
        y[i][iN] = y[i][iS]+y[i][iI]+y[i][iR];//iN

    }
    // recalculating the matrix "COME from the city k
    calculerProbVENIRk(nbVilles,y,arr_probVISITER,arr_probVENIRk);
    //recalculing the nbCONATCT forcing
    //calculateNbCONTACTForcing(nbVilles,arr_nbCONTACT0,arr_nbCONTACT1,periDISE,m_T,phi,arr_nbCONTACTFORCING);
    //Calculating Value of \lamda at time t
    calLamdaNewInfec(nbVilles,arr_probVISITER,y,arr_nbCONTACTFORCING,arr_probINFECTER,arr_probVENIRk,lamda);

    //Calculating the values of all propensity functions
    double **f;
    f = new double*[nbVilles];
    for (int i = 0; i <  nbVilles; i++)
    {
        f[i] = new double[nevent];
        for(int j = 0; j< nevent ; j++) f[i][j] = 0.0;
    }
    for(int i=0; i< nbVilles;i++){
        //S born
        f[i][0] = mu*y[i][iN];
        //S die
        f[i][1] = mu*y[i][iS];
        //I die
        f[i][2] = mu*y[i][iI];
        //R die
        f[i][3] = mu*y[i][iR];
        //S-> I
        f[i][4] = lamda[i]*y[i][iS];
        //I->R
        f[i][5] = gamma*y[i][iI];
    }
    // vector of result
    vector<double> vecRes;
    for(int j = 0; j< nevent ; j++)
    {
        for (int i = 0; i <  nbVilles; i++)
            vecRes.push_back(f[i][j]);
    }
    int nbvecRes=vecRes.size();
    double *resRates = new double[nbvecRes];
    for(int i=0; i<vecRes.size();i++){
        resRates[i] = vecRes[i];
    }
    // Freeing the memory of all dynnamic pointers
    delete []lamda;
    for (int i = 0; i<nbVilles; i++)
    {
        delete []y[i];
        delete []f[i];
    }
    delete []y; delete[]f;
    return(resRates);
}
/****/
}

/************** code TRAN Thi Cam Giang*****************/
///////////////////////////////////////////////////
class CStochasticEqns {
public:
    CStochasticEqns(SEXP initVal, int *nu, unsigned int numTrans,
                    //
                    SEXP rateFunc, SEXP rateJacobianFunc,
                    // for model SEIR
                    int nbVilles,double sigma,double gamma,double mu,
                    // for metapopulation
                    double** arr_probVISITER,double** arr_nbCONTACT0,double** arr_nbCONTACT1, double** arr_probINFECTER, double** arr_probVENIRk,
                    double *phiPHASE,double periDISE,
                    //
                    double* changeBound, SEXP maxTauFunc, SEXP detTrans) {
        // copy initial values into new vector (keeping in SEXP vector
        // allows easy calling of R function to calculate rates)
        m_NumStates = length(initVal);
        SEXP x; double *x_copy;
        PROTECT(x = NEW_NUMERIC(m_NumStates));
        x_copy = NUMERIC_POINTER(x);
        double *rinitVal = REAL(initVal);
        double* cinitVal; create1D(m_NumStates,0.0,cinitVal);
        for(int i=0; i<m_NumStates; i++) {cinitVal[i] = rinitVal[i];x_copy[i]=cinitVal[i];}
        //copyVector(x, initVal);
        //copyVector(x, coerceVector(initVal,REALSXP));

        if (isNull(getAttrib(initVal, R_NamesSymbol))) {
            m_VarNames = NULL;
        } else {
            SEXP initNames = getAttrib(initVal, R_NamesSymbol);
            int m_NumNames=length(initNames);
            PROTECT(m_VarNames = allocVector(STRSXP,m_NumNames));
            for(int i=0; i<m_NumNames; i++) {
                char* cnames=(char *) CHAR(STRING_ELT(initNames,i));
                SET_STRING_ELT(m_VarNames, i,mkChar(cnames));
            }
            //copyVector(m_VarNames, initNames);
            setAttrib(x, R_NamesSymbol,m_VarNames);
        }
        m_X = REAL(x);
        // copy full-size Nu matrix into sparse matrix
        m_Nu.resize(numTrans);
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            for (unsigned int j = 0;  j < numTrans;  ++j) {
                if (nu[j*m_NumStates + i] != 0) {
                    m_Nu[j].push_back(SChange(i, nu[j*m_NumStates + i]));
                }
            }
        }
        m_TransCats.resize(numTrans, eNoncritical);
        // potentially flag some transitions as "deterministic"
        if (detTrans  &&  !isNull(detTrans)) {
            x_SetDeterministic(LOGICAL(detTrans), length(detTrans));
        }
        // needed for ITL
        x_IdentifyBalancedPairs();
        x_IdentifyRealValuedVariables();
        // prepare R function for evaluation by setting up arguments
        // (current X values, parameters, current time)
        SEXP s_time;
        PROTECT(s_time = allocVector(REALSXP, 1));
        m_T = REAL(s_time);
        // TTCGIang: for SEIR and SIR models
        GnbVilles=nbVilles;  // number of subpopulations
        //Model
        Gsigma=sigma; Ggamma=gamma;Gmu=mu; GperiDISE=periDISE;
        // for simulation
        // for metapopulation
        //1D
        GphiPHASE=new double[GnbVilles];
        for (int i = 0; i < GnbVilles; i++) GphiPHASE[i]=phiPHASE[i];

        //2D
        //Garr_nbCONTACTFORCING
        Garr_probVISITER = new double*[GnbVilles];
        Garr_nbCONTACT0 = new double*[GnbVilles];
        Garr_nbCONTACT1 = new double*[GnbVilles];
        Garr_probINFECTER = new double*[GnbVilles];
        Garr_probVENIRk = new double*[GnbVilles];
        Garr_nbCONTACTFORCING=new double*[GnbVilles];
        for (int i = 0; i < GnbVilles; i++)
        {
            Garr_probVISITER[i] = new double[GnbVilles];
            Garr_nbCONTACT0[i] = new double[GnbVilles];
            Garr_nbCONTACT1[i] = new double[GnbVilles];
            Garr_probINFECTER[i] = new double[GnbVilles];
            Garr_probVENIRk[i] = new double[GnbVilles];
            Garr_nbCONTACTFORCING[i]=new double[GnbVilles];
            for(int j = 0; j< GnbVilles ; j++){
                Garr_probVISITER[i][j] = arr_probVISITER[i][j];
                Garr_nbCONTACT0[i][j] = arr_nbCONTACT0[i][j];
                Garr_nbCONTACT1[i][j] = arr_nbCONTACT1[i][j];
                Garr_probINFECTER[i][j] = arr_probINFECTER[i][j];
                Garr_probVENIRk[i][j] = arr_probVENIRk[i][j];
                Garr_nbCONTACTFORCING[i][j]=0.0;
            }
        }
        //
        calculateNbCONTACTForcing(GnbVilles,Garr_nbCONTACT0,Garr_nbCONTACT1,GperiDISE,0,GphiPHASE,Garr_nbCONTACTFORCING);
        //
        PROTECT(m_RateFunc);// = lang4(rateFunc, x, params, s_time));
        if (!rateJacobianFunc  ||  isNull(rateJacobianFunc)) {
            m_RateJacobianFunc = NULL;
        }/*
        else {
          PROTECT(m_RateJacobianFunc = lang4(rateJacobianFunc, x,params, s_time));
        }*/
        //Initializing the values for the flags of the variables SEIR
        flagSEIR = TRUE;
        if(Gsigma==INFINITY) flagSEIR=FALSE;
        //
        m_Rates = NULL;
        //default parameters to adaptive tau leaping algorithm
        m_Epsilon = 0.01;
        m_Ncritical = 10;
        m_Nstiff = 100;
        m_ExactThreshold = 10;
        m_Delta = 0.01;
        m_NumExactSteps[eExact] = 100;
        m_NumExactSteps[eExplicit] = 100;
        m_NumExactSteps[eImplicit] = 10;
        m_ITLConvergenceTol = 0.01;
        m_MaxTau = numeric_limits<double>::infinity();
        m_MaxSteps = 0; // special case 0 == no limit

        //useful additional parameters
        m_ExtraChecks = true;
        m_RateChangeBound = changeBound;
        if (!maxTauFunc  ||  isNull(maxTauFunc)) {
            m_MaxTauFunc = NULL;
        }/*
        else {
            PROTECT(m_MaxTauFunc = lang4(maxTauFunc, x, params,s_time));
        }*/


        //check initial conditions to make sure legit
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            if (m_X[i] < 0) {
                throwError("initial value for variable " << i+1 <<
                           " must be positive (currently " << m_X[i] << ")");
            }
            if (!m_RealValuedVariables[i]  &&  (m_X[i] - trunc(m_X[i]) > 1e5)) {
                throwError("initial value for variable " << i+1 <<
                           " must be an integer (currently " << m_X[i] << ")");
            }
        }
        
        *m_T = 0;
        m_PrevStepType = eExact;
        GetRNGstate();
        //Free pointers
        //delete []x_copy;
        //delete []cinitVal;
    }
    /********************************/

    ~CStochasticEqns(void) {
        int cnt = 3;
        if (m_RateJacobianFunc != NULL) {
            ++cnt;
        }
        if (m_Rates != NULL) {
            ++cnt;
        }
        if (m_MaxTauFunc != NULL) {
            ++cnt;
        }
        if (m_VarNames != NULL) {
            ++cnt;
        }
        UNPROTECT(cnt);
        //delete []m_Rates;
        //1D
        delete []GphiPHASE;
        //2D
        for (int i = 0; i<GnbVilles; i++){
            delete [] Garr_probVISITER[i];
            delete [] Garr_nbCONTACT0[i];
            delete [] Garr_nbCONTACT1[i];
            delete [] Garr_probINFECTER[i];
            delete [] Garr_probVENIRk[i];
            delete [] Garr_nbCONTACTFORCING[i];
        }
        delete  [] Garr_probVISITER; delete [] Garr_nbCONTACT0;
        delete [] Garr_nbCONTACT1; delete [] Garr_probINFECTER;
        delete [] Garr_probVENIRk; delete [] Garr_nbCONTACTFORCING;


    }
    void SetTLParams(SEXP list) {
        SEXP names = getAttrib(list, R_NamesSymbol);

        for (int i = 0;  i < length(names);  ++i) {
            if (strcmp("epsilon", CHAR(STRING_PTR(names)[i])) == 0) {
                if (!isReal(VECTOR_ELT(list, i))  ||
                    length(VECTOR_ELT(list, i)) != 1) {
                    throwError("invalid value for parameter '" <<
                               CHAR(STRING_PTR(names)[i]) << "'");
                }
                m_Epsilon = REAL(VECTOR_ELT(list, i))[0];
            } else if (strcmp("delta", CHAR(STRING_PTR(names)[i])) == 0) {
                if (!isReal(VECTOR_ELT(list, i))  ||
                    length(VECTOR_ELT(list, i)) != 1) {
                    throwError("invalid value for parameter '" <<
                               CHAR(STRING_PTR(names)[i]) << "'");
                }
                m_Delta = REAL(VECTOR_ELT(list, i))[0];
            } else if (strcmp("maxtau", CHAR(STRING_PTR(names)[i])) == 0) {
                if (!isReal(VECTOR_ELT(list, i))  ||
                    length(VECTOR_ELT(list, i)) != 1) {
                    throwError("invalid value for parameter '" <<
                               CHAR(STRING_PTR(names)[i]) << "'");
                }
                m_MaxTau = REAL(VECTOR_ELT(list, i))[0];
            } else if (strcmp("extraChecks",
                              CHAR(STRING_PTR(names)[i])) == 0) {
                if (!isLogical(VECTOR_ELT(list, i))  ||
                    length(VECTOR_ELT(list, i)) != 1) {
                    throwError("invalid value for parameter '" <<
                               CHAR(STRING_PTR(names)[i]) << "'");
                }
                m_ExtraChecks = LOGICAL(VECTOR_ELT(list, i))[0];
            } else {
                warning("ignoring unknown parameter '%s'",
                        CHAR(STRING_PTR(names)[i]));
            }
        }
    }

    // Functions below are a hack suggested by Simon Urbanek (although
    // he "would not recommend for general use") to check if the user has
    // asked to interrupt execution.  The problem with calling
    // R_CheckUserInterrupt directly is it longjmps and we don't have
    // a chance to free memory from the heap.
    // http://tolstoy.newcastle.edu.au/R/e13/devel/11/04/1049.html
    static void chkIntFn(void*) { R_CheckUserInterrupt(); }
    bool checkUserInterrupt(void) {
        return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
    }

    void EvaluateATLUntil(double tF) {
        unsigned int c = 0;
        //add initial conditions to time series
        m_TimeSeries.push_back(STimePoint(0, m_X, m_NumStates));
        //main loop
        while (*m_T < tF  &&  (m_MaxSteps == 0 || c < m_MaxSteps)) {
            x_SingleStepATL(tF);
            if (++c % 10 == 0  &&  checkUserInterrupt()) {
                throwError("simulation interrupted by user at time " << *m_T
                           << " after " << c << " time steps.");
            }
        }
        //save RNG state back to R (could also do in destructor, but
        //no harm in extra calls to PutRNGstate and avoids potential
        //problems with PROTECTing return value from GetTimeSeriesSEXP
        PutRNGstate();
    }
    void EvaluateExactUntil(double tF) {
        unsigned int c = 0;
        //add initial conditions to time series
        m_TimeSeries.push_back(STimePoint(0, m_X, m_NumStates));
        //main loop
        while (*m_T < tF  &&  (m_MaxSteps == 0 || c < m_MaxSteps)) {
            x_UpdateRates();
            x_SingleStepExact(tF);
            m_TimeSeries.push_back(STimePoint(*m_T, m_X, m_NumStates));
            if (++c % 10 == 0  &&  checkUserInterrupt()) {
                throwError("simulation interrupted by user at time " << *m_T
                           << " after " << c << " time steps.");
            }
        }
        //save RNG state back to R (could also do in destructor, but
        //no harm in extra calls to PutRNGstate and avoids potential
        //problems with PROTECTing return value from GetTimeSeriesSEXP
        PutRNGstate();
    }

    SEXP GetTimeSeriesSEXP(void) const {
        SEXP res;
        PROTECT(res = allocMatrix(REALSXP, m_TimeSeries.size(), m_NumStates+1));
        double *rvals = REAL(res);
        for (unsigned int t = 0;  t < m_TimeSeries.size();  ++t) {
            rvals[t] = m_TimeSeries[t].m_T;
            for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                rvals[(i+1) * m_TimeSeries.size() + t] = m_TimeSeries[t].m_X[i];
            }
        }

        SEXP dimnames, colnames;
        PROTECT(dimnames = allocVector(VECSXP, 2));
        PROTECT(colnames = allocVector(VECSXP, m_NumStates+1));
        SET_VECTOR_ELT(dimnames, 1, colnames);
        SET_VECTOR_ELT(colnames, 0, mkChar("time"));
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            if (m_VarNames  &&  (unsigned int)length(m_VarNames) > i) {
                SET_VECTOR_ELT(colnames, i+1,
                               STRING_PTR(m_VarNames)[i]);
            } else {
                char name[10];
                snprintf(name, 10, "x%i", i+1);
                SET_VECTOR_ELT(colnames, i+1, mkChar(name));
            }
        }
        setAttrib(res, R_DimNamesSymbol, dimnames);

        UNPROTECT(3);
        return res;
    }

protected:
    enum ETransCat {
        eCritical,
        eNoncritical,
        eDeterministic
    };
    typedef vector<ETransCat> TTransCats;
    typedef vector<pair<unsigned int, unsigned int> > TBalancedPairs;
    typedef vector<bool> TBools;
    typedef double* TStates;
    typedef double* TRates;
    struct SChange {
        SChange(short int s, short int m) : m_State(s), m_Mag(m) {}
        short int m_State;
        short int m_Mag;
    };
    typedef vector< vector<SChange> > TTransitions;
    struct STimePoint {
        STimePoint(double t, double *x, int n) {
            m_T = t;
            m_X = new double[n];
            memcpy(m_X, x, n*sizeof(*x));
        }
        double m_T;
        double *m_X;
    };
    class CTimeSeries : public vector<STimePoint> {
    public:
        ~CTimeSeries(void) {
            for (iterator i = begin();  i != end();  ++i) {
                delete[] i->m_X; i->m_X = NULL;
            }
        }
    };

protected:
    void x_IdentifyBalancedPairs(void);
    void x_IdentifyRealValuedVariables(void);
    void x_SetDeterministic(int *det, unsigned int n);

    void x_AdvanceDeterministic(double deltaT, bool clamp = false);
    void x_SingleStepExact(double tf);
    bool x_SingleStepETL(double tau);
    bool x_SingleStepITL(double tau);
    void x_SingleStepATL(double tf);

    void x_UpdateRates(void) {
        if (m_ExtraChecks) {
            for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                if (m_X[i] < 0) {
                    cout<<"i="<<i<<endl;
                    cout<<"m_X[i]==="<<m_X[i]<<endl;
                    throwError("negative variable: " << i+1 << " is " <<
                               m_X[i] << " (check rate function "
                               "and/or transition matrix)");
                } else if (isnan(m_X[i])) {
                    throwError("NaN variable: " << i+1 << " is " <<
                               m_X[i] << " (check rate function "
                               "and/or transition matrix)");
                }
            }
        }

        //not sure if this protect/unprotect block is necessary, but
        //seems better to err on the safe side
        if (m_Rates != NULL) { 
            UNPROTECT(1);
        }
        SEXP res;// = evalWithoutInterrupts(m_RateFunc);
        PROTECT(res);
        //TTCGiang for the SEIR and SIR models
        // Here, we call the function 'seirratefunc' that calculates the propensity functions
        if(flagSEIR){//SEIR model
            /*
            FILE *pFile; char line[255] = "";
            pFile = fopen("forcing.txt","a+");
            for(int i=0;i<GnbVilles;i++){
                for(int j=0;j<GnbVilles;j++){
                    cout<<"  "<<*m_T<<"  "<<Garr_nbCONTACTFORCING[i][j]<<endl;
                    sprintf(line, "%.3f\t\t%.3f \n",*m_T ,Garr_nbCONTACTFORCING[i][j]);
                }
                fputs (line,pFile);
                fclose (pFile);
            }*/
            calculateNbCONTACTForcing(GnbVilles,Garr_nbCONTACT0,Garr_nbCONTACT1,GperiDISE,*m_T,GphiPHASE,Garr_nbCONTACTFORCING);
            m_Rates  = seirratefunc(m_X,GnbVilles,Garr_probVISITER,Garr_probVENIRk,Garr_probINFECTER,
        Garr_nbCONTACT0,Garr_nbCONTACT1,GperiDISE,*m_T,GphiPHASE,Garr_nbCONTACTFORCING,Gmu,Gsigma,Ggamma);
        }
        else{//SIR model
             calculateNbCONTACTForcing(GnbVilles,Garr_nbCONTACT0,Garr_nbCONTACT1,GperiDISE,*m_T,GphiPHASE,Garr_nbCONTACTFORCING);
            m_Rates = sirratefunc(m_X,GnbVilles,Garr_probVISITER,Garr_probVENIRk,Garr_probINFECTER,
		Garr_nbCONTACT0,Garr_nbCONTACT1,GperiDISE,*m_T,GphiPHASE,Garr_nbCONTACTFORCING,Gmu,Ggamma);
          //  for(int i=0; i<6; i++) cout <<"ZZZZ = "<<m_Rates[i]<<endl;
            //cout<<"m_T=="<<*m_T<<endl;
        }

        if (m_ExtraChecks) {
            for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
                if (isnan(m_Rates[j])) {
                    cout<<"j ="<<j<<endl;
                    cout<<"m_Rates[j]==="<<m_Rates[j]<<endl;
                    throwError("invalid rate function -- rate for transition "
                               << j+1 << "is not a number (NA/NaN)! (check "
                               "for divison by zero or similar)");
                }
                if (m_Rates[j] < 0) {
                    throwError("invalid rate function -- rate for transition "
                               << j+1 << "is negative!");
                }
            }
        }
    }
    /********************************/
    double* x_CalcJacobian(void) {
        SEXP res = evalWithoutInterrupts(m_RateJacobianFunc);
        if (!isMatrix(res)) {
            throwError("invalid Jacobian function -- should return a " <<
                       m_NumStates << " by " << m_Nu.size() << " matrix");
        }
        unsigned int nrow = INTEGER(getAttrib(res, R_DimSymbol))[0];
        unsigned int ncol = INTEGER(getAttrib(res, R_DimSymbol))[1];
        if (nrow != m_NumStates  ||  ncol != m_Nu.size()) {
            throwError ("invalid Jacobian function -- returned a " << nrow
                        << " by " << ncol << " matrix instead of the expected "
                        << m_NumStates << " by " << m_Nu.size() <<
                        " (variables by transitions)");
        }
        return REAL(res);
    }
      /********************************/
    double x_CalcUserMaxTau(void) {
        if (!m_MaxTauFunc) { throwError("logic error at line " << __LINE__) }
        SEXP res = evalWithoutInterrupts(m_MaxTauFunc);
        if (length(res) != 1  || !isReal(res)) {
            throwError("invalid return value from maxTau function (should be "
                       "a single real number)");
        }
        return REAL(res)[0];
    }
     /********************************/
    unsigned int x_PickCritical(double prCrit) const;

    double x_TauEx(void) const {
        double tau = numeric_limits<double>::infinity();
        vector <double> mu(m_NumStates, 0);
        vector <double> sigma(m_NumStates, 0);

        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] != eCritical) {
                for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                    const SChange &c = m_Nu[j][i];
                    mu[c.m_State] += c.m_Mag * m_Rates[j];
                    sigma[c.m_State] += c.m_Mag * c.m_Mag * m_Rates[j];
                }
            }
        }
        //cerr << "-=| mu:";
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            //cerr << "\t" << mu[i];
            double val = max(m_Epsilon * m_X[i] / m_RateChangeBound[i],
                             1.)/fabs(mu[i]);
            //cerr << "/" << val;
            if (val < tau) {
                tau = val;
                if (tau < 0) {
                    throwError("tried to select tau < 0; most likely means "
                               "your rate function generated a negative rate");
                }
            }
            val = pow(max(m_Epsilon * m_X[i] / m_RateChangeBound[i],
                          1.),2) / sigma[i];
            if (val < tau) {
                tau = val;
                if (tau < 0) {
                    throwError("tried to select tau < 0; most likely means "
                               "your rate function generated a negative rate");
                }
            }
        }
        //cerr << endl;

        return tau;
    }
 /********************************/
    double x_TauIm(void) const {
        if (!m_RateJacobianFunc) {
            return 0;
        }
        vector<bool> equil(m_TransCats.size(), false);
        for (TBalancedPairs::const_iterator i = m_BalancedPairs.begin();
             i != m_BalancedPairs.end();  ++i) {
            if (fabs(m_Rates[i->first] - m_Rates[i->second]) <=
                m_Delta * min(m_Rates[i->first], m_Rates[i->second])) {
                equil[i->first] = true;
                equil[i->second] = true;
            }
        }

        vector<double> mu(m_NumStates, 0);
        vector<double> sigma(m_NumStates, 0);
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] != eCritical  &&  !equil[j]) {
                for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                    const SChange &c = m_Nu[j][i];
                    mu[c.m_State] += c.m_Mag * m_Rates[j];
                    sigma[c.m_State] += c.m_Mag * c.m_Mag * m_Rates[j];
                }
            }
        }

        double tau = numeric_limits<double>::infinity();
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            double val = max(m_Epsilon * m_X[i] / m_RateChangeBound[i],
                             1.)/fabs(mu[i]);
            if (val < tau) {
                tau = val;
            }
            val = pow(max(m_Epsilon * m_X[i] / m_RateChangeBound[i],
                          1.),2) / sigma[i];
            if (val < tau) {
                tau = val;
            }
        }
        
        return tau;
    }
 /********************************/
private:
    bool m_ExtraChecks; //turns on extra checks on rates returned by
                        //user-supplied rate function. Slower, but if
                        //the rate function does have a bug, this will
                        //give a more meaningful error message.

    unsigned int m_Ncritical;
    double m_Nstiff;
    double m_Epsilon;
    double m_ExactThreshold;
    double m_Delta;
    unsigned int m_NumExactSteps[3];
    double m_ITLConvergenceTol;
    double m_MaxTau;
    unsigned int m_MaxSteps;

    TRates m_Rates; // *current* rates (must be updated if m_X changes!)   
    bool flagSEIR;  // current flag for the SIER model
    double *m_T;    // *current* time
    int GnbVilles;  // number of subpopulations
    double GT;
    //Model
    double Gsigma, Ggamma, Gmu;
    double GperiDISE;
    double ** Garr_nbCONTACTFORCING;
    // for simulation
    // for metapopulation
    double** Garr_probVISITER;
    double** Garr_nbCONTACT0;
    double** Garr_nbCONTACT1;
    double** Garr_probINFECTER; 
    double** Garr_probVENIRk;
    double* GphiPHASE;
	
    TBalancedPairs m_BalancedPairs;
    TBools m_RealValuedVariables;
    EStepType m_PrevStepType;

    TStates m_X;              //current state
    unsigned int m_NumStates; //total number of states
    SEXP m_VarNames;          //variable names (if any)
    TTransitions m_Nu;        //state changes caused by transition
    TTransCats m_TransCats; //inc. whether transition is deterministic
    SEXP m_RateFunc; //R function to calculate rates as f(m_X)
    SEXP m_RateJacobianFunc; //R function to calculate Jacobian of rates as f(m_X) [optional!]
    double *m_RateChangeBound; //see Cao (2006) for details
    SEXP m_MaxTauFunc; //R function to calculate maximum leap given curr. state

    CTimeSeries m_TimeSeries;
};


/*---------------------------------------------------------------------------*/
// PRE : m_Nu initialized
// POST: all balanced pairs of transitions identified & saved
void CStochasticEqns::x_IdentifyBalancedPairs(void) {
    for (unsigned int j1 = 0;  j1 < m_Nu.size();  ++j1) {
        for (unsigned int j2 = j1 + 1;  j2 < m_Nu.size();  ++j2) {
            if (m_Nu[j1].size() != m_Nu[j2].size()) {
                continue;
            }
            unsigned int i;
            for (i = 0;  i < m_Nu[j1].size()  &&
                     m_Nu[j1][i].m_State == m_Nu[j2][i].m_State  &&
                     m_Nu[j1][i].m_Mag == -m_Nu[j2][i].m_Mag;  ++i);
            if (i == m_Nu[j1].size()) {
                m_BalancedPairs.push_back(TBalancedPairs::value_type(j1, j2));
                if (debug) {
                    cerr << "balanced pair " << j1 << " and " << j2 << endl;
                }
            }
        }
    }
}

/*---------------------------------------------------------------------------*/
// PRE : m_Nu initialized, boolean vector flagging transitions as
// deterministic or not (could be just FALSE and n=1)
// POST: deterministic flag set where appropriate
void CStochasticEqns::x_SetDeterministic(int *det, unsigned int n) {
    if (n != m_Nu.size()  &&  n > 1) {
        throwError("mismatch between length of logical vector specifying "
                   "deterministic transitions and total number of transitions");
    }
    bool atLeastOneStochastic = false;
    for (unsigned int i = 0;  i < n;  ++i) {
        if (det[i]) {
            m_TransCats[i] = eDeterministic;
        } else {
            atLeastOneStochastic = true;
        }
    }
    if (!atLeastOneStochastic) {
        throwError("At least one transition must be stochastic (all "
                   "transitions are currently flagged as deterministic).");
    }
    x_IdentifyRealValuedVariables();
}

/*---------------------------------------------------------------------------*/
// PRE : m_Nu initialized, deterministic transition set (if any)
// POST: all variables identified will take real values
// (i.e. either non-integer nu or modified by a deterministic transition)
void CStochasticEqns::x_IdentifyRealValuedVariables(void) {
    m_RealValuedVariables.clear();
    m_RealValuedVariables.resize(m_NumStates, false);

    for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
        if (m_TransCats[j] == eDeterministic) {
            for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                m_RealValuedVariables[m_Nu[j][i].m_State] = true;
            }
        } else {
            //code below not used since m_Nu matrix is forced to be integers
            for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                if (m_Nu[j][i].m_Mag - trunc(m_Nu[j][i].m_Mag) > 1e-5) {
                    m_RealValuedVariables[m_Nu[j][i].m_State] = true;
                }
            }
        }
    }
                    
}

/*---------------------------------------------------------------------------*/
// PRE : list of critical transitions & their total rate
// POST: one picked according to probability
unsigned int CStochasticEqns::x_PickCritical(double critRate) const {
    double r = runif(0,1);
    double d = 0;
    unsigned int j;
    for (j = 0;  j < m_Nu.size()  &&  d < r;  ++j) {
        if (m_TransCats[j] == eCritical) {
            d += m_Rates[j]/critRate;
        }
    }
    if (!(d >= r)) { throwError("logic error at line " << __LINE__) }
    return j-1;
}

/*---------------------------------------------------------------------------*/
// PRE : time period to step; whether to clamp variables at 0
// POST: all determinisitic transitions updated by the expected amount
// (i.e. Euler method); if clamping, then negative variables set to 0.
void CStochasticEqns::x_AdvanceDeterministic(double deltaT, bool clamp) {
    if (clamp) {
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] == eDeterministic) {
                for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                    m_X[m_Nu[j][i].m_State] += m_Nu[j][i].m_Mag * m_Rates[j] *
                        deltaT;
                    //clamp at zero if specified
                    if (m_X[m_Nu[j][i].m_State] < 0) {
                        m_X[m_Nu[j][i].m_State] = 0;
                    }
                }
            }
        }
    } else {
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] == eDeterministic) {
                for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                    m_X[m_Nu[j][i].m_State] += m_Nu[j][i].m_Mag * m_Rates[j] *
                        deltaT;
                }
            }
        }
    }
}

/*---------------------------------------------------------------------------*/
// PRE : simulation end time; **transition rates already updated**
// POST: a *single* transition taken (no approximation necessary)
void CStochasticEqns::x_SingleStepExact(double tf) {
    double stochRate = 0;
    double detRate = 0;
    for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
        if (m_TransCats[j] != eDeterministic) {
            stochRate += m_Rates[j];
        } else {
            detRate += m_Rates[j];
        }
    }
    if (stochRate == 0) {
        if (detRate == 0) {
            *m_T = numeric_limits<double>::infinity();
        } else {
            double tau = min(1/detRate, tf-*m_T);
            x_AdvanceDeterministic(tau, true);
            *m_T += tau;
        }
        return;
    }

    double tau = rexp(1./stochRate);
    if (tau > tf - *m_T) {
        tau = tf - *m_T;
    } else { // only take step if don't go off end
        double r = runif(0,1);
        double d = 0;
        unsigned int j;
        for (j = 0;  j < m_Nu.size()  &&  d < r;  ++j) {
            if (m_TransCats[j] != eDeterministic) {
                d += m_Rates[j]/stochRate;
            }
        }
        if (!(d >= r)) { throwError("logic error at line " << __LINE__) }
        --j;

        //take transition "j"
        for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
            m_X[m_Nu[j][i].m_State] += m_Nu[j][i].m_Mag;
        }
    }

    //clamp deterministic at 0, assuming that it is unreasonable to
    //take a smaller step then exact.
    x_AdvanceDeterministic(tau, true);
    *m_T += tau;
}

/*---------------------------------------------------------------------------*/
// PRE : tau value to use for step, list of "critical" transitions
// POST: whether single IMPLICIT tau step was successfully taken (m_X
// updated if so)
// NOTE: See equation (7) in Cao et al. (2007)
bool CStochasticEqns::x_SingleStepITL(double tau) {
    if (!m_RateJacobianFunc) { throwError("logic error at line " << __LINE__) }
    double *origX = new double[m_NumStates];
    double *origRates = new double[m_NumStates];
    memcpy(origX, m_X, sizeof(double)*m_NumStates);
    memcpy(origRates, m_Rates, sizeof(double)*m_NumStates);

    if (debug) {
        cerr << " origX: ";
        for (unsigned int i =0; i < m_NumStates;  ++i) {
            cerr << origX[i] << "\t";
        }
        cerr << endl;
    }

    // draw (stochastic) number of times each transition will occur
    vector<int> numTransitions(m_Nu.size(), 0);
    for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
        if (m_TransCats[j] == eNoncritical) {
            if (m_Rates[j]*tau > 1e8) {
                //for high rate, use normal to approx poisson.
                //should basically never yield negative, but just to
                //be sure, cap at 0
                numTransitions[j] = max(0.,floor(rnorm(m_Rates[j]*tau, sqrt(m_Rates[j]*tau))));
            } else {
                numTransitions[j] = rpois(m_Rates[j]*tau);
                //cerr << "nt[" << j << "] " << numTransitions[j] << endl;
            }
        }
    }

    // Calculate equation (7) terms not involving x[t+tau] and call this alpha:
    //   alpha = x + nu.(P - tau/2 R(x))
    // Also initialize iterative search for x[t+tau] at expectation (reset m_X)
    double* alpha = new double[m_NumStates];
    memcpy(alpha, m_X, sizeof(double)*m_NumStates);
    for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
        if (m_TransCats[j] == eNoncritical) {
            for (unsigned int k = 0;  k < m_Nu[j].size();  ++k) {
                alpha[m_Nu[j][k].m_State] += m_Nu[j][k].m_Mag * 
                    (numTransitions[j] - (tau/2)*m_Rates[j]);
                //reset m_X to expectation as our initial guess
                m_X[m_Nu[j][k].m_State] += m_Nu[j][k].m_Mag *
                    (tau/2)*m_Rates[j];
            }
        }
    }
    //expectations may send states negative; clamp!
    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
        if (m_X[i] < 0) {
            m_X[i] = 0;
        }
    }

    if (debug) {
        cerr << " alpha:";
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            cerr << " " << alpha[i];
        }
        cerr << endl;
        cerr << "it " << 0 << " newX: ";
        for (unsigned int i =0; i < m_NumStates;  ++i) {
            cerr << m_X[i] << "\t";
        }
        cerr << endl;
    }

    //a few variables needed by LAPACK
    int N = m_NumStates;
    int nrhs = 1;
    int *ipiv = new int[m_NumStates];
    int info;
    double *matrixA = new double[m_NumStates*m_NumStates];
    double *matrixB = new double[m_NumStates];

    
    //Use Newton's method to solve implicit equation:
    //  Let Y = x[t+tau]
    //  F(Y) = Y - alpha - nu.((tau/2)*R(Y))
    //Solve Jacobian(F(Y0)) Y1 = -F(Y0) for Y1 to iteratively approach solution
    //This eqn expands to (I - nu.((tau/2)Jacobian(R(Y0)))) Y1 = -F(Y0) where
    //the Jacobian of rates is supplied by the user.  The term to the
    //left of Y1 is called matrix A by LAPACK and the term on right is
    //matrix B.
    //
    //Perhaps should adjust max # of iterations..
    bool converged = false;
    unsigned int c = 0;
    while (++c <= 20  &&  !converged) {
        // Check to make sure we haven't taken too big a step --
        // i.e. no state variables should go negative
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            if (m_X[i] < 0) {
                delete[] origRates;
                delete[] alpha;
                delete[] ipiv;
                delete[] matrixA;
                delete[] matrixB;
                memcpy(m_X, origX, sizeof(double)*m_NumStates);
                delete[] origX;
                return false;
            }
        }

        // define matrix A
        double* rateJacobian = x_CalcJacobian();
        memset(matrixA, 0, m_NumStates*m_NumStates*sizeof(double));
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] == eNoncritical) {
                for (unsigned int k = 0;  k < m_Nu[j].size();  ++k) {
                    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                        //R stores matrices column-wise
                        //LAPACK stores matrices row-wise
                        matrixA[i*m_NumStates + m_Nu[j][k].m_State] +=
                            m_Nu[j][k].m_Mag * rateJacobian[j*m_NumStates + i];
                    }
                }
            }
        }
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            for (unsigned int i2 = 0;  i2 < m_NumStates;  ++i2) {
                matrixA[i*m_NumStates + i2] *= -tau/2;
            }
            matrixA[i*m_NumStates + i] += 1;
        }

        // define matrix B
        // m_X is now our proposed x[t+tau].  Note that m_X has changed
        // even in our first iteration (initialized to expected value).
        x_UpdateRates();
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            matrixB[i] = alpha[i] - m_X[i];
        }
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] == eNoncritical) {
                for (unsigned int k = 0;  k < m_Nu[j].size();  ++k) {
                    matrixB[m_Nu[j][k].m_State] += 
                        m_Nu[j][k].m_Mag * (tau/2) * m_Rates[j];
                }
            }
        }


    if (debug) {
        cerr << "A:" << endl;
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            for (unsigned int i2 = 0;  i2 < m_NumStates;  ++i2) {
                cerr << matrixA[i2*m_NumStates + i] << "\t";
            }
            cerr << endl;
        }

        cerr << "B:" << endl;
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            cerr << matrixB[i] << "\t";
        }
        cerr << endl;

        cerr << "a:" << endl;
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            cerr << m_Rates[j] << "\t";
        }
        cerr << endl;
    }

        //solve linear eqn
        F77_NAME(dgesv)(&N, &nrhs, matrixA, &N, ipiv, matrixB, &N, &info); 
        if (info != 0) {
            warning("warning: lapack ran into trouble solving implicit equation");
            break;
        }
        //matrixB now contains solution (change in X)
        double normDelta = 0, normX = 0;
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            m_X[i] += matrixB[i];
            normDelta += matrixB[i]*matrixB[i];
            normX += m_X[i] * m_X[i];
        }
        //cerr << "\tNorms: " << normDelta << "\t" << normX << endl;
        converged = (normDelta < normX * m_ITLConvergenceTol);
        if (debug) {
            /*
            cerr << "Delta: ";
            for (unsigned int i =0; i < m_NumStates;  ++i) {
                cerr << matrixB[i] << "\t";
            }
            cerr << endl;
            */
            cerr << "it " << c << " newX: ";
            for (unsigned int i =0; i < m_NumStates;  ++i) {
                cerr << m_X[i] << "\t";
            }
            cerr << endl;
            /*
            x_UpdateRates();
            double t[m_NumStates];
            memcpy(t, alpha, sizeof(double)*m_NumStates);
            for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
                if (m_TransCats[j] == eNoncritical) {
                    for (unsigned int k = 0;  k < m_Nu[j].size();  ++k) {
                        t[m_Nu[j][k].m_State] +=
                            m_Nu[j][k].m_Mag * (tau/2) * m_Rates[j];
                    }
                }
            }
            cerr << "     newX: ";
            for (unsigned int i =0; i < m_NumStates;  ++i) {
                cerr << t[i] << "\t";
            }
            cerr << endl;
            */
        }

    } // end of iterating for Newton's method
    if (!converged) {
        warning("ITL solution did not converge!");
    }

    //restore original rates to execute deterministic transitions
    memcpy(m_Rates, origRates, sizeof(double)*m_NumStates);
    x_AdvanceDeterministic(tau);

    delete[] origRates;
    delete[] alpha;
    delete[] ipiv;
    delete[] matrixA;
    delete[] matrixB;

    bool tauTooBig = false;
    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
        if (m_X[i] < 0) {
            tauTooBig = true;
            break;
        }
        if (!m_RealValuedVariables[i]) {
            m_X[i] = round(m_X[i]);
        }
    }
    if (tauTooBig) {
        memcpy(m_X, origX, sizeof(double)*m_NumStates);
        delete[] origX;
        return false;
    }
    delete[] origX;
    *m_T += tau;
    return true;
}

/*---------------------------------------------------------------------------*/
// PRE : tau value to use for step, list of "critical" transitions
// POST: whether single EXPLICIT tau step was successfully taken (m_X
// updated if so)
bool CStochasticEqns::x_SingleStepETL(double tau) {
    double *origX = new double[m_NumStates];
    memcpy(origX, m_X, sizeof(double)*m_NumStates);
    for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
        if (m_TransCats[j] == eNoncritical) {
            double k;
            if (m_Rates[j]*tau > 1e8) {
                //for high rate, use normal to approx poisson.
                //should basically never yield negative, but just to
                //be sure, cap at 0
                k = max(0.,floor(rnorm(m_Rates[j]*tau, sqrt(m_Rates[j]*tau))));
            } else {
                k = rpois(m_Rates[j]*tau);
            }
            for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                m_X[m_Nu[j][i].m_State] +=  k * m_Nu[j][i].m_Mag;
            }
        } else if (m_TransCats[j] == eDeterministic) {
            for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                m_X[m_Nu[j][i].m_State] += m_Nu[j][i].m_Mag * m_Rates[j] *
                    tau;
            }
        }
    }

    bool tauTooBig = false;
    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
        if (m_X[i] < 0) {
            tauTooBig = true;
            break;
        }
    }
    if (tauTooBig) {
        memcpy(m_X, origX, sizeof(double)*m_NumStates);
        delete[] origX;
        return false;
    }

    *m_T += tau;
    delete[] origX;
    return true;
}

/*---------------------------------------------------------------------------*/
// PRE : time at which to end simulation
// POST: single adaptive tau leaping step taken & time series updated.
// Implemented from Cao Y, Gillespie DT, Petzold LR. The Journal of Chemical
// Physics (2007).
void CStochasticEqns::x_SingleStepATL(double tf) {
    x_UpdateRates();
    EStepType stepType;

    //identify "critical" transitions
    double criticalRate = 0;
    double noncritRate = 0;
    {
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] == eDeterministic) {
                noncritRate += m_Rates[j];
                continue;
            }
            unsigned int minTimes = numeric_limits<unsigned int>::max();
            for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                if (m_Nu[j][i].m_Mag < 0  &&
                    m_X[m_Nu[j][i].m_State]/abs(m_Nu[j][i].m_Mag) < minTimes) {
                    minTimes = m_X[m_Nu[j][i].m_State]/abs(m_Nu[j][i].m_Mag);
                }
            }
            if (minTimes < m_Ncritical) {
                m_TransCats[j] = eCritical;
                criticalRate += m_Rates[j];
            } else {
                m_TransCats[j] = eNoncritical;
                noncritRate += m_Rates[j];
            }
        }
    }

    if (debug) {
        cerr << "critical rate: " << criticalRate << "\t" << "noncrit rate: " << noncritRate << endl;
    }
    if (criticalRate + noncritRate == 0) {
        *m_T = tf;//numeric_limits<double>::infinity();
        m_TimeSeries.push_back(STimePoint(*m_T, m_X, m_NumStates));
        return;
    }

    // calc explicit & implicit taus
    double tau1, tau2;
    double tauEx = x_TauEx();
    double tauIm = x_TauIm();
    if (debug) {
        cerr << "tauEx: " << tauEx << "  tauIm:" << tauIm << endl;
    }
    if (tauEx*m_Nstiff < tauIm) {
        stepType = eImplicit;
        tau1 = tauIm;
    } else {
        stepType = eExplicit;
        tau1 = tauEx;
    }
    if (tau1 > tf - *m_T) { //cap at the final simulation time
        tau1 = tf - *m_T;
    }
    if (tau1 > m_MaxTau) {
        tau1 = m_MaxTauFunc ? min(tau1, x_CalcUserMaxTau()) : m_MaxTau;
        if (debug) {
            cerr << "maxtau: " << tau1 << " (" <<
                (m_MaxTauFunc ? x_CalcUserMaxTau() : m_MaxTau) << ")" << endl;
        }
    }

    bool tauTooBig;
    do {
        tauTooBig = false;
        if (!(tau1 > 0)) { throwError("logic error at line " << __LINE__) }
        if (tau1 < m_ExactThreshold / (criticalRate + noncritRate)) {
            if (debug) {
                cerr << "Taking exact steps.. (tau1 = " << tau1 << ")" << endl;
            }
            stepType = eExact;
            for (unsigned int i = 0;
                 i < m_NumExactSteps[m_PrevStepType]  &&  *m_T < tf;  ++i) {
                if (i > 0) {
                    x_UpdateRates();
                }
                x_SingleStepExact(tf);
                if (*m_T == numeric_limits<double>::infinity()) {
                    //signal that rates = 0
                    *m_T = tf;
                }
                if (debug) {
                    cerr << *m_T << " -- ";
                    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                        cerr << m_X[i] << " ";
                    }
                    cerr << endl;
                }
                m_TimeSeries.push_back(STimePoint(*m_T, m_X, m_NumStates));
            }
        } else {
            tau2 = (criticalRate == 0) ? numeric_limits<double>::infinity() :
                rexp(1./criticalRate);
            if (stepType == eExplicit  ||
                (tau1 > tau2  &&  stepType == eImplicit  &&  tau2 <= tauEx)) {
                if (debug) {
                    cerr << "going explicit w/ tau = " << min(tau1, tau2) << endl;
                }
                tauTooBig = !x_SingleStepETL(min(tau1, tau2));
            } else {
                if (debug) {
                    cerr << "going implicit w/ tau = " << tau1 << endl;
                }
                tauTooBig = !x_SingleStepITL(tau1);
            }
            if (!tauTooBig) {
                if (tau1 > tau2) { //pick one critical transition
                    unsigned int j = x_PickCritical(criticalRate);
                    if (debug) {
                        cerr << "hittin' the critical (" << j << ")" << endl;
                    }
                    for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                        m_X[m_Nu[j][i].m_State] +=  m_Nu[j][i].m_Mag;
                        if (m_X[m_Nu[j][i].m_State] < 0) {
                            throwError("variable " << m_Nu[j][i].m_State+1 <<
                                       " went negative after executing "
                                       "transition " << j+1 << ".  Most likely "
                                       "either your rate calculation or "
                                       "transition matrix is flawed.");
                        }
                    }
                }

                m_TimeSeries.push_back(STimePoint(*m_T, m_X, m_NumStates));
            }
            if (debug) {
                cerr << *m_T << " -- ";
                for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                    cerr << m_X[i] << " ";
                }
                cerr << endl;
            }
        }
        if (tauTooBig) {
            if (debug) {
                cerr << "whoa.. knock that tau back down" << endl;
            }
            tau1 /= 2;
        }
    } while (tauTooBig);

    m_PrevStepType = stepType;
}

/*---------------------------------------------------------------------------*/
// Exported C entrypoints for calling from R
extern "C" {
    // Doing the event that occurs

	    // for the SEIR model
		int* seirtransitions(int nbVilles){
	        //Initializing
	        //4 for number of states, 8 for number of events
	        int *transitions = new int[nbVilles*4*nbVilles*8];
	        for(int i=0; i<nbVilles*4*nbVilles*8;i++) transitions[i] =0;

	        for (unsigned int i = 0;  i < nbVilles;  ++i) {
	            //susceptible birth; S<- S+1
	            transitions[0*nbVilles*4*nbVilles + i*4*nbVilles+0*nbVilles+i] = 1;
	            //susceptible death : S <-S-1
	            transitions[1*nbVilles*4*nbVilles + i*4*nbVilles+0*nbVilles+i] = -1;
	            //exposed death : E <- E-1
	            transitions[2*nbVilles*4*nbVilles + i*4*nbVilles+1*nbVilles+i] = -1;
	            //infectious death : I<- I-1
	            transitions[3*nbVilles*4*nbVilles + i*4*nbVilles+2*nbVilles+i] = -1;
	            //recovered death : R<-R-1
	            transitions[4*nbVilles*4*nbVilles + i*4*nbVilles+3*nbVilles+i] = -1;
	            // infection : S <-S-1; E<- E+1
	            transitions[5*nbVilles*4*nbVilles + i*4*nbVilles+0*nbVilles+i] = -1;
	            transitions[5*nbVilles*4*nbVilles + i*4*nbVilles+1*nbVilles+i] = 1;
	            // becoming infectious: E<-E-1; I<-I+1
	            transitions[6*nbVilles*4*nbVilles + i*4*nbVilles+1*nbVilles+i] = -1;
	            transitions[6*nbVilles*4*nbVilles + i*4*nbVilles+2*nbVilles+i] = 1;
	            // recovery : I<- I-1; R <-R+1
	            transitions[7*nbVilles*4*nbVilles + i*4*nbVilles+2*nbVilles+i] = -1;
	            transitions[7*nbVilles*4*nbVilles + i*4*nbVilles+3*nbVilles+i] = 1;
	        }

	        // Printing the result on the screen to check the correction of the results
	       /* 
	        for(int i=0; i<nbVilles*4*nbVilles*8;i++)
	            cout<<"i= "<<i<<"  trans="<<transitions[i]<<"   ";
	        cout<<endl;
	    */
	      return(transitions);
	    }
//-------------------------------------------------------------//
// Transitions for the SEIR model
    int* sirtransitions(int nbVilles){
        //Initializing
        // 3 for number of states, 6 for number of events
        int *transitions = new int[nbVilles*3*nbVilles*6];
        for(int i=0; i<nbVilles*3*nbVilles*6;i++) transitions[i] =0;

        // event that occurs
        for (unsigned int i = 0;  i < nbVilles;  ++i) {
            //susceptible birth; S<- S+1
            transitions[0*nbVilles*3*nbVilles + i*3*nbVilles+0*nbVilles+i] = 1;
            //susceptible death : S <-S-1
            transitions[1*nbVilles*3*nbVilles + i*3*nbVilles+0*nbVilles+i] = -1;
            //infectious death : I<- I-1
            transitions[2*nbVilles*3*nbVilles + i*3*nbVilles+1*nbVilles+i] = -1;
            //recovered death : R<-R-1
            transitions[3*nbVilles*3*nbVilles + i*3*nbVilles+2*nbVilles+i] = -1;
            // infection : S <-S-1; I<- I+1
            transitions[4*nbVilles*3*nbVilles + i*3*nbVilles+0*nbVilles+i] = -1;
            transitions[4*nbVilles*3*nbVilles + i*3*nbVilles+1*nbVilles+i] = 1;
            // recovery : I<- I-1; R <-R+1
            transitions[5*nbVilles*3*nbVilles + i*3*nbVilles+1*nbVilles+i] = -1;
            transitions[5*nbVilles*3*nbVilles + i*3*nbVilles+2*nbVilles+i] = 1;
        }
        // Printing the result on the screen to check the correction of the results
       /* 
        for(int i=0; i<nbVilles*3*nbVilles*6;i++)
            cout<<"i= "<<i<<"  trans="<<transitions[i]<<"   ";
        cout<<endl;        
        */
      return(transitions);
    }
    //---------------------------------------------------------------//
    /*
      Doing stochastic simulation for the SEIR and SIR models
      by using the "adaptive tau-leaping" algorithm described by Cao Y, Gillespie DT, Petzold LR. The Journal of Chemical Physics (2007)
      and by integrating the R and the C++.

      In order to the values of parameters and of variables from R.
        Step 1: we get these values from R that are in the R type.
        Step 2: we convert the values of the R type to these of the C++ type
      */
    SEXP ssesAdaptiveTau(SEXP s_x0, // initial values of variables
                         SEXP s_f, SEXP s_fJacob,
                         // values of parameters
                         SEXP sigma, SEXP gamma,SEXP mu,
                         SEXP nbCONTACT0, SEXP nbCONTACT1,
                         SEXP phiPHASE, SEXP probVISITER, SEXP probINFECTER, SEXP nbVilles,
                         SEXP periDISE,SEXP tmax,// time of simulation
                         SEXP s_deterministic, SEXP s_changebound,
                         SEXP s_tlparams, SEXP s_fMaxtau) {
        try{
            if (!isVector(s_x0)  ||  !isReal(s_x0)) {
                error("invalid vector of initial values");
            }
            if (!isNull(s_fJacob)  &&  !isFunction(s_fJacob)) {
                error("invalid Jacobian function");
            }
            if (length(tmax) != 1) {
                error("invalid final time");
            }

            if (!isVector(s_deterministic)  ||  !isLogical(s_deterministic)) {
                error("invalid deterministic parameter -- must be logical vector");
            }

            if (!isVector(s_changebound)  ||  !isReal(s_changebound)  || length(s_changebound) != length(s_x0)) {
                error("invalid relratechange");
            }

            if (!isNull(s_tlparams)  &&  !isVector(s_tlparams)) {
                error("tl.params must be a list");
            }

            if (!isNull(s_fMaxtau)  &&  !isFunction(s_fMaxtau)) {
                error("invalid maxTau function");
            }

            //Getting the values of parameters and variables
            //number of subpopulations
            //Number of cities in the metapopulation
            if(!isReal(nbVilles) || length(nbVilles)!=1)
                error("invalid number of cities.");
            int cnbVilles = INTEGER_VALUE(nbVilles);

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

            //simulation time tmax
            if (!isReal(tmax)  ||  length(tmax) != 1) {
                error("invalid value of simulation time,it should be a real number");
            }
            double ctmax=NUMERIC_VALUE(tmax);

            // Disease period
            if (!isReal(periDISE)  ||  length(periDISE) != 1) {
                error("invalid value of period, it should be a real number");
            }
            double cperiDISE=NUMERIC_VALUE(periDISE);

            // for metapopulation
            double ** arr_probVISITER = calculerProbVISITER(cnbVilles,cprobVISITER);

            // xay dung ma tran INFECTER
            double **arr_probINFECTER; create2D(cnbVilles,cnbVilles,cprobINFECTER,arr_probINFECTER);

            double **arrNbCONTACT0 = calculerNbCONTACT(cnbVilles,cnbCONTACT0,1);
            double **arrNbCONTACT1 = calculerNbCONTACT(cnbVilles,cnbCONTACT1,1.0);
            // xay dung ma tran VENIR de la ville k
            // initialiser les valeurs des variables
            //separate initial values of variables
            unsigned int nbTotalVar= length(s_x0); //total number of states
            double *rs_x0 = REAL(s_x0);
            double* cs_x0; create1D(nbTotalVar,0.0,cs_x0);
            for(int i=0; i<nbTotalVar; i++) { cs_x0[i] = rs_x0[i];}
            int nvar= 6;
            unsigned long **valSEIR0;
            valSEIR0 = new unsigned long*[cnbVilles];
            for (int i = 0; i <  cnbVilles; i++)
            {
                valSEIR0[i] = new unsigned long[nvar];
                valSEIR0[i][0]= 0;
                valSEIR0[i][iS]= (unsigned long)cs_x0[0*cnbVilles+i];//iS
                valSEIR0[i][iE]= (unsigned long)cs_x0[1*cnbVilles+i];//iE
                valSEIR0[i][iI]= (unsigned long)cs_x0[2*cnbVilles+i];//iI
                valSEIR0[i][iR]= (unsigned long)cs_x0[3*cnbVilles+i];//iR
                valSEIR0[i][iN]= valSEIR0[i][iS]+ valSEIR0[i][iE]+ valSEIR0[i][iI]+valSEIR0[i][iR];
            }
            /*
            cout<<"valSEIR0"<<endl;
            for(int i=0; i<cnbVilles; i++){
                for (int j=0; j<nvar; j++)
                    cout<<"    " << valSEIR0[i][j];
                cout<<endl;
            }
            cout<<"valSEIRAZZZZZZZ"<<endl;
            */


            double **arr_probVENIRk;  create2D(cnbVilles,cnbVilles,0, arr_probVENIRk);
            calculerProbVENIRk(cnbVilles,valSEIR0,arr_probVISITER,arr_probVENIRk);

            // vector of transitions (event occurs)
            int* transitions;
            unsigned int nbrate = 0;
            // check here is SEIR model or SIR model
            //-sigma 0.125   //sigma
            if (!isReal(sigma) || length(sigma)!=1) {
                error("invalid parameter sigma");
            }
            double csigma = NUMERIC_VALUE(sigma);
            // select Model
            bool flagSEIR=TRUE;
            if(csigma==INFINITY) flagSEIR=FALSE;
            // Doing simulation for SEIR/SIR models
            if(flagSEIR){//SEIR model
                transitions = seirtransitions(cnbVilles);
                nbrate=cnbVilles*8; //number of events for all subpopulations
            }
            else{//SIR model
                transitions = sirtransitions(cnbVilles);
                nbrate=cnbVilles*6; //number of events for all subpopulations
            }
            // stochastic simulation
            CStochasticEqns eqns(s_x0,transitions,nbrate,
                                 s_f, s_fJacob,
                                 //parameters
                                 cnbVilles,csigma,cgamma,cmu,arr_probVISITER,arrNbCONTACT0,arrNbCONTACT1,
                                 arr_probINFECTER,arr_probVENIRk,cphiPHASE,cperiDISE,
                                 //
                                 REAL(s_changebound),
                                 s_fMaxtau, s_deterministic);
            delete []transitions;
            for (int i = 0; i <  cnbVilles; i++){
                delete []arr_probVISITER[i];
                delete []arr_probINFECTER[i];
                delete []arrNbCONTACT0[i];
                delete []arrNbCONTACT1[i];
                delete []valSEIR0[i];
                delete []arr_probVENIRk[i];
            }
            delete []arr_probVISITER; delete[]arr_probINFECTER;delete []arrNbCONTACT0; delete[]arrNbCONTACT1;delete []valSEIR0;
            delete []arr_probVENIRk;
            delete []cphiPHASE; delete []cs_x0;
            if (!isNull(s_tlparams)) {
                eqns.SetTLParams(s_tlparams);
            }
            eqns.EvaluateATLUntil(REAL(tmax)[0]);
            return eqns.GetTimeSeriesSEXP();
        } catch (exception &e) {
            error(e.what());
            return R_NilValue;
        }
    }
}
//-----------------------------------------------------------------------//
