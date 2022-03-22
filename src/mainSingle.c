#include "R.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <R_ext/Utils.h>
#include<time.h>
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h>

#include "allBF.h"
//#include "allBF.c"

void gBF (int *pn, int *pk2, int *pk0, double *pQ, double *B21)
{
	void R_CheckUserInterrupt(void);
	gsl_set_error_handler_off();
	
	//PARAMETERS: (R version)
	int n=*pn;
	int k2=*pk2;
	int k0=*pk0;
	double Q=*pQ;
	*B21=gBF21fun(n, k2, k0, Q);

}

void flsBF (int *pp, int *pn, int *pk2, int *pk0, double *pQ, double *B21)
{
    void R_CheckUserInterrupt(void);
    gsl_set_error_handler_off();
    
    //PARAMETERS: (R version)
    int n=*pn;
    int k2=*pk2;
    int k0=*pk0;
    int p=*pp;
    double Q=*pQ;
    *B21=flsBF21fun(p, n, k2, k0, Q);
    
}


void RobustBF (int *pn, int *pk2, int *pk0, double *pQ, double *B21)
{
	void R_CheckUserInterrupt(void);
	gsl_set_error_handler_off();
	
	//PARAMETERS: (R version)
	int n=*pn;
	int k2=*pk2;
	int k0=*pk0;	
	double Q=*pQ;
	*B21=RobustBF21fun(n, k2, k0, Q);

}


void LiangBF (int *pn, int *pk2, int *pk0, double *pQ, double *B21)
{
	void R_CheckUserInterrupt(void);
	gsl_set_error_handler_off();
	
	//PARAMETERS: (R version)
	int n=*pn;
	int k2=*pk2;
	int k0=*pk0;		
	double Q=*pQ;
	*B21=LiangBF21fun(n, k2, k0, Q);

}


void ZSBF (int *pn, int *pk2, int *pk0, double *pQ, double *B21)
{
	void R_CheckUserInterrupt(void);
	gsl_set_error_handler_off();
	
	//PARAMETERS: (R version)
	int n=*pn;
	int k2=*pk2;
	int k0=*pk0;			
	double Q=*pQ;
	*B21=ZSBF21fun(n, k2, k0, Q);

}

void intrinsicBF (int *pn, int *pk2, int *pk0, double *pQ, double *B21)
{
	void R_CheckUserInterrupt(void);
	gsl_set_error_handler_off();
	
	//PARAMETERS: (R version)
	int n=*pn;
	int k2=*pk2;
	int k0=*pk0;			
	double Q=*pQ;
	*B21=intrinsicBF21fun(n, k2, k0, Q);

}

void geointrinsicBF (int *pn, int *pk2, int *pk0, double *pQ, double *B21)
{
	void R_CheckUserInterrupt(void);
	gsl_set_error_handler_off();
	
	//PARAMETERS: (R version)
	int n=*pn;
	int k2=*pk2;
	int k0=*pk0;			
	double Q=*pQ;
	*B21=geointrinsicBF21fun(n, k2, k0, Q);

}

void geointrinsic2BF (int *pn, int *pk2, int *pk0, double *pQ, double *B21)
{
	void R_CheckUserInterrupt(void);
	gsl_set_error_handler_off();
	
	//PARAMETERS: (R version)
	int n=*pn;
	int k2=*pk2;
	int k0=*pk0;			
	double Q=*pQ;
	*B21=geointrinsic2BF21fun(n, k2, k0, Q);

}
