//the prior probability of a model (depending on the dimension and the dimension of the full)
#include <R.h>
#include <gsl/gsl_sf_gamma.h>
#include<math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "Gibbsauxiliaryfuncs.h"



//Constant prior:
double Constpriorprob(int p, int dimensioni)
{
	//use this for Uniform:
	return(1.0);
}

//Scott and Berger prior
double SBpriorprob(int p, int dimensioni)
{
	//use this for Scott and Berger prior:
	return(exp(-gsl_sf_lnchoose(p, dimensioni)));
}

//Double Scott-Berger prior
double SBSBpriorprob(gsl_vector * indexfr, gsl_vector * positionsx, gsl_matrix * positions, int nofvars, gsl_vector * levels,
	                   int p, gsl_vector * isfactor)
{
 	//v is a vector that contains positions*index
	gsl_vector * v = gsl_vector_calloc(nofvars);
	double suma=0.0, m1=0.0, m1plusm2=0.0;
	int j=0;
	
  //m1 is the number of numerical vars in the model:		 
  for (int i = 0; i < p; i++){
  	m1=m1+gsl_vector_get(positionsx, i)*gsl_vector_get(indexfr, i);
	}
	
  for (int i = 0; i < nofvars; i++)
   {
		 suma=0.0; j=0;
		 while(suma<gsl_vector_get(levels, i) & j<p){
       suma=suma+gsl_matrix_get(positions,i,j)*gsl_vector_get(indexfr,j);	
			 			 	if (suma==gsl_vector_get(levels, i) & gsl_vector_get(isfactor,i)==1.0){
								gsl_vector_set(indexfr, j, 0);
			 	}
			 j++;			
			}
		 if (suma>0){
			 m1plusm2++;
		 	}
			
		 gsl_vector_set(v, i, suma);
   }
  
	 		
	//m2 is the number of factors in the model:
	double m2=m1plusm2-m1;
	double PrMg=log(nofvars+1.0)+gsl_sf_lnchoose(nofvars, m1+m2);
	for (int i=0; i<nofvars; i++){
		if (gsl_vector_get(v,i)>0){
		PrMg=PrMg+log(gsl_vector_get(levels, i))+gsl_sf_lnchoose(gsl_vector_get(levels, i), gsl_vector_get(v, i));
	}
	}
	
	PrMg=exp(-PrMg);
	
	//Rprintf("Prior probability is %.20f \n", PrMg);
	return(PrMg);
	
}

//Double Constant prior
double ConstConstpriorprob(gsl_vector * indexfr, gsl_vector * positionsx, gsl_matrix * positions, int nofvars, gsl_vector * levels,
	                   int p, gsl_vector * isfactor)
{
 	//v is a vector that contains positions*index
	gsl_vector * v = gsl_vector_calloc(nofvars);
	double suma=0.0, m1=0.0, m1plusm2=0.0;
	int j=0;
	
  //m1 is the number of numerical vars in the model:		 
  for (int i = 0; i < p; i++){
  	m1=m1+gsl_vector_get(positionsx, i)*gsl_vector_get(indexfr, i);
	}
	
  for (int i = 0; i < nofvars; i++)
   {
		 suma=0.0; j=0;
		 while(suma<gsl_vector_get(levels, i) & j<p){
       suma=suma+gsl_matrix_get(positions,i,j)*gsl_vector_get(indexfr,j);	
			 			 	if (suma==gsl_vector_get(levels, i) & gsl_vector_get(isfactor,i)==1.0){
								gsl_vector_set(indexfr, j, 0);
			 	}
			 j++;			
			}
		 if (suma>0){
			 m1plusm2++;
		 	}
			
		 gsl_vector_set(v, i, suma);
   }
  
	 		
	//m2 is the number of factors in the model:
	//double m2=m1plusm2-m1;
	double PrMg=nofvars*log(2);
	for (int i=0; i<nofvars; i++){
		if (gsl_vector_get(v,i)>0){
		PrMg=PrMg+log(pow(2.0, gsl_vector_get(levels, i))-1.0);
	}
	}
	
	PrMg=exp(-PrMg);
	
	//Rprintf("Prior probability is %.20f \n", PrMg);
	return(PrMg);
	
}

double SBConstpriorprob(gsl_vector * indexfr, gsl_vector * positionsx, gsl_matrix * positions, int nofvars, gsl_vector * levels,
	                   int p, gsl_vector * isfactor)
{
    //A prior that is SB at the level of variables and Constant within the factors (over the levels of the factors)
 	//v is a vector that contains positions*index
	gsl_vector * v = gsl_vector_calloc(nofvars);
	double suma=0.0, m1=0.0, m1plusm2=0.0;
	int j=0;
	
  //m1 is the number of numerical vars in the model:		 
  for (int i = 0; i < p; i++){
  	m1=m1+gsl_vector_get(positionsx, i)*gsl_vector_get(indexfr, i);
	}
	
  for (int i = 0; i < nofvars; i++)
   {
		 suma=0.0; j=0;
		 while(suma<gsl_vector_get(levels, i) & j<p){
       suma=suma+gsl_matrix_get(positions,i,j)*gsl_vector_get(indexfr,j);	
			 			 	if (suma==gsl_vector_get(levels, i) & gsl_vector_get(isfactor,i)==1.0){
								gsl_vector_set(indexfr, j, 0);
			 	}
			 j++;			
			}
		 if (suma>0){
			 m1plusm2++;
		 	}
			
		 gsl_vector_set(v, i, suma);
   }
  
	double m2=m1plusm2-m1;
	double PrMg=log(nofvars+1.0)+gsl_sf_lnchoose(nofvars, m1+m2);
	for (int i=0; i<nofvars; i++){
		if (gsl_vector_get(v,i)>0){
		PrMg=PrMg+log(pow(2.0, gsl_vector_get(levels, i))-1.0);
	}
	}
	
	PrMg=exp(-PrMg);

	return(PrMg);
	
}