#include <R.h>
//#include <stdio.h>
//#include<math.h>
//#include<stdlib.h>
//#include<string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_heapsort.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

//Next function in charge of the calculations previous to obtain the Bayes factor
//Returns the long double Q=SSE/SSEnull
double Gibbsstatistics(int p, int n, double SSEnull, gsl_matrix * X, 
				  gsl_vector * y, gsl_vector * index, int *k2,
				  gsl_vector *hatbetap)
{
	
	
	//hatbeta will store the mle of the k2-dimensional beta
	//gsl_vector * hatbeta=gsl_vector_calloc(*k2);
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	//gsl_vector_set_zero(hatbetap);
	
	//we construct the matriz design corresponding to the covariates in index
	//as Xindex=X*L
	//gsl_matrix * L = gsl_matrix_calloc(p,*k2);
	gsl_matrix * Xindex = gsl_matrix_calloc(n,*k2);
	//gsl_matrix_set_zero(L);
	//the vector who has the indexes of the covariates
	int cont=0; 
    int otrocont=0;
	int i=0;
	double SSE=0.0;

//  for (cont=0; cont<(*k2); cont++){
    while (cont < (*k2)) {
        if (gsl_vector_get(index,otrocont)==1)
        {
            for (i=0; i<n; i++){
                gsl_matrix_set(Xindex, i, cont, gsl_matrix_get(X, i, otrocont));
            }
            cont++;
        }
        otrocont++;
    }
    /*
     for (i=0; i<p; i++)
	{
		if (gsl_vector_get(index,i)==1)
		{
			gsl_matrix_set(L,i,otrocont,1.0);
			//gsl_vector_int_set(who, otrocont, i);
			otrocont++;
		}
	}
	
	//Now multiply X*L->Xindex
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X, L, 0.0, Xindex);
    */
	
	//Now the QR decomposition
	gsl_vector * tau = gsl_vector_calloc(*k2);
	gsl_linalg_QR_decomp(Xindex, tau); //Notice that now Xindex contains the QR decomposition

    gsl_vector * residual=gsl_vector_alloc(n);
    gsl_vector_memcpy (residual, y);
    
    gsl_linalg_QR_QTvec(Xindex, tau, residual);
    
	
    
	//gsl_linalg_QR_lssolve (Xindex, tau, y, hatbeta, residual);
	//gsl_vector_view wres = gsl_vector_subvector(residual, 0, n);
    
	gsl_vector_view wres = gsl_vector_subvector(residual, *k2, n-(*k2));
    gsl_blas_ddot(&wres.vector, &wres.vector, &SSE);
	
	//To compute the hatbetap
	//gsl_blas_dgemv(CblasNoTrans, 1.0, L, hatbeta, 0.0, hatbetap);
	
	
	
	double Q=exp(log(SSE)-log(SSEnull));
	
	//gsl_matrix_free(L);
	gsl_matrix_free(Xindex);
	gsl_vector_free(tau);
	//gsl_vector_free(hatbeta);
	gsl_vector_free(residual);
	
	return(Q);
}


