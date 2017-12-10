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

#include "priorprob.h"

int gsl_vector_mifrac(gsl_vector *v, const double x){
	//v=v/x via logs: v=exp(log(v)/log(x))
	//OJO: only works for x>0
	int i=0;
	double sign=1.0, temp=1.0;
	double *p_v=v->data;
	for (i=0; i<v->size; i++)
	{
		temp=*p_v;
		sign=1.0;
		if (*p_v<0){
			sign=-1.0;
		}
        *p_v=sign*exp(log(fabs(temp))-log(fabs(x)));
        p_v++;
	}
	free(p_v);
	return 0;
}


//Next function in charge of the calculations previous to obtain the Bayes factor
//Returns the long double Q=SSE/SSEnull
double statistics(int model, int p, int n, double SSEnull, gsl_matrix * X, 
				  gsl_vector * y, gsl_vector * index, int *k2,
				  gsl_vector *hatbetap)
{
	
	*k2=0;
	int cont=0;
	int i=0;
	double SSE=0.0;
	
	
	int binmodel=model;
	//index is the binary expression of variable model
	while (binmodel != 0)
	{
		gsl_vector_set(index, cont, (binmodel%2));
		*k2+=binmodel%2;
		binmodel = binmodel/2;
		cont++;
	}	
	
	/*
	//This new piece of code by Enric Lazcorreta (August 2014):
	for (int j=0; j<p; j++){
		if (model && pow(2,j)) {
			gsl_vector_set(index, j, 1);
			*k2+=1;
		}
	}
	*/
	
	//hatbeta will store the mle of the k2-dimensional beta
	gsl_vector * hatbeta=gsl_vector_calloc(*k2);
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector_set_zero(hatbetap);
	
	//we construct the matriz design corresponding to the covariates in index
	//as Xindex=X*L
	gsl_matrix * L = gsl_matrix_calloc(p,*k2);
	gsl_matrix * Xindex = gsl_matrix_calloc(n,*k2);
	gsl_matrix_set_zero(L);
	//the vector who has the indexes of the covariates
	int otrocont=0;
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
	
	//Now the QR decomposition
	gsl_vector * tau = gsl_vector_calloc(*k2);
	gsl_linalg_QR_decomp(Xindex, tau); //Notice that now Xindex contains the QR decomposition
	
	
	gsl_vector * residual=gsl_vector_calloc(n);
	gsl_linalg_QR_lssolve (Xindex, tau, y, hatbeta, residual);
	
	
	gsl_vector_view wres = gsl_vector_subvector(residual, 0, n);	
	gsl_blas_ddot(&wres.vector, &wres.vector, &SSE);
	
	//To compute the hatbetap
	gsl_blas_dgemv(CblasNoTrans, 1.0, L, hatbeta, 0.0, hatbetap);
	
	
	
	double Q=exp(log(SSE)-log(SSEnull));
	
	gsl_matrix_free(L);
	gsl_matrix_free(Xindex);
	gsl_vector_free(tau);
	gsl_vector_free(hatbeta);
	gsl_vector_free(residual);
	
	return(Q);
}


//the function in charge of the main calculations. It depends on BF21,
//p and index (the binary expression of the model). Then updates NormConstant, NormConstantPrior, incl_prob, dimension_prob 
//. It returns unnormPostProb, used in the loop.

//Note: there are three mainalgebraics. Constmainalgebraics (uses Pr(Mi)=cte), SBmainalgebraics (uses Pr(Mi)=Scott and Berger) 
//and Usermainalgebraics (uses a vector defined by the user)

double Constmainalgebraics(double BF21, int p, gsl_vector * index, double *NormConstant, 
						double *NormConstantPrior, gsl_vector *incl_prob, 
						gsl_vector *dimension_prob,
						int k2, gsl_vector *hatbetap, gsl_vector *meanhatbetap, 
						gsl_matrix *joint_incl_prob)
{
	
	
	double acdimensioni=0.0;
	
	//Multiply the BF by the prior probability to obtain the unnormalized posterior probability
	//(this is the one that is going to be used in the calculations)
	double unnormPostProb=BF21*Constpriorprob(p,k2);
	
	//update the normalizing constants:
	//i)of the posterior
	*NormConstant+=unnormPostProb;
	
	//ii)and of the prior
	*NormConstantPrior+=Constpriorprob(p,k2);
	
	//update the inclusion probs:
	gsl_blas_daxpy(unnormPostProb, index, incl_prob);
	
	
	//update the joint inclusion probs:
	gsl_blas_dger(unnormPostProb, index, index, joint_incl_prob);
	
	//update the mean of betahats
	gsl_blas_daxpy(unnormPostProb, hatbetap, meanhatbetap);
	//ii)sum with the accumulated (unnormalized) probability of each dimension 
	acdimensioni=gsl_vector_get(dimension_prob,k2)+unnormPostProb;
	//iii)fill the corresponding new value with the updated one
	gsl_vector_set(dimension_prob,k2,acdimensioni);
	
    return(unnormPostProb);	
	
}


double SBmainalgebraics(double BF21, int p, gsl_vector * index, double *NormConstant, 
							 double *NormConstantPrior, gsl_vector *incl_prob, 
							 gsl_vector *dimension_prob,
							 int k2, gsl_vector *hatbetap, gsl_vector *meanhatbetap, 
							 gsl_matrix *joint_incl_prob)
{
	
	
	double acdimensioni=0.0;
	
	//Multiply the BF by the prior probability to obtain the unnormalized posterior probability
	//(this is the one that is going to be used in the calculations)
	double unnormPostProb=BF21*SBpriorprob(p,k2);
	
	//update the normalizing constants:
	//i)of the posterior
	*NormConstant+=unnormPostProb;
	
	//ii)and of the prior
	*NormConstantPrior+=SBpriorprob(p,k2);
	
	//update the inclusion probs:
	gsl_blas_daxpy(unnormPostProb, index, incl_prob);
	
	
	//update the joint inclusion probs:
	gsl_blas_dger(unnormPostProb, index, index, joint_incl_prob);
	
	//update the mean of betahats
	gsl_blas_daxpy(unnormPostProb, hatbetap, meanhatbetap);
	//ii)sum with the accumulated (unnormalized) probability of each dimension 
	acdimensioni=gsl_vector_get(dimension_prob,k2)+unnormPostProb;
	//iii)fill the corresponding new value with the updated one
	gsl_vector_set(dimension_prob,k2,acdimensioni);
	
    return(unnormPostProb);	
	
}

double Usermainalgebraics(gsl_vector *priorvector, double BF21, int p, gsl_vector * index, double *NormConstant, 
						double *NormConstantPrior, gsl_vector *incl_prob, 
						gsl_vector *dimension_prob,
						int k2, gsl_vector *hatbetap, gsl_vector *meanhatbetap, 
						gsl_matrix *joint_incl_prob)
{
	
	
	double acdimensioni=0.0;
	
	//Multiply the BF by the prior probability to obtain the unnormalized posterior probability
	//(this is the one that is going to be used in the calculations)
	double unnormPostProb=BF21*gsl_vector_get(priorvector,k2);
	
	//update the normalizing constants:
	//i)of the posterior
	*NormConstant+=unnormPostProb;
	
	//ii)and of the prior
	*NormConstantPrior+=gsl_vector_get(priorvector,k2);
	
	//update the inclusion probs:
	gsl_blas_daxpy(unnormPostProb, index, incl_prob);
	
	
	//update the joint inclusion probs:
	gsl_blas_dger(unnormPostProb, index, index, joint_incl_prob);
	
	//update the mean of betahats
	gsl_blas_daxpy(unnormPostProb, hatbetap, meanhatbetap);
	//ii)sum with the accumulated (unnormalized) probability of each dimension 
	acdimensioni=gsl_vector_get(dimension_prob,k2)+unnormPostProb;
	//iii)fill the corresponding new value with the updated one
	gsl_vector_set(dimension_prob,k2,acdimensioni);
	
    return(unnormPostProb);	
	
}



//A function given a vector v of size N (and a vector of indexes w), ordered in descending order except the last one
//which resturn the vectors reordered incorporating the last one 
void recompute(gsl_vector * v, gsl_vector * w, int N)
{
	int i=N-1;
	double copiaBF=0.0;
	long double copiaInd=0;
	while ((gsl_vector_get(v, i)>gsl_vector_get(v,i-1)) && (i>=2))
	{
		copiaBF=gsl_vector_get(v,i-1);
		gsl_vector_set(v,i-1,gsl_vector_get(v,i));
		gsl_vector_set(v,i,copiaBF);
		copiaInd=gsl_vector_get(w,i-1);
		gsl_vector_set(w,i-1,gsl_vector_get(w,i));
		gsl_vector_set(w,i,copiaInd);
		i=i-1;
	}
	if (i==1 && gsl_vector_get(v, i)>gsl_vector_get(v,i-1))
	{
		copiaBF=gsl_vector_get(v,i-1);
		gsl_vector_set(v,i-1,gsl_vector_get(v,i));
		gsl_vector_set(v,i,copiaBF);
		copiaInd=gsl_vector_get(w,i-1);
		gsl_vector_set(w,i-1,gsl_vector_get(w,i));
		gsl_vector_set(w,i,copiaInd);
	}
	
}


int my_gsl_matrix_fprintf(FILE *stream,gsl_matrix *m,char *fmt)
{
	size_t rows=m->size1;
	size_t cols=m->size2;
	size_t row,col,ml;
	int fill;
	char buf[100];
	gsl_vector *maxlen;
	
	maxlen=gsl_vector_alloc(cols);
	for (col=0;col<cols;++col) {
		ml=0;
		for (row=0;row<rows;++row) {
			sprintf(buf,fmt,gsl_matrix_get(m,row,col));
			if (strlen(buf)>ml)
				ml=strlen(buf);
		}
		gsl_vector_set(maxlen,col,ml);
	}
	
	for (row=0;row<rows;++row) {
		for (col=0;col<cols;++col) {
			sprintf(buf,fmt,gsl_matrix_get(m,row,col));
			fprintf(stream,"%s",buf);
			fill=gsl_vector_get(maxlen,col)+1-strlen(buf);
			while (--fill>=0)
				fprintf(stream," ");
		}
		fprintf(stream,"\n");
	}
	gsl_vector_free(maxlen);
	return 0;
}

//A function to print a vector in a file by rows
int my_gsl_vector_fprintf(FILE *stream, gsl_vector *v, char *fmt)
{
	size_t length=v->size;
	size_t i=0;
	
	for (i=0;i<length;++i) {
		fprintf(stream, fmt, gsl_vector_get(v, i));
		fprintf(stream, " ");
	}
	
	fprintf(stream,"\n");
	
	return 0;
}

