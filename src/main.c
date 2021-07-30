//To be included as the first lines in main.c
#include "R.h"
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
#include<time.h>
#include <R_ext/Utils.h>


#include "auxiliaryfuncs.h"
//#include "auxiliaryfuncs.c"
#include "allBF.h"
//#include "allBF.c"
#include "priorprob.h"
//#include "priorprob.c"
void gConst (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*Constpriorprob(p,dimensionNull);
	NormConstantPrior+=Constpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*Constpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*Constpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= gBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= gBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void gSB (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*SBpriorprob(p,dimensionNull);
	NormConstantPrior+=SBpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*SBpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*SBpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= gBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= gBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void gUser (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*gsl_vector_get(priorvector,dimensionNull);
	NormConstantPrior+=gsl_vector_get(priorvector,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= gBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= gBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void RobustConst (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*Constpriorprob(p,dimensionNull);
	NormConstantPrior+=Constpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*Constpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*Constpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= RobustBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= RobustBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void RobustSB (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*SBpriorprob(p,dimensionNull);
	NormConstantPrior+=SBpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*SBpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*SBpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= RobustBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= RobustBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void RobustUser (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*gsl_vector_get(priorvector,dimensionNull);
	NormConstantPrior+=gsl_vector_get(priorvector,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= RobustBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= RobustBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void LiangConst (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*Constpriorprob(p,dimensionNull);
	NormConstantPrior+=Constpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*Constpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*Constpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= LiangBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= LiangBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void LiangSB (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*SBpriorprob(p,dimensionNull);
	NormConstantPrior+=SBpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*SBpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*SBpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= LiangBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= LiangBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void LiangUser (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*gsl_vector_get(priorvector,dimensionNull);
	NormConstantPrior+=gsl_vector_get(priorvector,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= LiangBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= LiangBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void ZSConst (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*Constpriorprob(p,dimensionNull);
	NormConstantPrior+=Constpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*Constpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*Constpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= ZSBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= ZSBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void ZSSB (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*SBpriorprob(p,dimensionNull);
	NormConstantPrior+=SBpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*SBpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*SBpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= ZSBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= ZSBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void ZSUser (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*gsl_vector_get(priorvector,dimensionNull);
	NormConstantPrior+=gsl_vector_get(priorvector,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= ZSBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= ZSBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void flsConst (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*Constpriorprob(p,dimensionNull);
	NormConstantPrior+=Constpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*Constpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*Constpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= flsBF21fun(p, n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= flsBF21fun(p, n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void flsSB (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*SBpriorprob(p,dimensionNull);
	NormConstantPrior+=SBpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*SBpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*SBpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= flsBF21fun(p, n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= flsBF21fun(p, n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void flsUser (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*gsl_vector_get(priorvector,dimensionNull);
	NormConstantPrior+=gsl_vector_get(priorvector,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= flsBF21fun(p, n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= flsBF21fun(p, n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void intrinsicConst (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*Constpriorprob(p,dimensionNull);
	NormConstantPrior+=Constpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*Constpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*Constpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= intrinsicBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= intrinsicBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void intrinsicSB (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*SBpriorprob(p,dimensionNull);
	NormConstantPrior+=SBpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*SBpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*SBpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= intrinsicBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= intrinsicBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void intrinsicUser (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*gsl_vector_get(priorvector,dimensionNull);
	NormConstantPrior+=gsl_vector_get(priorvector,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= intrinsicBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= intrinsicBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void Robust2Const (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*Constpriorprob(p,dimensionNull);
	NormConstantPrior+=Constpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*Constpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*Constpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= Robust2BF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= Robust2BF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void Robust2SB (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*SBpriorprob(p,dimensionNull);
	NormConstantPrior+=SBpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*SBpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*SBpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= Robust2BF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= Robust2BF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void Robust2User (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*gsl_vector_get(priorvector,dimensionNull);
	NormConstantPrior+=gsl_vector_get(priorvector,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= Robust2BF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= Robust2BF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void geointrinsicConst (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*Constpriorprob(p,dimensionNull);
	NormConstantPrior+=Constpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*Constpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*Constpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= geointrinsicBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= geointrinsicBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void geointrinsicSB (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*SBpriorprob(p,dimensionNull);
	NormConstantPrior+=SBpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*SBpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*SBpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= geointrinsicBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= geointrinsicBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void geointrinsicUser (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*gsl_vector_get(priorvector,dimensionNull);
	NormConstantPrior+=gsl_vector_get(priorvector,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= geointrinsicBF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= geointrinsicBF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void geointrinsic2Const (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*Constpriorprob(p,dimensionNull);
	NormConstantPrior+=Constpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*Constpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*Constpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= geointrinsic2BF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= geointrinsic2BF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Constmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void geointrinsic2SB (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*SBpriorprob(p,dimensionNull);
	NormConstantPrior+=SBpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*SBpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*SBpriorprob(p,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= geointrinsic2BF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= geointrinsic2BF21fun(n, k2e ,knull,Q);
				
		//Now the main calculations:
		unnormPostProb=SBmainalgebraics(BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}

void geointrinsic2User (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time, int *pknull)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int knull=*pknull;
		
	int SAVE=*pSAVE;
	int StartAt= *pinicio;
	int FinishAt= *pfinal; //where to finish
	
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	//double info=pow(2,p);

	//Data files: (R version)
	char strtmp[100] = ""; 
	
	char nfileR1[100] = "/Dependent.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR1);
	strcpy(nfileR1,strtmp);
	FILE * fResponse = fopen(nfileR1, "r");
	
	char nfileR2[100] = "/Design.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR2);
	strcpy(nfileR2,strtmp);
	FILE * fDesign = fopen(nfileR2, "r");
	//-----------
	
	//The prior probabilities:
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
		
	
	//int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the null model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2e=0;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs:
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension
	//position 0 contains contains Pr(M|dim=knull+0,data)(=Pr(MNull|data)), position 1 contains Pr(M|dim=knull+1,data),
	//...position p contains Pr(M|dim=knull+p,data)(=Pr(Mfull|data)).
	gsl_vector * dimension_prob=gsl_vector_calloc(p+1);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;
	
	//the prior probabilities are read once in a vector
	//this vector is accesed to compute them
	//this vector only used in the "User"-type routines (are in all for
	//simplicity in maintenance)
	//priorvector which is read from a vector located in the working directory
	//call knull the dimension of the ORIGINAL null model
	//is a (p+1)-vector: position 0 contains Pr(M|dim=knull+0)(=Pr(MNull)), position 1 contains Pr(M|dim=knull+1),
	//...position p contains Pr(M|dim=knull+p)(=Pr(Mfull)).
	gsl_vector *priorvector = gsl_vector_calloc(p+1);
	gsl_vector_fscanf(fPriorProb, priorvector);
	fclose(fPriorProb);

	//Variables to store the results: We save the SAVE maximum probable models, stored in
	//Who_Max_SAVE with unnormalized post. probabilities in Max_SAVE_BF
	gsl_vector * Who_Max_SAVE = gsl_vector_calloc(SAVE);
	gsl_vector * Max_SAVE_BF=gsl_vector_calloc(SAVE);
	gsl_permutation * initperm = gsl_permutation_calloc(SAVE);
	gsl_permutation_init(initperm);
	//the vector with the posterior probs of the most probable models
	gsl_vector * postprob=gsl_vector_calloc(SAVE);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//Complete the results for the null model:
	int dimensionNull=0;
	NormConstant=1.0*gsl_vector_get(priorvector,dimensionNull);
	NormConstantPrior+=gsl_vector_get(priorvector,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*gsl_vector_get(priorvector,dimensionNull));
    
    
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model)*/
                

		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= geointrinsic2BF21fun(n, k2e, knull,Q);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				
		//save the model visited and the BF's
		gsl_vector_set(Who_Max_SAVE, (model-StartAt+1), model);
		gsl_vector_set(Max_SAVE_BF, (model-StartAt+1), unnormPostProb);
			
		gsl_vector_set_zero(index);
		}
	
		//Outside the LOOP I: order the values
		gsl_sort_vector_index(initperm, Max_SAVE_BF);
		gsl_permutation_reverse(initperm);
		gsl_permute_vector(initperm, Max_SAVE_BF);
		gsl_permute_vector(initperm, Who_Max_SAVE);


	// //////////////////////////////////////////////
	//LOOP II: Explore the rest of the models
		for (model=(SAVE+StartAt-1); model<(FinishAt+1); model++)
		{
		//printf("%Ld\n", model);
		//The progress of the calculations:
		R_CheckUserInterrupt();
		//if (model%1000000==0) Rprintf("Done(percentage): %- 3.5f\n", 
		//	(model-SAVE-StartAt)*100.0/(FinishAt-SAVE-StartAt));

		//Obtain Q=SSE/SSEnull and other calculations
		//also: i)index, which is returned, contains
		//the binary expression of the model;
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
				
		//the bayes factor in favor of modeli and against M0 
		k2e=k2+knull;
		BF21= geointrinsic2BF21fun(n, k2e ,knull,Q);
		//printf("The Bayes factor is: %.6f\n", BF21);
				
		//Now the main calculations:
		unnormPostProb=Usermainalgebraics(priorvector, BF21, p, index, 
                               &NormConstant, &NormConstantPrior, incl_prob, dimension_prob, 
                               k2, hatbetap, meanhatbetap, joint_incl_prob);
				


		//Check if the set of most probable models have to be recomputed
		if (unnormPostProb>gsl_vector_get(Max_SAVE_BF,(SAVE-1)))
				{
				gsl_vector_set(Max_SAVE_BF, (SAVE-1), unnormPostProb);
				gsl_vector_set(Who_Max_SAVE, (SAVE-1), model);
				recompute(Max_SAVE_BF, Who_Max_SAVE, SAVE);
				}

		//Reset
		gsl_vector_set_zero(index);

		}

	// //////////////////////////////////////////////
	//III Reescale the results and write the files:

	//reescale the inclusion probs: (divide incl_prob/NormConstant): 
	//PrintVector(incl_prob, p);
		gsl_vector_scale(incl_prob, 1.0/NormConstant);
	
		gsl_matrix_scale(joint_incl_prob, 1.0/NormConstant);
	
		//reescale the dimension probs:
		gsl_vector_scale(dimension_prob, 1.0/NormConstant);
		//reescale the posterior probs of the most probable models into the vector postprob:
		gsl_vector_memcpy(postprob, Max_SAVE_BF);
		gsl_vector_scale(postprob, 1.0/NormConstant);
		//reescale the betahat:
		gsl_vector_scale(meanhatbetap, 1.0/NormConstant);
	

		NormConstant=exp(log(NormConstant)-log(NormConstantPrior));

	
	//write the results to files (R version)
	char nfile1[100]="/PostProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fProb = fopen(strcat(nfile1,subindex), "w");
	
	char nfile2[100]="/MostProbModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile2);
	strcpy(nfile2,strtmp);
	FILE * fModels = fopen(strcat(nfile2,subindex), "w");
	
	char nfile3[100]="/InclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile3);
	strcpy(nfile3,strtmp);
	FILE * fInclusion = fopen(strcat(nfile3,subindex), "w");
	
	char nfile4[100]="/NormConstant";
	strcpy(strtmp,home);
	strcat(strtmp,nfile4);
	strcpy(nfile4,strtmp);	
	FILE * fTotalBF = fopen(strcat(nfile4,subindex), "w");
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
	
	char nfile6[100]="/StartEnd";
	strcpy(strtmp,home);
	strcat(strtmp,nfile6);
	strcpy(nfile6,strtmp);
	FILE * fStEnd = fopen(strcat(nfile6,subindex), "w");
	
	char nfile7[100]="/NormConstantPrior";
	strcpy(strtmp,home);
	strcat(strtmp,nfile7);
	strcpy(nfile7,strtmp);
	FILE * fNormCPrior = fopen(strcat(nfile7,subindex), "w");
	
	char nfile8[100]="/JointInclusionProb";
	strcpy(strtmp,home);
	strcat(strtmp,nfile8);
	strcpy(nfile8,strtmp);
	FILE * fJointInclusion = fopen(strcat(nfile8,subindex), "w");
	
	char nfile9[100]="/betahat";
	strcpy(strtmp,home);
	strcat(strtmp,nfile9);
	strcpy(nfile9,strtmp);
	FILE * fbetahat = fopen(strcat(nfile9,subindex), "w");
	//---------
	
		gsl_vector_fprintf(fbetahat, meanhatbetap, "%.20f");
		gsl_vector_fprintf(fInclusion, incl_prob, "%.20f");
		gsl_vector_fprintf(fProb, postprob, "%.20f");
		gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
		my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
		gsl_vector_fprintf(fModels, Who_Max_SAVE, "%f");
	

		fprintf(fTotalBF, "%.20f", NormConstant);
		fprintf(fNormCPrior, "%.20f", NormConstantPrior);
		fprintf(fStEnd, "%d %d", StartAt, FinishAt);
	

		fclose(fProb);
		fclose(fModels);	
		fclose(fInclusion);
		fclose(fTotalBF);
		fclose(fDim);
		fclose(fJointInclusion);
		fclose(fStEnd);
		fclose(fNormCPrior);
		fclose(fbetahat);
	
		gsl_vector_free (y);
		gsl_matrix_free (X);	
		gsl_vector_free (index);
	
		gsl_vector_free(incl_prob);
		gsl_matrix_free(joint_incl_prob);
		gsl_vector_free(dimension_prob);
		gsl_vector_free(Who_Max_SAVE);
		gsl_vector_free(Max_SAVE_BF);
		gsl_permutation_free(initperm);
		gsl_vector_free(postprob);
		gsl_vector_free(hatbetap);
		gsl_vector_free(meanhatbetap);
		
	
//	Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;
}
