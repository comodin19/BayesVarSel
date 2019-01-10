#include "R.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <R_ext/Utils.h>
//#include <gsl/gsl_sf_hyperg.h>
//#include <gsl/gsl_permutation.h>
//#include <gsl/gsl_permute_vector.h>
//#include <gsl/gsl_heapsort.h>
//#include <gsl/gsl_sort.h>
//#include <gsl/gsl_sort_vector.h>
#include<time.h>
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h>


#include "auxiliaryfuncs.h"
//#include "auxiliaryfuncs.c"
#include "Gibbsauxiliaryfuncs.h"
//#include "Gibbsauxiliaryfuncs.c"
#include "allBF.h"
//#include "allBF.c"
#include "priorprob.h"

void GibbsFSBSB (char *pI[], int *pn, int *pp, int *pSAVE, char *homePath[], int *pBurnin, double *time, int *pknull, int *pnthin, int *pseed)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	//18-12-18 version of Gibbs with Factors and prior probabilities obtained with SB-SB
	//the type of Bayes factor is defined in a file called typeofBF.txt 
	
	void R_CheckUserInterrupt(void);
	
	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();
	
	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int SAVE=*pSAVE;
	int knull=*pknull;
	int nthin=*pnthin;
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
    gsl_rng * ran = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(ran, *pseed);
	//double info=pow(2,p);
	
	Rprintf("The full design matrix has %d columns and n=%d rows\n", p, n);
	//Rprintf("The whole problem would have %f different competing models\n", info);
	
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
	
	
	//Factors:
	//positions is a file that contains, in binary, the positions of each of the vars
	//(either factor and its levels or numerical vars) in the design matrix
	char nfileR12[100] = "/positions.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR12);
	strcpy(nfileR12,strtmp);
	FILE * fpositions = fopen(nfileR12, "r");
	int i;
	int nofvars=0; //Number of different explanatory variables (either factors or numerical)
  while((i=fgetc(fpositions))!=EOF)
    {
        if (i == '\n')
            nofvars++;
    }
		
  Rprintf("Number of explanatory variables: %d\n", nofvars);
	
	//Factors: I close and open to start the reading from the beginning
	fclose(fpositions);		
	fpositions = fopen(nfileR12, "r");
		
		
	//Factors:
	//positionsx is a file that contains, in binary, the positions of numerical vars) in the design matrix
	char nfileR13[100] = "/positionsx.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR13);
	strcpy(nfileR13,strtmp);
	FILE * fpositionsx = fopen(nfileR13, "r");
			
	
	//The initial model:
	char nfileR3[100] = "/initialmodel.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR3);
	strcpy(nfileR3,strtmp);
	FILE * finitM = fopen(nfileR3, "r");
		
	//The prior probabilities:
	//Factor:
	/*
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
	*/
		
	//File that contain all visited models after burnin
	char nfile10[100]="/AllModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile10);
	strcpy(nfile10,strtmp);
	FILE * fAllModels = fopen(strcat(nfile10,subindex), "w");
	
	//File that contain all BF's of previuos models
	char nfile11[100]="/AllBF";
	strcpy(strtmp,home);
	strcat(strtmp,nfile11);
	strcpy(nfile11,strtmp);
	FILE * fAllBF = fopen(strcat(nfile11,subindex), "w");
	
	
	int cont=0; //cont is used for the thining
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);
	
	
	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	
	
  //Factors:
	gsl_vector * positionsx = gsl_vector_alloc(p);
  gsl_matrix * positions = gsl_matrix_calloc(nofvars, p);	
	gsl_vector_fscanf(fpositionsx, positionsx);
	fclose(fpositionsx);
	gsl_matrix_fscanf(fpositions, positions);
	fclose(fpositions);	
	
	//Rprintf("Matrix positions:\n");
	//PrintMatrix(positions, nofvars, p);
	
	//Rprintf("Vector positions of x:\n");
	//PrintVector(positionsx, p);
	
	//Factors:
	//A vector with the number of "levels" for each variable:
	int lev=0;
	gsl_vector * levels=gsl_vector_calloc(nofvars);
  for (i=0; i<nofvars; i++){
		lev=0;
		for (int j=0; j<p; j++){
			lev=lev+gsl_matrix_get(positions, i, j);
		}
		gsl_vector_set(levels, i, lev);		
  }

	//Rprintf("Vector number of levels:\n");
	//PrintVector(levels, nofvars);
	
	//Factors:
	//A vector of the same length as the number of variables (factors or nums)
	//with a one in the position of a vector (its corresponding row in positionsx
	//has more than one element being 1)
	gsl_vector * isfactor = gsl_vector_calloc(nofvars);
	double suma=0.0;
	  for (int i = 0; i < nofvars; i++)
	   {
			 suma=0.0;	 
			 for (int j=0; j<p; j++){
	       suma=suma+gsl_matrix_get(positions,i,j);				
			 }
			 if (suma>1){
				 gsl_vector_set(isfactor, i, 1.0);
			 }
	   }
	 
	//PrintVector(isfactor, nofvars);
	
	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);

	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2fr=0; //The number of covs in the full rank representation
	
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
	//indexfr is a full rank copy of index
	gsl_vector * indexfr = gsl_vector_calloc(p);
	
	//We fill it with the initial model
	gsl_vector_fscanf(finitM, index);
	fclose(finitM);

	//Factors:
	/*This is not used:
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
	*/
	
	//The HPM of the models visited:
	gsl_vector * HPM = gsl_vector_calloc(p);
	gsl_vector_memcpy(HPM, index);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	
	//this file contains the type of Bayes factor used:
	//0 is unitary; 1 is Robust; 2 is gprior; 4 is Liang and 5 is Zellner-Siow
	//fls is not yet implemented because it depends on an extra parameter
	char nfileP1[100] = "/typeofBF.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileP1);
	strcpy(nfileP1,strtmp);
	FILE * fpriorconfig = fopen(nfileP1, "r");
	//it is read and then written in an integer variable called typeofBF
	int typeofBF=0, ttt=0; //in some compilers we need to assign the value of fscanf (which is an integer)
	//correspoding to the number of input items successfully matched and assigned
	ttt = fscanf(fpriorconfig, "%d", &typeofBF);
	if (ttt != 1) {
		Rprintf("Error reading file containing the type of Bayes factor\n");
	}
	fclose(fpriorconfig);	
	//Now define the function for Bayes factors which is a pointer
	double (*BF21fun)(int, int, int, double)=NULL;
	if (typeofBF==0) BF21fun=unitBF21fun;
	if (typeofBF==1) BF21fun=RobustBF21fun;
	if (typeofBF==2) BF21fun=gBF21fun;
	if (typeofBF==4) BF21fun=LiangBF21fun;
	if (typeofBF==5) BF21fun=ZSBF21fun;
	
	
	
	// //////////////////////////////////////////////
	int iter=1, component=1, oldcomponent=1, newcomponent=1;
	
	//In what follows we call, for a given model Mi, PBF=BF(Mi vs M0)*Pr(Mi)
	double oldPBF=0.0, newPBF=0.0;
	//PrMg will contain the probability of each sampled model
	double PrMg=0.0;
	
    k2=(int) gsl_blas_dasum(index);
    if (k2>0){
			//Copy the sampled model (index) to a new one (indexfr) that will contain its full rank dimension copy
			//to compute Q and the Bayes factor. The real model, index, is used for the rest of purposes
			//but be CAREFUL as each time the SBSBpriorprob is called the "model" is changed
		  gsl_vector_memcpy(indexfr, index);
			PrMg=SBSBpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
	    k2fr=(int) gsl_blas_dasum(indexfr);			
      Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
	    oldPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
    }
    else{
		  gsl_vector_memcpy(indexfr, index);	
			PrMg=SBSBpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
			gsl_vector_set_zero(hatbetap);
      oldPBF= 1.0*PrMg;
 	  }
    
	double HPMBF=oldPBF;
	double ratio=0.0;

	//Burnin
	//Interpret old as current and new as proposal
	for (iter=1; iter<(*pBurnin+1); iter++){
		for (component=0; component<p; component++)
		{
			oldcomponent=gsl_vector_get(index, component);
			gsl_vector_set(index, component, 1-oldcomponent);
            k2=(int) gsl_blas_dasum(index);
            if (k2>0){
        			//Copy the sampled model (index) to a new one (indexfr) that will contain its full rank dimension copy
			        //to compute the Bayes factor. The real model, index, is used for the rest of purposes
						  gsl_vector_memcpy(indexfr, index);
							PrMg=SBSBpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
							k2fr=(int) gsl_blas_dasum(indexfr);			
              Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
		          newPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
            }
            else{
						  gsl_vector_memcpy(indexfr, index);	
							PrMg=SBSBpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
							gsl_vector_set_zero(hatbetap);
				      newPBF= 1.0*PrMg;
            }
						
            ratio=(oldcomponent*(oldPBF-newPBF)+newPBF)/(newPBF+oldPBF);
            newcomponent=gsl_ran_bernoulli(ran, ratio);
		      	if (newcomponent==oldcomponent){
	       	 		gsl_vector_set(index, component, newcomponent);
			      	k2=k2-1+2*oldcomponent;
			      }
		      	else
			      	oldPBF=newPBF;
	        	}
		R_CheckUserInterrupt();
	}
	
	//main loop
	//Interpret old as current and new as proposal
	for (iter=1; iter<(SAVE+1); iter++){
		cont=0;
		while (cont<nthin){
		for (component=0; component<p; component++)
		{
			oldcomponent=gsl_vector_get(index, component);
			gsl_vector_set(index, component, 1-oldcomponent);
            k2=(int) gsl_blas_dasum(index);
						if (k2>0){
						  gsl_vector_memcpy(indexfr, index);
							PrMg=SBSBpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
							k2fr=(int) gsl_blas_dasum(indexfr);			
              Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
		          newPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
            }
            else{
						  gsl_vector_memcpy(indexfr, index);	
							PrMg=SBSBpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
							gsl_vector_set_zero(hatbetap);
				      newPBF= 1.0*PrMg;
            }

              ratio=(oldcomponent*(oldPBF-newPBF)+newPBF)/(newPBF+oldPBF);
			        newcomponent=gsl_ran_bernoulli(ran, ratio);
		        	if (newcomponent==oldcomponent){
			          	gsl_vector_set(index, component, newcomponent);
				          k2=k2-1+2*oldcomponent;
		        	}
	        		else
	          			oldPBF=newPBF;
          		}
			cont++;
			R_CheckUserInterrupt();
		}
		
		/*
		Rprintf("\n\n\n");
		Rprintf("Accepted model:\n");
		PrintVector(index, p);
		Rprintf("Su dimension:%d", k2);
		Rprintf("Bg0*PrMg=%.10f\n", oldPBF);
		*/
		
		//update the summaries
		//inclusion probs:
		gsl_blas_daxpy(1.0, index, incl_prob);		
		//joint inclusion probs:
		gsl_blas_dger(1.0, index, index, joint_incl_prob);		
		//mean of betahats
		gsl_blas_daxpy(1.0, hatbetap, meanhatbetap);		
		//probability of dimensions:
		gsl_vector_set(dimension_prob,k2,gsl_vector_get(dimension_prob,k2)+1.0);
		//HPM
		if (oldPBF>HPMBF)
			{
			HPMBF=oldPBF;
			gsl_vector_memcpy(HPM, index);
			}
		//Write to the file the visited model
		my_gsl_vector_fprintf(fAllModels, index, "%f");
		//and the BF's
		fprintf(fAllBF, "%.20f\n", oldPBF/PrMg);  
      

	}
	
	//normalize:
	gsl_vector_scale(incl_prob, 1.0/SAVE);
	gsl_vector_scale(meanhatbetap, 1.0/SAVE);
	gsl_vector_scale(dimension_prob, 1.0/SAVE);
	gsl_matrix_scale(joint_incl_prob, 1.0/SAVE);
	
		
	//Write the results:
	char nfile1[100]="/LastModel";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fLastModel = fopen(strcat(nfile1,subindex), "w");
	
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
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
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
	gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
	my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
	gsl_vector_fprintf(fModels, HPM, "%f");
	gsl_vector_fprintf(fLastModel, index, "%f");
	
	
	fclose(fLastModel);
	fclose(fModels);	
	fclose(fInclusion);
	fclose(fDim);
	fclose(fJointInclusion);
	fclose(fbetahat);
	fclose(fAllModels);
	fclose(fAllBF);	

	
	gsl_vector_free (y);
	gsl_matrix_free (X);	
	gsl_vector_free (index);
	
	gsl_vector_free(positionsx);
	gsl_vector_free(levels);
	gsl_vector_free(isfactor);
	gsl_vector_free(indexfr);
	gsl_matrix_free(positions);
	
	
	gsl_vector_free(incl_prob);
	gsl_matrix_free(joint_incl_prob);
	gsl_vector_free(dimension_prob);
	gsl_vector_free(hatbetap);
	gsl_vector_free(meanhatbetap);
	
	free(BF21fun);
		
	gsl_rng_free (ran);
		
	//Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);	
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;

}

void GibbsFSBConst (char *pI[], int *pn, int *pp, int *pSAVE, char *homePath[], int *pBurnin, double *time, int *pknull, int *pnthin, int *pseed)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	//9-5-18 version of GibbsRobust with Factors and prior probabilities obtained with SB-SB
  //18-12-18 
  //the type of Bayes factor is defined in a file called typeofBF.txt 
	
	void R_CheckUserInterrupt(void);
	
	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();
	
	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int SAVE=*pSAVE;
	int knull=*pknull;
	int nthin=*pnthin;
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
    gsl_rng * ran = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(ran, *pseed);
	//double info=pow(2,p);
	
	Rprintf("The full design matrix has %d columns and n=%d rows\n", p, n);
	//Rprintf("The whole problem would have %f different competing models\n", info);
	
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
	
	
	//Factors:
	//positions is a file that contains, in binary, the positions of each of the vars
	//(either factor and its levels or numerical vars) in the design matrix
	char nfileR12[100] = "/positions.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR12);
	strcpy(nfileR12,strtmp);
	FILE * fpositions = fopen(nfileR12, "r");
	int i;
	int nofvars=0; //Number of different explanatory variables (either factors or numerical)
  while((i=fgetc(fpositions))!=EOF)
    {
        if (i == '\n')
            nofvars++;
    }
		
  Rprintf("Number of explanatory variables: %d\n", nofvars);
	
	//Factors: I close and open to start the reading from the beginning
	fclose(fpositions);		
	fpositions = fopen(nfileR12, "r");
		
		
	//Factors:
	//positionsx is a file that contains, in binary, the positions of numerical vars) in the design matrix
	char nfileR13[100] = "/positionsx.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR13);
	strcpy(nfileR13,strtmp);
	FILE * fpositionsx = fopen(nfileR13, "r");
			
	
	//The initial model:
	char nfileR3[100] = "/initialmodel.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR3);
	strcpy(nfileR3,strtmp);
	FILE * finitM = fopen(nfileR3, "r");
		
	//The prior probabilities:
	//Factor:
	/*
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
	*/
		
	//File that contain all visited models after burnin
	char nfile10[100]="/AllModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile10);
	strcpy(nfile10,strtmp);
	FILE * fAllModels = fopen(strcat(nfile10,subindex), "w");
	
	//File that contain all BF's of previuos models
	char nfile11[100]="/AllBF";
	strcpy(strtmp,home);
	strcat(strtmp,nfile11);
	strcpy(nfile11,strtmp);
	FILE * fAllBF = fopen(strcat(nfile11,subindex), "w");
	
	
	int cont=0; //cont is used for the thining
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);
	
	
	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	
	
  //Factors:
	gsl_vector * positionsx = gsl_vector_alloc(p);
  gsl_matrix * positions = gsl_matrix_calloc(nofvars, p);	
	gsl_vector_fscanf(fpositionsx, positionsx);
	fclose(fpositionsx);
	gsl_matrix_fscanf(fpositions, positions);
	fclose(fpositions);	
	
	//Rprintf("Matrix positions:\n");
	//PrintMatrix(positions, nofvars, p);
	
	//Rprintf("Vector positions of x:\n");
	//PrintVector(positionsx, p);
	
	//Factors:
	//A vector with the number of "levels" for each variable:
	int lev=0;
	gsl_vector * levels=gsl_vector_calloc(nofvars);
  for (i=0; i<nofvars; i++){
		lev=0;
		for (int j=0; j<p; j++){
			lev=lev+gsl_matrix_get(positions, i, j);
		}
		gsl_vector_set(levels, i, lev);		
  }

	//Rprintf("Vector number of levels:\n");
	//PrintVector(levels, nofvars);
	
	//Factors:
	//A vector of the same length as the number of variables (factors or nums)
	//with a one in the position of a vector (its corresponding row in positionsx
	//has more than one element being 1)
	gsl_vector * isfactor = gsl_vector_calloc(nofvars);
	double suma=0.0;
	  for (int i = 0; i < nofvars; i++)
	   {
			 suma=0.0;	 
			 for (int j=0; j<p; j++){
	       suma=suma+gsl_matrix_get(positions,i,j);				
			 }
			 if (suma>1){
				 gsl_vector_set(isfactor, i, 1.0);
			 }
	   }
	 
	//PrintVector(isfactor, nofvars);
	
	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);

	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2fr=0; //The number of covs in the full rank representation
	
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
	//indexfr is a full rank copy of index
	gsl_vector * indexfr = gsl_vector_calloc(p);
	
	//We fill it with the initial model
	gsl_vector_fscanf(finitM, index);
	fclose(finitM);

	//Factors:
	/*This is not used:
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
	*/
	
	//The HPM of the models visited:
	gsl_vector * HPM = gsl_vector_calloc(p);
	gsl_vector_memcpy(HPM, index);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	//this file contains the type of Bayes factor used:
	//0 is unitary; 1 is Robust; 2 is gprior; 4 is Liang and 5 is Zellner-Siow
	//fls is not yet implemented because it depends on an extra parameter
	char nfileP1[100] = "/typeofBF.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileP1);
	strcpy(nfileP1,strtmp);
	FILE * fpriorconfig = fopen(nfileP1, "r");
	//it is read and then written in an integer variable called typeofBF
	int typeofBF=0, ttt=0; //in some compilers we need to assign the value of fscanf (which is an integer)
	//correspoding to the number of input items successfully matched and assigned
	ttt = fscanf(fpriorconfig, "%d", &typeofBF);
	if (ttt != 1) {
		Rprintf("Error reading file containing the type of Bayes factor\n");
	}
	fclose(fpriorconfig);	
	//Now define the function for Bayes factors which is a pointer
	double (*BF21fun)(int, int, int, double)=NULL;
	if (typeofBF==0) BF21fun=unitBF21fun;
	if (typeofBF==1) BF21fun=RobustBF21fun;
	if (typeofBF==2) BF21fun=gBF21fun;
	if (typeofBF==4) BF21fun=LiangBF21fun;
	if (typeofBF==5) BF21fun=ZSBF21fun;

	
	
	// //////////////////////////////////////////////
	int iter=1, component=1, oldcomponent=1, newcomponent=1;
	
	//In what follows we call, for a given model Mi, PBF=BF(Mi vs M0)*Pr(Mi)
	double oldPBF=0.0, newPBF=0.0;
	//PrMg will contain the probability of each sampled model
	double PrMg=0.0;
	
    k2=(int) gsl_blas_dasum(index);
    if (k2>0){
			//Copy the sampled model (index) to a new one (indexfr) that will contain its full rank dimension copy
			//to compute Q and the Bayes factor. The real model, index, is used for the rest of purposes
			//but be CAREFUL as each time the SBSBpriorprob is called the "model" is changed
		  gsl_vector_memcpy(indexfr, index);
			PrMg=SBConstpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
	    k2fr=(int) gsl_blas_dasum(indexfr);			
      Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
	    oldPBF=BF21fun(n,k2fr+knull,knull,Q)*PrMg;
    }
    else{
		  gsl_vector_memcpy(indexfr, index);	
			PrMg=SBConstpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
			gsl_vector_set_zero(hatbetap);
      oldPBF= 1.0*PrMg;
 	  }
    
	double HPMBF=oldPBF;
	double ratio=0.0;

	//Burnin
	//Interpret old as current and new as proposal
	for (iter=1; iter<(*pBurnin+1); iter++){
		for (component=0; component<p; component++)
		{
			oldcomponent=gsl_vector_get(index, component);
			gsl_vector_set(index, component, 1-oldcomponent);
            k2=(int) gsl_blas_dasum(index);
            if (k2>0){
        			//Copy the sampled model (index) to a new one (indexfr) that will contain its full rank dimension copy
			        //to compute the Bayes factor. The real model, index, is used for the rest of purposes
						  gsl_vector_memcpy(indexfr, index);
							PrMg=SBConstpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
							k2fr=(int) gsl_blas_dasum(indexfr);			
              Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
		          newPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
            }
            else{
						  gsl_vector_memcpy(indexfr, index);	
							PrMg=SBConstpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
							gsl_vector_set_zero(hatbetap);
				      newPBF= 1.0*PrMg;
            }
						
            ratio=(oldcomponent*(oldPBF-newPBF)+newPBF)/(newPBF+oldPBF);
            newcomponent=gsl_ran_bernoulli(ran, ratio);
		      	if (newcomponent==oldcomponent){
	       	 		gsl_vector_set(index, component, newcomponent);
			      	k2=k2-1+2*oldcomponent;
			      }
		      	else
			      	oldPBF=newPBF;
	        	}
		R_CheckUserInterrupt();
	}
	
	//main loop
	//Interpret old as current and new as proposal
	for (iter=1; iter<(SAVE+1); iter++){
		cont=0;
		while (cont<nthin){
		for (component=0; component<p; component++)
		{
			oldcomponent=gsl_vector_get(index, component);
			gsl_vector_set(index, component, 1-oldcomponent);
            k2=(int) gsl_blas_dasum(index);
						if (k2>0){
						  gsl_vector_memcpy(indexfr, index);
							PrMg=SBConstpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
							k2fr=(int) gsl_blas_dasum(indexfr);			
              Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
		          newPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
            }
            else{
						  gsl_vector_memcpy(indexfr, index);	
							PrMg=SBConstpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
							gsl_vector_set_zero(hatbetap);
				      newPBF= 1.0*PrMg;
            }

              ratio=(oldcomponent*(oldPBF-newPBF)+newPBF)/(newPBF+oldPBF);
			        newcomponent=gsl_ran_bernoulli(ran, ratio);
		        	if (newcomponent==oldcomponent){
			          	gsl_vector_set(index, component, newcomponent);
				          k2=k2-1+2*oldcomponent;
		        	}
	        		else
	          			oldPBF=newPBF;
          		}
			cont++;
			R_CheckUserInterrupt();
		}
		
		/*
		Rprintf("\n\n\n");
		Rprintf("Accepted model:\n");
		PrintVector(index, p);
		Rprintf("Su dimension:%d", k2);
		Rprintf("Bg0*PrMg=%.10f\n", oldPBF);
		*/
		
		//update the summaries
		//inclusion probs:
		gsl_blas_daxpy(1.0, index, incl_prob);		
		//joint inclusion probs:
		gsl_blas_dger(1.0, index, index, joint_incl_prob);		
		//mean of betahats
		gsl_blas_daxpy(1.0, hatbetap, meanhatbetap);		
		//probability of dimensions:
		gsl_vector_set(dimension_prob,k2,gsl_vector_get(dimension_prob,k2)+1.0);
		//HPM
		if (oldPBF>HPMBF)
			{
			HPMBF=oldPBF;
			gsl_vector_memcpy(HPM, index);
			}
		//Write to the file the visited model
		my_gsl_vector_fprintf(fAllModels, index, "%f");
		//and the BF's
		fprintf(fAllBF, "%.20f\n", oldPBF/PrMg);  
      

	}
	
	//normalize:
	gsl_vector_scale(incl_prob, 1.0/SAVE);
	gsl_vector_scale(meanhatbetap, 1.0/SAVE);
	gsl_vector_scale(dimension_prob, 1.0/SAVE);
	gsl_matrix_scale(joint_incl_prob, 1.0/SAVE);
	
		
	//Write the results:
	char nfile1[100]="/LastModel";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fLastModel = fopen(strcat(nfile1,subindex), "w");
	
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
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
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
	gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
	my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
	gsl_vector_fprintf(fModels, HPM, "%f");
	gsl_vector_fprintf(fLastModel, index, "%f");
	
	
	fclose(fLastModel);
	fclose(fModels);	
	fclose(fInclusion);
	fclose(fDim);
	fclose(fJointInclusion);
	fclose(fbetahat);
	fclose(fAllModels);
	fclose(fAllBF);
	

	
	gsl_vector_free (y);
	gsl_matrix_free (X);	
	gsl_vector_free (index);
	
	gsl_vector_free(positionsx);
	gsl_vector_free(levels);
	gsl_vector_free(isfactor);
	gsl_vector_free(indexfr);
	gsl_matrix_free(positions);
	
	
	gsl_vector_free(incl_prob);
	gsl_matrix_free(joint_incl_prob);
	gsl_vector_free(dimension_prob);
	gsl_vector_free(hatbetap);
	gsl_vector_free(meanhatbetap);
	
	free(BF21fun);
		
	gsl_rng_free (ran);
		
	//Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);	
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;

}



void GibbsFConst (char *pI[], int *pn, int *pp, int *pSAVE, char *homePath[], int *pBurnin, double *time, int *pknull, int *pnthin, int *pseed)
{
	//Version where the null model is only the error term and vs is performed over the whole design
	//matrix. It keeps track of the possibility that this is used in combination with a previous
	//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
	//the original null model is only the intercept
	
	//16-5-18 version of GibbsRobust with Factors and prior probabilities obtained with a pure constant prior
  //18-12-18 
  //the type of Bayes factor is defined in a file called typeofBF.txt 
	
	void R_CheckUserInterrupt(void);
	
	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();
	
	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int SAVE=*pSAVE;
	int knull=*pknull;
	int nthin=*pnthin;
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
    gsl_rng * ran = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(ran, *pseed);
	//double info=pow(2,p);
	
	Rprintf("The full design matrix has %d columns and n=%d rows\n", p, n);
	//Rprintf("The whole problem would have %f different competing models\n", info);
	
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
	
	
	//Factors:
	//positions is a file that contains, in binary, the positions of each of the vars
	//(either factor and its levels or numerical vars) in the design matrix
	char nfileR12[100] = "/positions.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR12);
	strcpy(nfileR12,strtmp);
	FILE * fpositions = fopen(nfileR12, "r");
	int i;
	int nofvars=0; //Number of different explanatory variables (either factors or numerical)
  while((i=fgetc(fpositions))!=EOF)
    {
        if (i == '\n')
            nofvars++;
    }
		
  Rprintf("Number of explanatory variables: %d\n", nofvars);
	
	//Factors: I close and open to start the reading from the beginning
	fclose(fpositions);		
	fpositions = fopen(nfileR12, "r");
		
		
	//Factors:
	//positionsx is a file that contains, in binary, the positions of numerical vars) in the design matrix
	char nfileR13[100] = "/positionsx.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR13);
	strcpy(nfileR13,strtmp);
	FILE * fpositionsx = fopen(nfileR13, "r");
			
	
	//The initial model:
	char nfileR3[100] = "/initialmodel.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR3);
	strcpy(nfileR3,strtmp);
	FILE * finitM = fopen(nfileR3, "r");
		
	//The prior probabilities:
	//Factor:
	/*
	char nfileR4[100] = "/priorprobs.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR4);
	strcpy(nfileR4,strtmp);
	FILE * fPriorProb = fopen(nfileR4, "r");	
	*/
		
	//File that contain all visited models after burnin
	char nfile10[100]="/AllModels";
	strcpy(strtmp,home);
	strcat(strtmp,nfile10);
	strcpy(nfile10,strtmp);
	FILE * fAllModels = fopen(strcat(nfile10,subindex), "w");
	
	//File that contain all BF's of previuos models
	char nfile11[100]="/AllBF";
	strcpy(strtmp,home);
	strcat(strtmp,nfile11);
	strcpy(nfile11,strtmp);
	FILE * fAllBF = fopen(strcat(nfile11,subindex), "w");
	
	
	int cont=0; //cont is used for the thining
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);
	
	
	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	
	
  //Factors:
	gsl_vector * positionsx = gsl_vector_alloc(p);
  gsl_matrix * positions = gsl_matrix_calloc(nofvars, p);	
	gsl_vector_fscanf(fpositionsx, positionsx);
	fclose(fpositionsx);
	gsl_matrix_fscanf(fpositions, positions);
	fclose(fpositions);	
	
	//Rprintf("Matrix positions:\n");
	//PrintMatrix(positions, nofvars, p);
	
	//Rprintf("Vector positions of x:\n");
	//PrintVector(positionsx, p);
	
	//Factors:
	//A vector with the number of "levels" for each variable:
	int lev=0;
	gsl_vector * levels=gsl_vector_calloc(nofvars);
  for (i=0; i<nofvars; i++){
		lev=0;
		for (int j=0; j<p; j++){
			lev=lev+gsl_matrix_get(positions, i, j);
		}
		gsl_vector_set(levels, i, lev);		
  }

	//Rprintf("Vector number of levels:\n");
	//PrintVector(levels, nofvars);
	
	//Factors:
	//A vector of the same length as the number of variables (factors or nums)
	//with a one in the position of a vector (its corresponding row in positionsx
	//has more than one element being 1)
	gsl_vector * isfactor = gsl_vector_calloc(nofvars);
	double suma=0.0;
	  for (int i = 0; i < nofvars; i++)
	   {
			 suma=0.0;	 
			 for (int j=0; j<p; j++){
	       suma=suma+gsl_matrix_get(positions,i,j);				
			 }
			 if (suma>1){
				 gsl_vector_set(isfactor, i, 1.0);
			 }
	   }
	 
	//PrintVector(isfactor, nofvars);
	
	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the model with no covariates
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;

	SSEnull = gsl_blas_dnrm2(y);
	SSEnull = pow(SSEnull,2);

	//k2 will contain the number of covariates in each of the models visited
	int k2=0;
	int k2fr=0; //The number of covs in the full rank representation
	
	
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
	//indexfr is a full rank copy of index
	gsl_vector * indexfr = gsl_vector_calloc(p);
	
	//We fill it with the initial model
	gsl_vector_fscanf(finitM, index);
	fclose(finitM);

	//Factors:
	/*This is not used:
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
	*/
	
	//The HPM of the models visited:
	gsl_vector * HPM = gsl_vector_calloc(p);
	gsl_vector_memcpy(HPM, index);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	//this file contains the type of Bayes factor used:
	//0 is unitary; 1 is Robust; 2 is gprior; 4 is Liang and 5 is Zellner-Siow
	//fls is not yet implemented because it depends on an extra parameter
	char nfileP1[100] = "/typeofBF.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileP1);
	strcpy(nfileP1,strtmp);
	FILE * fpriorconfig = fopen(nfileP1, "r");
	//it is read and then written in an integer variable called typeofBF
	int typeofBF=0, ttt=0; //in some compilers we need to assign the value of fscanf (which is an integer)
	//correspoding to the number of input items successfully matched and assigned
	ttt = fscanf(fpriorconfig, "%d", &typeofBF);
	if (ttt != 1) {
		Rprintf("Error reading file containing the type of Bayes factor\n");
	}
	fclose(fpriorconfig);	
	//Now define the function for Bayes factors which is a pointer
	double (*BF21fun)(int, int, int, double)=NULL;
	if (typeofBF==0) BF21fun=unitBF21fun;
	if (typeofBF==1) BF21fun=RobustBF21fun;
	if (typeofBF==2) BF21fun=gBF21fun;
	if (typeofBF==4) BF21fun=LiangBF21fun;
	if (typeofBF==5) BF21fun=ZSBF21fun;

	
	
	// //////////////////////////////////////////////
	int iter=1, component=1, oldcomponent=1, newcomponent=1;
	
	//In what follows we call, for a given model Mi, PBF=BF(Mi vs M0)*Pr(Mi)
	double oldPBF=0.0, newPBF=0.0;
	//PrMg will contain the probability of each sampled model
	double PrMg=0.0;
	
    k2=(int) gsl_blas_dasum(index);
    if (k2>0){
		  gsl_vector_memcpy(indexfr, index);
			//I use the SBSBpriorprob function only to get indexfr
			PrMg=SBSBpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
			PrMg=Constpriorprob(p,k2);
	    k2fr=(int) gsl_blas_dasum(indexfr);			
      Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
      oldPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
    }
    else{
			PrMg=Constpriorprob(p,k2);
			gsl_vector_set_zero(hatbetap);
      oldPBF= 1.0*PrMg;
		}
    
	double HPMBF=oldPBF;
	double ratio=0.0;

	//Burnin
	//Interpret old as current and new as proposal
	for (iter=1; iter<(*pBurnin+1); iter++){
		for (component=0; component<p; component++)
		{
			oldcomponent=gsl_vector_get(index, component);
			gsl_vector_set(index, component, 1-oldcomponent);
            k2=(int) gsl_blas_dasum(index);
            if (k2>0){
							
						  gsl_vector_memcpy(indexfr, index);
							//I use the SBSBpriorprob function only to get indexfr
							PrMg=SBSBpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
							PrMg=Constpriorprob(p,k2);
					    k2fr=(int) gsl_blas_dasum(indexfr);			
				      Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
				      newPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
	          }
            else{
							PrMg=Constpriorprob(p,k2);
							gsl_vector_set_zero(hatbetap);
				      newPBF= 1.0*PrMg;
            }
            ratio=(oldcomponent*(oldPBF-newPBF)+newPBF)/(newPBF+oldPBF);
            newcomponent=gsl_ran_bernoulli(ran, ratio);
			if (newcomponent==oldcomponent){
				gsl_vector_set(index, component, newcomponent);
				k2=k2-1+2*oldcomponent;
			}
			else
				oldPBF=newPBF;
		}
		R_CheckUserInterrupt();
	}
	
	//main loop
	//Interpret old as current and new as proposal
	for (iter=1; iter<(SAVE+1); iter++){
		cont=0;
		while (cont<nthin){
		for (component=0; component<p; component++)
		{
			oldcomponent=gsl_vector_get(index, component);
			gsl_vector_set(index, component, 1-oldcomponent);
            k2=(int) gsl_blas_dasum(index);

            if (k2>0){
							
						  gsl_vector_memcpy(indexfr, index);
							//I use the SBSBpriorprob function only to get indexfr
							PrMg=SBSBpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
							PrMg=Constpriorprob(p,k2);
					    k2fr=(int) gsl_blas_dasum(indexfr);			
				      Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
				      newPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
	          }
            else{
							PrMg=Constpriorprob(p,k2);
							gsl_vector_set_zero(hatbetap);
				      newPBF= 1.0*PrMg;
            }

            ratio=(oldcomponent*(oldPBF-newPBF)+newPBF)/(newPBF+oldPBF);
			newcomponent=gsl_ran_bernoulli(ran, ratio);
			if (newcomponent==oldcomponent){
				gsl_vector_set(index, component, newcomponent);
				k2=k2-1+2*oldcomponent;
			}
			else
				oldPBF=newPBF;
		}
			cont++;
			R_CheckUserInterrupt();
		}
		
		//update the summaries
		//inclusion probs:
		gsl_blas_daxpy(1.0, index, incl_prob);		
		//joint inclusion probs:
		gsl_blas_dger(1.0, index, index, joint_incl_prob);		
		//mean of betahats
		gsl_blas_daxpy(1.0, hatbetap, meanhatbetap);		
		//probability of dimensions:
		gsl_vector_set(dimension_prob,k2,gsl_vector_get(dimension_prob,k2)+1.0);
		//HPM
		if (oldPBF>HPMBF)
			{
			HPMBF=oldPBF;
			gsl_vector_memcpy(HPM, index);
			}
		//Write to the file the visited model
		my_gsl_vector_fprintf(fAllModels, index, "%f");
		//and the BF's
		fprintf(fAllBF, "%.20f\n", oldPBF/PrMg); 
      

	}
	
	//normalize:
	gsl_vector_scale(incl_prob, 1.0/SAVE);
	gsl_vector_scale(meanhatbetap, 1.0/SAVE);
	gsl_vector_scale(dimension_prob, 1.0/SAVE);
	gsl_matrix_scale(joint_incl_prob, 1.0/SAVE);
	
		
	//Write the results:
	char nfile1[100]="/LastModel";
	strcpy(strtmp,home);
	strcat(strtmp,nfile1);
	strcpy(nfile1,strtmp);
	FILE * fLastModel = fopen(strcat(nfile1,subindex), "w");
	
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
	
	char nfile5[100]="/ProbDimension";
	strcpy(strtmp,home);
	strcat(strtmp,nfile5);
	strcpy(nfile5,strtmp);
	FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
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
	gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
	my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
	gsl_vector_fprintf(fModels, HPM, "%f");
	gsl_vector_fprintf(fLastModel, index, "%f");
	
	
	fclose(fLastModel);
	fclose(fModels);	
	fclose(fInclusion);
	fclose(fDim);
	fclose(fJointInclusion);
	fclose(fbetahat);
	fclose(fAllModels);
	fclose(fAllBF);

	
	gsl_vector_free (y);
	gsl_matrix_free (X);	
	gsl_vector_free (index);
	
	gsl_vector_free(positionsx);
	gsl_vector_free(levels);
	gsl_vector_free(isfactor);
	gsl_vector_free(indexfr);
	gsl_matrix_free(positions);
	
	
	gsl_vector_free(incl_prob);
	gsl_matrix_free(joint_incl_prob);
	gsl_vector_free(dimension_prob);
	gsl_vector_free(hatbetap);
	gsl_vector_free(meanhatbetap);
	
	free(BF21fun);
		
	gsl_rng_free (ran);
		
	//Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);	
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;

}


void GibbsFSB (char *pI[], int *pn, int *pp, int *pSAVE, char *homePath[], int *pBurnin, double *time, int *pknull, int *pnthin, int *pseed)
		{
			//Version where the null model is only the error term and vs is performed over the whole design
			//matrix. It keeps track of the possibility that this is used in combination with a previous
			//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
			//the original null model is only the intercept
	
//18-12-18 version of Gibbs with Factors and prior probabilities obtained with SB-SB
//the type of Bayes factor is defined in a file called typeofBF.txt 

	
			void R_CheckUserInterrupt(void);
	
			clock_t tiempo_ejec=clock();
	
			gsl_set_error_handler_off();
	
			//PARAMETERS: (R version)
			int n=*pn;
			int p=*pp;
			int SAVE=*pSAVE;
			int knull=*pknull;
			int nthin=*pnthin;
			//The files have an index attached:
			char subindex[100];
			strcpy(subindex,*pI);
			//where the files are read and writen the results:
			char home[100];
			strcpy(home,*homePath);
			//-----------
	
		    gsl_rng * ran = gsl_rng_alloc(gsl_rng_mt19937);
		    gsl_rng_set(ran, *pseed);
			//double info=pow(2,p);
	
			Rprintf("The full design matrix has %d columns and n=%d rows\n", p, n);
			//Rprintf("The whole problem would have %f different competing models\n", info);
	
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
	
	
			//Factors:
			//positions is a file that contains, in binary, the positions of each of the vars
			//(either factor and its levels or numerical vars) in the design matrix
			char nfileR12[100] = "/positions.txt";
			strcpy(strtmp,home);
			strcat(strtmp,nfileR12);
			strcpy(nfileR12,strtmp);
			FILE * fpositions = fopen(nfileR12, "r");
			int i;
			int nofvars=0; //Number of different explanatory variables (either factors or numerical)
		  while((i=fgetc(fpositions))!=EOF)
		    {
		        if (i == '\n')
		            nofvars++;
		    }
		
		  Rprintf("Number of explanatory variables: %d\n", nofvars);
	
			//Factors: I close and open to start the reading from the beginning
			fclose(fpositions);		
			fpositions = fopen(nfileR12, "r");
		
		
			//Factors:
			//positionsx is a file that contains, in binary, the positions of numerical vars) in the design matrix
			char nfileR13[100] = "/positionsx.txt";
			strcpy(strtmp,home);
			strcat(strtmp,nfileR13);
			strcpy(nfileR13,strtmp);
			FILE * fpositionsx = fopen(nfileR13, "r");
			
	
			//The initial model:
			char nfileR3[100] = "/initialmodel.txt";
			strcpy(strtmp,home);
			strcat(strtmp,nfileR3);
			strcpy(nfileR3,strtmp);
			FILE * finitM = fopen(nfileR3, "r");
		
			//The prior probabilities:
			//Factor:
			/*
			char nfileR4[100] = "/priorprobs.txt";
			strcpy(strtmp,home);
			strcat(strtmp,nfileR4);
			strcpy(nfileR4,strtmp);
			FILE * fPriorProb = fopen(nfileR4, "r");	
			*/
		
			//File that contain all visited models after burnin
			char nfile10[100]="/AllModels";
			strcpy(strtmp,home);
			strcat(strtmp,nfile10);
			strcpy(nfile10,strtmp);
			FILE * fAllModels = fopen(strcat(nfile10,subindex), "w");
	
			//File that contain all BF's of previuos models
			char nfile11[100]="/AllBF";
			strcpy(strtmp,home);
			strcat(strtmp,nfile11);
			strcpy(nfile11,strtmp);
			FILE * fAllBF = fopen(strcat(nfile11,subindex), "w");
	
	
			int cont=0; //cont is used for the thining
	
			gsl_vector * y=gsl_vector_calloc(n);
			gsl_matrix * X=gsl_matrix_calloc(n,p);
	
	
			//Put the values in the response file into the vector y	
			gsl_vector_fscanf(fResponse, y);
			fclose(fResponse);
			//and those of the design
			gsl_matrix_fscanf(fDesign, X);
			fclose(fDesign);	
	
		  //Factors:
			gsl_vector * positionsx = gsl_vector_alloc(p);
		  gsl_matrix * positions = gsl_matrix_calloc(nofvars, p);	
			gsl_vector_fscanf(fpositionsx, positionsx);
			fclose(fpositionsx);
			gsl_matrix_fscanf(fpositions, positions);
			fclose(fpositions);	
	
			//Rprintf("Matrix positions:\n");
			//PrintMatrix(positions, nofvars, p);
	
			//Rprintf("Vector positions of x:\n");
			//PrintVector(positionsx, p);
	
			//Factors:
			//A vector with the number of "levels" for each variable:
			int lev=0;
			gsl_vector * levels=gsl_vector_calloc(nofvars);
		  for (i=0; i<nofvars; i++){
				lev=0;
				for (int j=0; j<p; j++){
					lev=lev+gsl_matrix_get(positions, i, j);
				}
				gsl_vector_set(levels, i, lev);		
		  }

			//Rprintf("Vector number of levels:\n");
			//PrintVector(levels, nofvars);
	
			//Factors:
			//A vector of the same length as the number of variables (factors or nums)
			//with a one in the position of a vector (its corresponding row in positionsx
			//has more than one element being 1)
			gsl_vector * isfactor = gsl_vector_calloc(nofvars);
			double suma=0.0;
			  for (int i = 0; i < nofvars; i++)
			   {
					 suma=0.0;	 
					 for (int j=0; j<p; j++){
			       suma=suma+gsl_matrix_get(positions,i,j);				
					 }
					 if (suma>1){
						 gsl_vector_set(isfactor, i, 1.0);
					 }
			   }
	 
			//PrintVector(isfactor, nofvars);
	
			//Get the Sum of squared errors for the null model:
			//SSEnull is the sum of squared errors of the model with no covariates
			//Q=SSE/SSEnull
			double SSEnull=0.0, Q=0.0;

			SSEnull = gsl_blas_dnrm2(y);
			SSEnull = pow(SSEnull,2);

			//k2 will contain the number of covariates in each of the models visited
			int k2=0;
			int k2fr=0; //The number of covs in the full rank representation
	
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
			//indexfr is a full rank copy of index
			gsl_vector * indexfr = gsl_vector_calloc(p);
	
			//We fill it with the initial model
			gsl_vector_fscanf(finitM, index);
			fclose(finitM);

			//Factors:
			/*This is not used:
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
			*/
	
			//The HPM of the models visited:
			gsl_vector * HPM = gsl_vector_calloc(p);
			gsl_vector_memcpy(HPM, index);
	
			//hatbetap will store the mle of the k2-dimensional beta but, inserted
			//in a p-dimensional vector (with corresponding positions)
			gsl_vector * hatbetap=gsl_vector_calloc(p);
			//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
			gsl_vector * meanhatbetap=gsl_vector_calloc(p);
			
			//this file contains the type of Bayes factor used:
			//0 is unitary; 1 is Robust; 2 is gprior; 4 is Liang and 5 is Zellner-Siow
			//fls is not yet implemented because it depends on an extra parameter
			char nfileP1[100] = "/typeofBF.txt";
			strcpy(strtmp,home);
			strcat(strtmp,nfileP1);
			strcpy(nfileP1,strtmp);
			FILE * fpriorconfig = fopen(nfileP1, "r");
			//it is read and then written in an integer variable called typeofBF
			int typeofBF=0, ttt=0; //in some compilers we need to assign the value of fscanf (which is an integer)
			//correspoding to the number of input items successfully matched and assigned
			ttt = fscanf(fpriorconfig, "%d", &typeofBF);
			if (ttt != 1) {
				Rprintf("Error reading file containing the type of Bayes factor\n");
			}
			fclose(fpriorconfig);	
			//Now define the function for Bayes factors which is a pointer
			double (*BF21fun)(int, int, int, double)=NULL;
			if (typeofBF==0) BF21fun=unitBF21fun;
			if (typeofBF==1) BF21fun=RobustBF21fun;
			if (typeofBF==2) BF21fun=gBF21fun;
			if (typeofBF==4) BF21fun=LiangBF21fun;
			if (typeofBF==5) BF21fun=ZSBF21fun;

	
			// //////////////////////////////////////////////
			int iter=1, component=1, oldcomponent=1, newcomponent=1;
	
			//In what follows we call, for a given model Mi, PBF=BF(Mi vs M0)*Pr(Mi)
			double oldPBF=0.0, newPBF=0.0;
			//PrMg will contain the probability of each sampled model
			double PrMg=0.0;
	
		    k2=(int) gsl_blas_dasum(index);
		    if (k2>0){
				  gsl_vector_memcpy(indexfr, index);
					//I use the SBSBpriorprob function only to get indexfr
					PrMg=SBSBpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
					PrMg=SBpriorprob(p,k2);
			    k2fr=(int) gsl_blas_dasum(indexfr);			
		      Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
		      oldPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
		    }
		    else{
					PrMg=SBpriorprob(p,k2);
					gsl_vector_set_zero(hatbetap);
		      oldPBF= 1.0*PrMg;
				}
    
			double HPMBF=oldPBF;
			double ratio=0.0;

			//Burnin
			//Interpret old as current and new as proposal
			for (iter=1; iter<(*pBurnin+1); iter++){
				for (component=0; component<p; component++)
				{
					oldcomponent=gsl_vector_get(index, component);
					gsl_vector_set(index, component, 1-oldcomponent);
		            k2=(int) gsl_blas_dasum(index);
		            if (k2>0){
								  gsl_vector_memcpy(indexfr, index);
									//I use the SBSBpriorprob function only to get indexfr
									PrMg=SBSBpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
									PrMg=SBpriorprob(p,k2);
							    k2fr=(int) gsl_blas_dasum(indexfr);			
						      Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
						      newPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
		            }
		            else{
									PrMg=SBpriorprob(p,k2);
									gsl_vector_set_zero(hatbetap);
						      newPBF= 1.0*PrMg;
		            }
		            ratio=(oldcomponent*(oldPBF-newPBF)+newPBF)/(newPBF+oldPBF);
		            newcomponent=gsl_ran_bernoulli(ran, ratio);
					if (newcomponent==oldcomponent){
						gsl_vector_set(index, component, newcomponent);
						k2=k2-1+2*oldcomponent;
					}
					else
						oldPBF=newPBF;
				}
				R_CheckUserInterrupt();
			}
	
			//main loop
			//Interpret old as current and new as proposal
			for (iter=1; iter<(SAVE+1); iter++){
				cont=0;
				while (cont<nthin){
				for (component=0; component<p; component++)
				{
					oldcomponent=gsl_vector_get(index, component);
					gsl_vector_set(index, component, 1-oldcomponent);
		            k2=(int) gsl_blas_dasum(index);
								//Rprintf("\n\n\n");
								//Rprintf("Model:\n");
								//PrintVector(index, p);
								//Rprintf("Su dimension:%d", k2);

		            if (k2>0){
							
								  gsl_vector_memcpy(indexfr, index);
									//I use the SBSBpriorprob function only to get indexfr
									PrMg=SBSBpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
									PrMg=SBpriorprob(p,k2);
							    k2fr=(int) gsl_blas_dasum(indexfr);			
						      Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
						      newPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
		            }
		            else{
									PrMg=SBpriorprob(p,k2);
									gsl_vector_set_zero(hatbetap);
						      newPBF= 1.0*PrMg;
		            }

		            ratio=(oldcomponent*(oldPBF-newPBF)+newPBF)/(newPBF+oldPBF);
					newcomponent=gsl_ran_bernoulli(ran, ratio);
					if (newcomponent==oldcomponent){
						gsl_vector_set(index, component, newcomponent);
						k2=k2-1+2*oldcomponent;
					}
					else
						oldPBF=newPBF;
				}
					cont++;
					R_CheckUserInterrupt();
				}
		
				//update the summaries
				//inclusion probs:
				gsl_blas_daxpy(1.0, index, incl_prob);		
				//joint inclusion probs:
				gsl_blas_dger(1.0, index, index, joint_incl_prob);		
				//mean of betahats
				gsl_blas_daxpy(1.0, hatbetap, meanhatbetap);		
				//probability of dimensions:
				gsl_vector_set(dimension_prob,k2,gsl_vector_get(dimension_prob,k2)+1.0);
				//HPM
				if (oldPBF>HPMBF)
					{
					HPMBF=oldPBF;
					gsl_vector_memcpy(HPM, index);
					}
				//Write to the file the visited model
				my_gsl_vector_fprintf(fAllModels, index, "%f");
				//and the BF's
				fprintf(fAllBF, "%.20f\n", oldPBF/PrMg); 
      

			}
	
			//normalize:
			gsl_vector_scale(incl_prob, 1.0/SAVE);
			gsl_vector_scale(meanhatbetap, 1.0/SAVE);
			gsl_vector_scale(dimension_prob, 1.0/SAVE);
			gsl_matrix_scale(joint_incl_prob, 1.0/SAVE);
	
		
			//Write the results:
			char nfile1[100]="/LastModel";
			strcpy(strtmp,home);
			strcat(strtmp,nfile1);
			strcpy(nfile1,strtmp);
			FILE * fLastModel = fopen(strcat(nfile1,subindex), "w");
	
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
	
			char nfile5[100]="/ProbDimension";
			strcpy(strtmp,home);
			strcat(strtmp,nfile5);
			strcpy(nfile5,strtmp);
			FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
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
			gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
			my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
			gsl_vector_fprintf(fModels, HPM, "%f");
			gsl_vector_fprintf(fLastModel, index, "%f");
	
	
			fclose(fLastModel);
			fclose(fModels);	
			fclose(fInclusion);
			fclose(fDim);
			fclose(fJointInclusion);
			fclose(fbetahat);
			fclose(fAllModels);
			fclose(fAllBF);

	
			gsl_vector_free (y);
			gsl_matrix_free (X);	
			gsl_vector_free (index);
			
			gsl_vector_free(positionsx);
			gsl_vector_free(levels);
			gsl_vector_free(isfactor);
			gsl_vector_free(indexfr);
			gsl_matrix_free(positions);
			
	
			gsl_vector_free(incl_prob);
			gsl_matrix_free(joint_incl_prob);
			gsl_vector_free(dimension_prob);
			gsl_vector_free(hatbetap);
			gsl_vector_free(meanhatbetap);
			
			free(BF21fun);
		
			gsl_rng_free (ran);
		
			//Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);	
			*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;

		}
		
		void GibbsFConstConst (char *pI[], int *pn, int *pp, int *pSAVE, char *homePath[], int *pBurnin, double *time, int *pknull, int *pnthin, int *pseed)
		{
			//Version where the null model is only the error term and vs is performed over the whole design
			//matrix. It keeps track of the possibility that this is used in combination with a previous
			//transformation of the data on which the null ORIGINAL model had knull covariates (eg. knull=1 if
			//the original null model is only the intercept
	
			//9-5-18 version of GibbsRobust with Factors and prior probabilities obtained with Const-Const
	
      //18-12-18 
      //the type of Bayes factor is defined in a file called typeofBF.txt 
	    
			void R_CheckUserInterrupt(void);
	
			clock_t tiempo_ejec=clock();
	
			gsl_set_error_handler_off();
	
			//PARAMETERS: (R version)
			int n=*pn;
			int p=*pp;
			int SAVE=*pSAVE;
			int knull=*pknull;
			int nthin=*pnthin;
			//The files have an index attached:
			char subindex[100];
			strcpy(subindex,*pI);
			//where the files are read and writen the results:
			char home[100];
			strcpy(home,*homePath);
			//-----------
	
		    gsl_rng * ran = gsl_rng_alloc(gsl_rng_mt19937);
		    gsl_rng_set(ran, *pseed);
			//double info=pow(2,p);
	
			Rprintf("The full design matrix has %d columns and n=%d rows\n", p, n);
			//Rprintf("The whole problem would have %f different competing models\n", info);
	
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
	
	
			//Factors:
			//positions is a file that contains, in binary, the positions of each of the vars
			//(either factor and its levels or numerical vars) in the design matrix
			char nfileR12[100] = "/positions.txt";
			strcpy(strtmp,home);
			strcat(strtmp,nfileR12);
			strcpy(nfileR12,strtmp);
			FILE * fpositions = fopen(nfileR12, "r");
			int i;
			int nofvars=0; //Number of different explanatory variables (either factors or numerical)
		  while((i=fgetc(fpositions))!=EOF)
		    {
		        if (i == '\n')
		            nofvars++;
		    }
		
		  Rprintf("Number of explanatory variables: %d\n", nofvars);
	
			//Factors: I close and open to start the reading from the beginning
			fclose(fpositions);		
			fpositions = fopen(nfileR12, "r");
		
		
			//Factors:
			//positionsx is a file that contains, in binary, the positions of numerical vars) in the design matrix
			char nfileR13[100] = "/positionsx.txt";
			strcpy(strtmp,home);
			strcat(strtmp,nfileR13);
			strcpy(nfileR13,strtmp);
			FILE * fpositionsx = fopen(nfileR13, "r");
			
	
			//The initial model:
			char nfileR3[100] = "/initialmodel.txt";
			strcpy(strtmp,home);
			strcat(strtmp,nfileR3);
			strcpy(nfileR3,strtmp);
			FILE * finitM = fopen(nfileR3, "r");
		
			//The prior probabilities:
			//Factor:
			/*
			char nfileR4[100] = "/priorprobs.txt";
			strcpy(strtmp,home);
			strcat(strtmp,nfileR4);
			strcpy(nfileR4,strtmp);
			FILE * fPriorProb = fopen(nfileR4, "r");	
			*/
		
			//File that contain all visited models after burnin
			char nfile10[100]="/AllModels";
			strcpy(strtmp,home);
			strcat(strtmp,nfile10);
			strcpy(nfile10,strtmp);
			FILE * fAllModels = fopen(strcat(nfile10,subindex), "w");
	
			//File that contain all BF's of previuos models
			char nfile11[100]="/AllBF";
			strcpy(strtmp,home);
			strcat(strtmp,nfile11);
			strcpy(nfile11,strtmp);
			FILE * fAllBF = fopen(strcat(nfile11,subindex), "w");
	
	
			int cont=0; //cont is used for the thining
	
			gsl_vector * y=gsl_vector_calloc(n);
			gsl_matrix * X=gsl_matrix_calloc(n,p);
	
	
			//Put the values in the response file into the vector y	
			gsl_vector_fscanf(fResponse, y);
			fclose(fResponse);
			//and those of the design
			gsl_matrix_fscanf(fDesign, X);
			fclose(fDesign);	
	
		  //Factors:
			gsl_vector * positionsx = gsl_vector_alloc(p);
		  gsl_matrix * positions = gsl_matrix_calloc(nofvars, p);	
			gsl_vector_fscanf(fpositionsx, positionsx);
			fclose(fpositionsx);
			gsl_matrix_fscanf(fpositions, positions);
			fclose(fpositions);	
	
			//Rprintf("Matrix positions:\n");
			//PrintMatrix(positions, nofvars, p);
	
			//Rprintf("Vector positions of x:\n");
			//PrintVector(positionsx, p);
	
			//Factors:
			//A vector with the number of "levels" for each variable:
			int lev=0;
			gsl_vector * levels=gsl_vector_calloc(nofvars);
		  for (i=0; i<nofvars; i++){
				lev=0;
				for (int j=0; j<p; j++){
					lev=lev+gsl_matrix_get(positions, i, j);
				}
				gsl_vector_set(levels, i, lev);		
		  }

			//Rprintf("Vector number of levels:\n");
			//PrintVector(levels, nofvars);
	
			//Factors:
			//A vector of the same length as the number of variables (factors or nums)
			//with a one in the position of a vector (its corresponding row in positionsx
			//has more than one element being 1)
			gsl_vector * isfactor = gsl_vector_calloc(nofvars);
			double suma=0.0;
			  for (int i = 0; i < nofvars; i++)
			   {
					 suma=0.0;	 
					 for (int j=0; j<p; j++){
			       suma=suma+gsl_matrix_get(positions,i,j);				
					 }
					 if (suma>1){
						 gsl_vector_set(isfactor, i, 1.0);
					 }
			   }
	 
			//PrintVector(isfactor, nofvars);
	
			//Get the Sum of squared errors for the null model:
			//SSEnull is the sum of squared errors of the model with no covariates
			//Q=SSE/SSEnull
			double SSEnull=0.0, Q=0.0;

			SSEnull = gsl_blas_dnrm2(y);
			SSEnull = pow(SSEnull,2);

			//k2 will contain the number of covariates in each of the models visited
			int k2=0;
			int k2fr=0; //The number of covs in the full rank representation
	
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
			//indexfr is a full rank copy of index
			gsl_vector * indexfr = gsl_vector_calloc(p);
	
			//We fill it with the initial model
			gsl_vector_fscanf(finitM, index);
			fclose(finitM);

			//Factors:
			/*This is not used:
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
			*/
	
			//The HPM of the models visited:
			gsl_vector * HPM = gsl_vector_calloc(p);
			gsl_vector_memcpy(HPM, index);
	
			//hatbetap will store the mle of the k2-dimensional beta but, inserted
			//in a p-dimensional vector (with corresponding positions)
			gsl_vector * hatbetap=gsl_vector_calloc(p);
			//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
			gsl_vector * meanhatbetap=gsl_vector_calloc(p);
			
			//this file contains the type of Bayes factor used:
			//0 is unitary; 1 is Robust; 2 is gprior; 4 is Liang and 5 is Zellner-Siow
			//fls is not yet implemented because it depends on an extra parameter
			char nfileP1[100] = "/typeofBF.txt";
			strcpy(strtmp,home);
			strcat(strtmp,nfileP1);
			strcpy(nfileP1,strtmp);
			FILE * fpriorconfig = fopen(nfileP1, "r");
			//it is read and then written in an integer variable called typeofBF
			int typeofBF=0, ttt=0; //in some compilers we need to assign the value of fscanf (which is an integer)
			//correspoding to the number of input items successfully matched and assigned
			ttt = fscanf(fpriorconfig, "%d", &typeofBF);
			if (ttt != 1) {
				Rprintf("Error reading file containing the type of Bayes factor\n");
			}
			fclose(fpriorconfig);	
			//Now define the function for Bayes factors which is a pointer
			double (*BF21fun)(int, int, int, double)=NULL;
			if (typeofBF==0) BF21fun=unitBF21fun;
			if (typeofBF==1) BF21fun=RobustBF21fun;
			if (typeofBF==2) BF21fun=gBF21fun;
			if (typeofBF==4) BF21fun=LiangBF21fun;
			if (typeofBF==5) BF21fun=ZSBF21fun;

			
	
			// //////////////////////////////////////////////
			int iter=1, component=1, oldcomponent=1, newcomponent=1;
	
			//In what follows we call, for a given model Mi, PBF=BF(Mi vs M0)*Pr(Mi)
			double oldPBF=0.0, newPBF=0.0;
			//PrMg will contain the probability of each sampled model
			double PrMg=0.0;
	
		    k2=(int) gsl_blas_dasum(index);
		    if (k2>0){
					//Copy the sampled model (index) to a new one (indexfr) that will contain its full rank dimension copy
					//to compute Q and the Bayes factor. The real model, index, is used for the rest of purposes
					//but be CAREFUL as each time the ConstConstriorprob is called the "model" is changed
				  gsl_vector_memcpy(indexfr, index);
					PrMg=ConstConstpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
			    k2fr=(int) gsl_blas_dasum(indexfr);			
		      Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
			    oldPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
		    }
		    else{
				  gsl_vector_memcpy(indexfr, index);	
					PrMg=ConstConstpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
					gsl_vector_set_zero(hatbetap);
		      oldPBF= 1.0*PrMg;
				}
    
			double HPMBF=oldPBF;
			double ratio=0.0;

			//Burnin
			//Interpret old as current and new as proposal
			for (iter=1; iter<(*pBurnin+1); iter++){
				for (component=0; component<p; component++)
				{
					oldcomponent=gsl_vector_get(index, component);
					gsl_vector_set(index, component, 1-oldcomponent);
		            k2=(int) gsl_blas_dasum(index);
						    if (k2>0){
									//Copy the sampled model (index) to a new one (indexfr) that will contain its full rank dimension copy
									//to compute Q and the Bayes factor. The real model, index, is used for the rest of purposes
									//but be CAREFUL as each time the ConstConstriorprob is called the "model" is changed
								  gsl_vector_memcpy(indexfr, index);
									PrMg=ConstConstpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
							    k2fr=(int) gsl_blas_dasum(indexfr);			
						      Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
							    newPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
						    }
						    else{
								  gsl_vector_memcpy(indexfr, index);	
									PrMg=ConstConstpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
									gsl_vector_set_zero(hatbetap);
						      newPBF= 1.0*PrMg;
								}
		            ratio=(oldcomponent*(oldPBF-newPBF)+newPBF)/(newPBF+oldPBF);
		            newcomponent=gsl_ran_bernoulli(ran, ratio);
					if (newcomponent==oldcomponent){
						gsl_vector_set(index, component, newcomponent);
						k2=k2-1+2*oldcomponent;
					}
					else
						oldPBF=newPBF;
				}
				R_CheckUserInterrupt();
			}
	
			//main loop
			//Interpret old as current and new as proposal
			for (iter=1; iter<(SAVE+1); iter++){
				cont=0;
				while (cont<nthin){
				for (component=0; component<p; component++)
				{
					oldcomponent=gsl_vector_get(index, component);
					gsl_vector_set(index, component, 1-oldcomponent);
		            k2=(int) gsl_blas_dasum(index);
								//Rprintf("\n\n\n");
								//Rprintf("Model:\n");
								//PrintVector(index, p);
								//Rprintf("Su dimension:%d", k2);

						    if (k2>0){
									//Copy the sampled model (index) to a new one (indexfr) that will contain its full rank dimension copy
									//to compute Q and the Bayes factor. The real model, index, is used for the rest of purposes
									//but be CAREFUL as each time the ConstConstriorprob is called the "model" is changed
								  gsl_vector_memcpy(indexfr, index);
									PrMg=ConstConstpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
							    k2fr=(int) gsl_blas_dasum(indexfr);			
						      Q=Gibbsstatistics(p, n, SSEnull, X, y, indexfr, &k2fr, hatbetap);
							    newPBF= BF21fun(n,k2fr+knull,knull,Q)*PrMg;
						    }
						    else{
								  gsl_vector_memcpy(indexfr, index);	
									PrMg=ConstConstpriorprob(indexfr, positionsx, positions, nofvars, levels, p, isfactor);
									gsl_vector_set_zero(hatbetap);
						      newPBF= 1.0*PrMg;
								}
	
		            ratio=(oldcomponent*(oldPBF-newPBF)+newPBF)/(newPBF+oldPBF);
					newcomponent=gsl_ran_bernoulli(ran, ratio);
					if (newcomponent==oldcomponent){
						gsl_vector_set(index, component, newcomponent);
						k2=k2-1+2*oldcomponent;
					}
					else
						oldPBF=newPBF;
				}
					cont++;
					R_CheckUserInterrupt();
				}
		
				//update the summaries
				//inclusion probs:
				gsl_blas_daxpy(1.0, index, incl_prob);		
				//joint inclusion probs:
				gsl_blas_dger(1.0, index, index, joint_incl_prob);		
				//mean of betahats
				gsl_blas_daxpy(1.0, hatbetap, meanhatbetap);		
				//probability of dimensions:
				gsl_vector_set(dimension_prob,k2,gsl_vector_get(dimension_prob,k2)+1.0);
				//HPM
				if (oldPBF>HPMBF)
					{
					HPMBF=oldPBF;
					gsl_vector_memcpy(HPM, index);
					}
				//Write to the file the visited model
				my_gsl_vector_fprintf(fAllModels, index, "%f");
				//and the BF's
				fprintf(fAllBF, "%.20f\n", oldPBF/PrMg); 
      

			}
	
			//normalize:
			gsl_vector_scale(incl_prob, 1.0/SAVE);
			gsl_vector_scale(meanhatbetap, 1.0/SAVE);
			gsl_vector_scale(dimension_prob, 1.0/SAVE);
			gsl_matrix_scale(joint_incl_prob, 1.0/SAVE);
	
		
			//Write the results:
			char nfile1[100]="/LastModel";
			strcpy(strtmp,home);
			strcat(strtmp,nfile1);
			strcpy(nfile1,strtmp);
			FILE * fLastModel = fopen(strcat(nfile1,subindex), "w");
	
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
	
			char nfile5[100]="/ProbDimension";
			strcpy(strtmp,home);
			strcat(strtmp,nfile5);
			strcpy(nfile5,strtmp);
			FILE * fDim = fopen(strcat(nfile5,subindex), "w");
	
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
			gsl_vector_fprintf(fDim, dimension_prob, "%.20f");
			my_gsl_matrix_fprintf(fJointInclusion, joint_incl_prob, "%.20f");	
			gsl_vector_fprintf(fModels, HPM, "%f");
			gsl_vector_fprintf(fLastModel, index, "%f");
	
	
			fclose(fLastModel);
			fclose(fModels);	
			fclose(fInclusion);
			fclose(fDim);
			fclose(fJointInclusion);
			fclose(fbetahat);
			fclose(fAllModels);
			fclose(fAllBF);

	
			gsl_vector_free (y);
			gsl_matrix_free (X);	
			gsl_vector_free (index);
			
			gsl_vector_free(positionsx);
			gsl_vector_free(levels);
			gsl_vector_free(isfactor);
			gsl_vector_free(indexfr);
			gsl_matrix_free(positions);
			
	
			gsl_vector_free(incl_prob);
			gsl_matrix_free(joint_incl_prob);
			gsl_vector_free(dimension_prob);
			gsl_vector_free(hatbetap);
			gsl_vector_free(meanhatbetap);
			
			free(BF21fun);
		
			gsl_rng_free (ran);
		
			//Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);	
			*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;

		}

		



