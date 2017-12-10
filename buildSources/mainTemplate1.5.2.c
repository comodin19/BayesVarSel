void gConst (char *pI[], int *pn, int *pp, int *pSAVE, int *pinicio, int *pfinal, char *homePath[], double *time)
{
	void R_CheckUserInterrupt(void);

	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();

	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
		
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
	
	double info=pow(2,p-1);

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
	
	int i=1;
	
	gsl_vector * y=gsl_vector_calloc(n);
	gsl_matrix * X=gsl_matrix_calloc(n,p);

	//Put the values in the response file into the vector y	
	gsl_vector_fscanf(fResponse, y);
	fclose(fResponse);
	//and those of the design
	gsl_matrix_fscanf(fDesign, X);
	fclose(fDesign);	

	//Get the Sum of squared errors for the null model:
	//SSEnull is the sum of squared errors of the model with only intercept
	//Q=SSE/SSEnull
	double SSEnull=0.0, Q=0.0;
	double mean=0.0;
	for (i=0; i<n; i++)
		{
		mean=mean+gsl_vector_get(y,i);
		}
	mean=mean/n;
		
	gsl_vector * y0 = gsl_vector_calloc(n);
	gsl_vector_memcpy(y0,y);
	gsl_vector_add_constant(y0,-1.0*mean);
	SSEnull = gsl_blas_dnrm2(y0);
	SSEnull = pow(SSEnull,2);
	double NormConstant=0.0,NormConstantPrior=0.0;
	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=1;

	double unnormPostProb=0.0;

	//the vector with the inclusion probs (including the intercept for easy of writing):
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension (between 1:just intercept and p:full):
	gsl_vector * dimension_prob=gsl_vector_calloc(p);

	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);

	//binmodel is an auxiliary covariate
	int model=0;

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
	int dimensionNull=1;
	NormConstant=1.0*Constpriorprob(p,dimensionNull);
	NormConstantPrior+=Constpriorprob(p,dimensionNull);

	gsl_vector_set(dimension_prob, 0, 1.0*Constpriorprob(p,dimensionNull));
	gsl_vector_set(incl_prob, 0, 1.0*Constpriorprob(p,dimensionNull));
	
	gsl_vector_set(Who_Max_SAVE, 0, 0);
	gsl_vector_set(Max_SAVE_BF, 0, 1.0*Constpriorprob(p,dimensionNull));
    
    
    int k0=1;/*The number of covariates in the most simple model... 1 means the most simple model being the NULL (just containing the intercept)*/
	// //////////////////////////////////////////////
	//LOOP I: First we fill the vectors with the first values of the Bayes factor
	for (model=StartAt; model<(SAVE+StartAt-1); model++)
			{
		Q=statistics(model, p, n, SSEnull, X, y, index, &k2, hatbetap);
        /*Returns Q_i0, index (binary expression of the model), and k2(number of covariates in the model including the intercept)*/
                
		//printf("Q= %.20f\n", Q);

		//the bayes factor in favor of modeli and against M0 
		BF21= gBF21fun(n,k2,k0,Q);
		//printf("The Bayes factor is: %.20f\n", BF21);

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
		BF21= gBF21fun(n,k2,k0,Q);
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
		gsl_vector_free (y0);
	
	
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

