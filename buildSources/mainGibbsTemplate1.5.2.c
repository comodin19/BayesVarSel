void GibbsgConst (char *pI[], int *pn, int *pp, int *pSAVE, char *homePath[], int *pBurnin, double *time)
{
	
	void R_CheckUserInterrupt(void);
	
	clock_t tiempo_ejec=clock();
	
	gsl_set_error_handler_off();
	
	//PARAMETERS: (R version)
	int n=*pn;
	int p=*pp;
	int SAVE=*pSAVE;
    int k0=1;//number of covariates in the most simple model
	//The files have an index attached:
	char subindex[100];
	strcpy(subindex,*pI);
	//where the files are read and writen the results:
	char home[100];
	strcpy(home,*homePath);
	//-----------
	
	const gsl_rng_type * T; 
	gsl_rng * ran;
	gsl_rng_env_setup();
	T = gsl_rng_default; 
	ran = gsl_rng_alloc (T);
	
	double info=pow(2,p-1);
	
	//Rprintf("The problem has %d variables and %d observations\n", p, n);
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
	
	//The initial model:
	char nfileR3[100] = "/initialmodel.txt";
	strcpy(strtmp,home);
	strcat(strtmp,nfileR3);
	strcpy(nfileR3,strtmp);
	FILE * finitM = fopen(nfileR3, "r");
	
	
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

	double BF21=0.0;
	//k2 will contain the number of covariates in each of the models visited
	int k2=1;
	
	//the vector with the inclusion probs (including the intercept for easy of writing):
	gsl_vector * incl_prob=gsl_vector_calloc(p);
	
	//the matrix with the joint inclusion probs:
	gsl_matrix * joint_incl_prob=gsl_matrix_calloc(p,p);
	
	//the vector with the posterior probs of each dimension (between 1:just intercept and p:full):
	gsl_vector * dimension_prob=gsl_vector_calloc(p);
	
	//index is a vector of 0 and 1's with the binary expression saying which covariates are active
	//(ie index is just the binary expression of model, an integer between 1 and N)
	gsl_vector * index = gsl_vector_calloc(p);
	
	//We fill it with the initial model
	gsl_vector_fscanf(finitM, index);
	fclose(finitM);
	
	int model=0;

	//The HPM of the models visited:
	gsl_vector * HPM = gsl_vector_calloc(p);
	gsl_vector_memcpy(HPM, index);
	
	//hatbetap will store the mle of the k2-dimensional beta but, inserted
	//in a p-dimensional vector (with corresponding positions)
	gsl_vector * hatbetap=gsl_vector_calloc(p);
	//meanhatbetap will contain the posterior mean of the betahats's (model averaged)
	gsl_vector * meanhatbetap=gsl_vector_calloc(p);
	
	// //////////////////////////////////////////////
	int iter=1, component=1, oldcomponent=1, newcomponent=1;
	
	//In what follows we call, for a given model Mi, PBF=BF(Mi vs M0)*UPr(Mi)
	//where UPr(Mi) is the unnormalized prior probability
	//this is the calculation for the initial model:
	double oldPBF=0.0, newPBF=0.0;
	
	Q=Gibbsstatistics(p, n, SSEnull, X, y, index, &k2, hatbetap);
	//printf("Q= %.20f\n", Q);
	
	//the bayes factor in favor of Mi and against M0 
	oldPBF= gBF21fun(n,k2,k0,Q)*Constpriorprob(p,k2);
	double HPMBF=oldPBF;
	double ratio=0.0;
	
	//Burnin	
	for (iter=1; iter<(*pBurnin+1); iter++){
		for (component=1; component<p; component++)
		{
			oldcomponent=gsl_vector_get(index, component);
			gsl_vector_set(index, component, 1-oldcomponent);
			Q=Gibbsstatistics(p, n, SSEnull, X, y, index, &k2, hatbetap);
			newPBF= gBF21fun(n,k2,k0,Q)*Constpriorprob(p,k2);
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
	
	
	for (iter=1; iter<(SAVE+1); iter++){
		for (component=1; component<p; component++)
		{
			oldcomponent=gsl_vector_get(index, component);
			gsl_vector_set(index, component, 1-oldcomponent);
			Q=Gibbsstatistics(p, n, SSEnull, X, y, index, &k2, hatbetap);
			newPBF= gBF21fun(n,k2,k0,Q)*Constpriorprob(p,k2);
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
		
		//update the summaries
		//inclusion probs:
		gsl_blas_daxpy(1.0, index, incl_prob);		
		//joint inclusion probs:
		gsl_blas_dger(1.0, index, index, joint_incl_prob);		
		//mean of betahats
		gsl_blas_daxpy(1.0, hatbetap, meanhatbetap);		
		//probability of dimensions:
		gsl_vector_set(dimension_prob,k2-1,gsl_vector_get(dimension_prob,k2-1)+1.0);
		//HPM
		if (newPBF>HPMBF)
			{
			HPMBF=newPBF;
			gsl_vector_memcpy(HPM, index);
			}
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
	
	gsl_vector_free (y);
	gsl_matrix_free (X);	
	gsl_vector_free (index);
	gsl_vector_free (y0);
	
	
	gsl_vector_free(incl_prob);
	gsl_matrix_free(joint_incl_prob);
	gsl_vector_free(dimension_prob);
	gsl_vector_free(hatbetap);
	gsl_vector_free(meanhatbetap);
		
	gsl_rng_free (ran);
		
	//Rprintf("Tiempo de ejecucion: %f seg\n", (double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC);	
	*time=(double) (clock()-tiempo_ejec)/CLOCKS_PER_SEC;

}
