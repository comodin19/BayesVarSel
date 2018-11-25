#include <R.h>
#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
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
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>


// Gonzalo Garcia-Donato and Anabel Forte, 2010
// 17-jan-14 Now with the null model not necesarily being the intercept
//The null model contains k0 parameters (k0=1 if only intercept)
//and the complex model contains k2 parameters



/* -----g-prior-------*/

// The Bayes factor obtained with the the g-prior (in favour of
// a model with k2 regressors and against the null with k0)

double gBF21fun(int n, int k2, int k0, double Q)
{

		if (k2>=n) return 1.0;
		double BF21=0.0;
    BF21 = exp(((n-k2)/2.0)*log(1.0+n)-((n-k0)/2.0)*log(1.0+n*Q));
    if (!R_FINITE(BF21)){error("A Bayes factor is infinite.");}
		return BF21;
	
}

/* -----Fernandez, Ley and Steel-prior-------*/

// The Bayes factor obtained with the the g-prior (in favour of
// a model with k2 regressors and against the null with k0)

double flsBF21fun(int p, int n, int k2, int k0, double Q)
{
    
  	if (k2>=n) return 1.0;
    int g=GSL_MAX(n, p*p);
    double BF21=0.0;
    if (!R_FINITE(BF21)){error("A Bayes factor is infinite.");}
    BF21 = exp(((n-k2)/2.0)*log(1.0+g)-((n-k0)/2.0)*log(1.0+g*Q));
    return BF21;
    
}



/*------Robust-------*/

//The Bayes Factor obtained with the robust prior of Bayarri et al 2012.


/*Structure for Parameters*/
struct parrob {
	double a;
	double b;
	double c;
	double z;
};


/*auxiliar functions for integration */

double robint_aux (double x, void *p){
    /* pointer to a structure of type par. */
	struct parrob * params=(struct parrob *)p;
    
    /*DDefine the parameters included at the structure*/
	double a=(params->a);
	double b=(params->b);
	double c=(params->c);
	double z=(params->z);
	
	/*return the argument for integration*/
	double l=pow(x,b-1.0)*pow((1.0-x),c-b-1.0)*pow((1.0-x*z),-a);
	return l;
}


/* Integrated functions the arguments will be n,k,Qi0 */
double robint (double a,double b, double c,double z){
	/*guardamos espacio de memoria para realizar la integracion, este 10000 es el que luego va en la función
	 de integracion*/
	gsl_integration_workspace * w=gsl_integration_workspace_alloc(10000);
	
	double result=0.0;
	double error=0.0;
	
	/*Ponemos los parametros en la forma que necesitamos para la funcion*/
	struct parrob  params={a,b,c,z};
	
	/*Definimos cual es la funcion que vamos a usar y le pasamos los parametros*/
	gsl_function F;
	F.function = &robint_aux;
	F.params = &params;
	
	/*integramos y guardamos el resultado en result y el error en error*/
	gsl_integration_qags(&F, 0.0,1.0, 0, 1e-9,10000,w,&result,&error);
	
	
	/*Liberamos el espacio de trabajo*/
	gsl_integration_workspace_free (w);
	
	/*devolvemos el resultado*/
	return result*gsl_sf_gamma(c)/(gsl_sf_gamma(b)*gsl_sf_gamma(c-b));
}


/* Robust Bayes Factor for main.c*/

double RobustBF21fun(int n, int k2, int k0, double Q)
{
//k2, total number of covariates in the model 
	
	if (k2>=n) return 1.0;	

	double  T1=0.0, T2=0.0, T3=0.0;
	double arg=0.0;
	double R1=0.0;
	double z=0.0;
	double rho=0.0;
  rho=pow((k2),-1.0);
  
  double k2aux=0.0;
  k2aux=k2-k0+1.0;//equivalently k2aux=ki+1
  
  double Qaux=0.0;
  Qaux=pow(Q,-1.0);
	
	// Qaux means Q_0i
	T1=log(Qaux)*((n-k0)/2.0);
	//T1=pow(Qaux,(n-k0)/2.0);
	T2=log(rho*(n+1))*(-((k2-k0))/2.0)-log(k2aux);
	//T2=pow(rho*(n+1),-((k2-k0)/2.0))*pow(k2aux,-1.0);
	
	
	//for the hypergeometric factor we distinguish whether the argument is
	//>1 or not
	
	
	arg=(1.0-Qaux)/(rho*(1.0+n));
	gsl_sf_result result;
	int STATUS=0;
	
	
	//if (arg>=-1.0)
	//{
       //STATUS=gsl_sf_hyperg_2F1_e(k2aux/2.0, (n-k0)/2.0, (k2aux/2.0)+1.0, arg,&result);
       //if(STATUS==0){
         //T3=log(result.val);
         //T3=result.val;
       //}else{
              //T3= log(robint((n-k0)/2.0,k2aux/2.0, (k2aux/2.0)+1.0, arg));
              //T3=robint((n-k0)/2.0,k2aux/2.0, (k2aux/2.0)+1.0, arg);
              //}
       //Rprintf("arg=%.20f,T3=%.20f \n",arg,T3);
	//}else {
		
		z=arg/(arg-1.0);
		STATUS=gsl_sf_hyperg_2F1_e(1,(n-k0)/2.0, (k2aux/2.0)+1.0, z,&result);
		  if (STATUS==0){ 
		      T3=((k0-n)/2)*log(1.0-arg)+log(result.val);
		      //T3=pow((1.0-arg),(k0-n)/2)*result.val;
          //Rprintf("arg=%.20f,T3=%.20f \n",arg,T3);
        }//succed
		else //gsl_hyper failed, then numerical approx of the log(2F1)
		{
		  T3=((k0-n)/2)*log(1.0-arg)+log(robint((n-k0)/2.0,1.0, (k2aux/2.0)+1.0, z));
		  //T3=pow((1.0-arg),(k0-n)/2)*robint((n-k0)/2.0,1.0, (k2aux/2.0)+1.0, z);
           //Rprintf("arg=%.20f,T3=%.20f \n",arg,T3);
        }
	//}
	
	
	R1=exp(T1+T2+T3);
	//R1=T1*T2*T3;
	
    if (!R_FINITE(R1)){error("A Bayes factor is infinite.");}
	
	return(R1);
}


/*--------LIANG-------*/

/*Structure for Parameters, for Liang and ZS*/
struct par {
	double n;
	double k_i;
	double k_0;
	double Q_i0;
};

/*auxiliar function for the argument of the integration*/
double liang_aux (double x, void *p){
	/* pointer to a structure of type par. */
    struct par * params=(struct par *)p;
	
	/*Define the parameters included at the structure */
	double n=(params->n);
	double k=(params->k_i);/*it will be k2*/
	double kk0=(params->k_0);	
	double Q=(params->Q_i0);
	
	/*Calculo el valor de la función y lo devuelvo*/
	double l=pow((1.0+x), (n-k)/2.0)*pow((1.0+Q*x), (kk0-n)/2)*(1.0/(2.0*n))*pow(1.0+x/n, -1.5);
	
	return l;
}

/*Function for integrating*/
double liang (double n, double k, double k0, double Q){
	gsl_integration_workspace * w=gsl_integration_workspace_alloc(10000);
	
	double result=0.0, error=0.0;
	
	struct par  params={n,k,k0,Q};
	
	gsl_function F;
	F.function = &liang_aux;
	F.params = &params;
	
	gsl_integration_qagiu(&F, 0, 0, 1e-9,10000,w,&result,&error);
	
	gsl_integration_workspace_free (w);
	
	return result;
}


/* Liang Bayes Factor for main.c*/
double LiangBF21fun(int n, int k2, int k0, double Q)
{
	
	  if (k2>=n) return 1.0;
    double LiangBF21=0.0;
    LiangBF21 = liang((double) n, (double) k2, (double) k0, Q);
    if (!R_FINITE(LiangBF21)){error("A Bayes factor is infinite.");}
    return LiangBF21;
    
}


/*-------JZS-------*/


/*auxiliar functions for integration */

double zell_aux (double x, void *p){
	struct par * params=(struct par *)p;/*Defino un puntero a una estructura del tipo par*/
	/*Inicializo el puntero en la dirección de memoria en la que están los parametros que le estoy
	 pasando, p obligando a que sean de tipo struct par*/
	
	/*Defino los parametros que son los que estaran en la estructura que le pasamos*/
	double n=(params->n);
	double k=(params->k_i);/*it will be k2*/
	double kk0=(params->k_0);	
	double Q=(params->Q_i0);
	
	/*Calculo el valor de la función y lo devuelvo*/
	double l=pow((1.0+x), (n-k)/2.0)*pow((1.0+Q*x), (kk0-n)/2)*pow((n/(2.0*M_PI)),0.5)*pow(x, -1.5)*exp(-n/(2.0*x));
	
	return l;
}


/*Integrated functions the arguments will be n,k,Qi0*/
double zell (double n,double k, double k0, double Q){
	/*allocate space por integration*/
	gsl_integration_workspace * w=gsl_integration_workspace_alloc(10000);
	
	double result=0.0, error=0.0;
	
	/*set parameters in the appropiate structure*/
	struct par  params={n,k,k0,Q};
	
	/*define the function and pass parameters*/
	gsl_function F;
	F.function = &zell_aux;
	F.params = &params;
	
	/*integrate and save result and error*/
	gsl_integration_qagiu(&F, 0, 0, 1e-9,10000,w,&result,&error);
	
	/*free space*/
	gsl_integration_workspace_free (w);
	
	return result;
}



/* FUNCION QUE USAREMOS EN EL main.c*/
double ZSBF21fun(int n, int k2, int k0, double Q)
{
  	if (k2>=n) return 1.0;	
    double ZSBF21=0.0;
    ZSBF21 = zell ((double) n,(double) k2, (double) k0, Q);
    if (!R_FINITE(ZSBF21)){error("A Bayes factor is infinite.");}
    return ZSBF21;
    
}
