//the prior probability of a model (depending on the dimension and the dimension of the full)
#include <gsl/gsl_sf_gamma.h>
#include <math.h>


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
