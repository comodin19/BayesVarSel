/*Header for mainGibbs.c
 * last changes Oct 10, 2014
 */

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
//#include "priorprob.c"
