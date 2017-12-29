
/*Header for main.c
 * last changes Jun 11, 2014
 */

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
