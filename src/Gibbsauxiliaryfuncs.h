double Gibbsstatistics(int p, int n, double SSEnull, gsl_matrix * X, 
					   gsl_vector * y, gsl_vector * index, int *k2,
					   gsl_vector *hatbetap);
void PrintMatrix(gsl_matrix * A, int nrow, int ncol);
void PrintVector(gsl_vector * v, int n);
void PrintIntMatrix (gsl_matrix_int * A, int nrow, int ncol);
void PrintIntVector (gsl_vector_int * v, int n);


