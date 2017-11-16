void recompute(gsl_vector * v, gsl_vector * w, int N);
double statistics(int model, int p, int n, double SSEnull, gsl_matrix * X, 
				  gsl_vector * y, gsl_vector * index, int *k2,
				  gsl_vector *hatbetap);
double Constmainalgebraics(double BF21, int p, gsl_vector * index, double *NormConstant, 
					  double *NormConstantPrior, gsl_vector *incl_prob, 
					  gsl_vector *dimension_prob,
					  int k2, gsl_vector *hatbetap, gsl_vector *meanhatbetap, 
					  gsl_matrix *joint_incl_prob);
double SBmainalgebraics(double BF21, int p, gsl_vector * index, double *NormConstant, 
						   double *NormConstantPrior, gsl_vector *incl_prob, 
						   gsl_vector *dimension_prob,
						   int k2, gsl_vector *hatbetap, gsl_vector *meanhatbetap, 
						   gsl_matrix *joint_incl_prob);
double Usermainalgebraics(gsl_vector *priorvector, double BF21, int p, gsl_vector * index, double *NormConstant, 
						double *NormConstantPrior, gsl_vector *incl_prob, 
						gsl_vector *dimension_prob,
						int k2, gsl_vector *hatbetap, gsl_vector *meanhatbetap, 
						gsl_matrix *joint_incl_prob);
int my_gsl_matrix_fprintf(FILE *stream,gsl_matrix *m,char *fmt);
int my_gsl_vector_fprintf(FILE *stream,gsl_vector *v,char *fmt);
int gsl_vector_mifrac(gsl_vector *v, const double x);

