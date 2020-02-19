#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

/* I wrote this code to give Jaki the isentropes from our CP EoS. 
 * I first spline the 2d 
 * surface of s/nB as a function of T and muB. Then, I loop over T,
 * and for each value I root find the value of muB corresponding to 
 * the isentrope. */

int main(int argc, char *argv[]){


	//freopen("/dev/null", "w", stderr);

	char file_snB_in[100];
	/* argv[1] is the s/nB file */
	strcpy(file_snB_in, argv[1]);

	double T_lo = 30, T_hi = 800;
	int Ti;

	/* argv [2] is the value of s/nB*/
	double snB_val = atof(argv[2]);
        
	/* we read through the s/nB matrix to determine the size */
	int nT, nmu;
	double T_in, mu_in, val_in;
	double T1, T2, mu1, mu2;
	T1 = T2 = mu1 = mu2 = - 1.0;
	FILE * file_in = fopen(file_snB_in,"r");
	for (nT = nmu = 0; fscanf(file_in,"%lf %lf %lf\n", &mu_in, &T_in, &val_in) != EOF; ){
		T1 = T_in;  mu1 = mu_in;

		if (T1 != T2) nT += 1;
		if (mu1 != mu2) nmu += 1;

		T2 = T1;  mu2 = mu1;
	}
	fclose(file_in);
	nmu = nmu/nT;

	fprintf(stderr,"we have %d values of T and %d values of mu\n", nT, nmu);

	/* we set up the spline(s) */
	/* for snB */
	const gsl_interp2d_type * T_spline_snB = gsl_interp2d_bicubic; /* choose interpolation type */
	double xa_snB[nmu]; /* vector of x */
	double ya_snB[nT]; /* vector of y */
	const size_t nx_snB = sizeof(xa_snB) / sizeof(double); /* x grid points */
	const size_t ny_snB = sizeof(ya_snB) / sizeof(double); /* y grid points */
	double *za_snB = malloc(nx_snB * ny_snB * sizeof(double));
	gsl_spline2d *spline_snB = gsl_spline2d_alloc(T_spline_snB, nx_snB, ny_snB); /* allocate the spline object */
	gsl_interp_accel *xacc_snB = gsl_interp_accel_alloc(); /* allocate the accelerators */
	gsl_interp_accel *yacc_snB = gsl_interp_accel_alloc();


	fprintf(stderr,"Now open the files to store.\n");

	int mu_idx, T_idx;
	/* we loop through the matrix and save it this time into za_snB */
	file_in = fopen(file_snB_in,"r");
	for (; fscanf(file_in,"%lf %lf %lf\n", &mu_in, &T_in, &val_in) != EOF; ){
		mu_idx = (int) (mu_in-1); 
		T_idx = (int) (T_in - 30);
		xa_snB[mu_idx] = mu_in;
		ya_snB[T_idx] = T_in;
		gsl_spline2d_set(spline_snB, za_snB, mu_in-1, T_in-30, val_in);  
		//printf("%d %f %d %f\n", mu_idx, xa_snB[(int) mu_in], T_idx, ya_snB[(int) T_in-30] );
	}
	fclose(file_in);
	
	fprintf(stderr,"actual splining now\n");	

	/* initialize interpolation */
	gsl_spline2d_init(spline_snB, xa_snB, ya_snB, za_snB, nx_snB, ny_snB);

	fprintf(stderr,"did we make it here?\n");
	
	/* now we define the root finder */
	fprintf(stderr,"Define root finder.\n");
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fdfsolver_type *T_root;
	gsl_root_fdfsolver *s;
	double x, x0;
	
	double y_fixed; /* will be looped over */

	/* define parameters and functions */
	struct snB_params{
		double a;
	};
	struct snB_params par = {snB_val};

	/* splined function for root finding */
	double splined_func (double x, void *params){
		struct snB_params *p = (struct snB_params *) params;
	
		double a = p->a;
		double func = gsl_spline2d_eval(spline_snB, x, y_fixed, xacc_snB, yacc_snB);
	
		return func - a;
	}
	/* derivative of spline function wrt x */
	double splined_func_deriv (double x, void *params){
		struct snB_params *p = (struct snB_params *) params;
	
		double a = p->a;
	
		double func_deriv = gsl_spline2d_eval_deriv_x(spline_snB, x, y_fixed, xacc_snB, yacc_snB);
	
		return func_deriv;
	}
	/* funcion and derivative together */
	void splined_func_fdf (double x, void *params, double *y, double *dy){
		struct snB_params *p = (struct snB_params *) params;
	
		double a = p->a;
	
		double func = gsl_spline2d_eval(spline_snB, x, y_fixed, xacc_snB, yacc_snB) - a;
		double func_deriv = gsl_spline2d_eval_deriv_x(spline_snB, x, y_fixed, xacc_snB, yacc_snB);
	
		*y = func;
		*dy = func_deriv;
	}

	/* assign to gsl_function_fdf*/
	gsl_function_fdf FDF;
	
	fprintf(stderr,"Assign function and parameters.\n");
	FDF.f = &splined_func;
	FDF.df = &splined_func_deriv;
	FDF.fdf = &splined_func_fdf;
	FDF.params = &par;
	
	fprintf(stderr,"Define root finding method.\n");
	T_root = gsl_root_fdfsolver_newton;
	s = gsl_root_fdfsolver_alloc (T_root);
	
	fprintf (stderr,"using %s method\n",gsl_root_fdfsolver_name (s));
	
	double val1, val2;
	/* loop over temperatures */
	for (Ti = 50; Ti < 800; Ti++ ){
		y_fixed = Ti;

		fprintf(stderr,"T is %d MeV\n",Ti);
		/* check if solution exists*/
		val1=splined_func(1,&par);
		val2=splined_func(450,&par);

		fprintf(stderr,"s/nB at mu=1 MeV: %lf, s/nB at mu=450 MeV: %lf \n",val1,val2);

		if (val1*val2 > 0){
			fprintf(stderr,"Skipping this T value\n");
			continue;
		}	

		//for (i=1;i<=450;i++) printf("%lf %lf\n",splined_func());

		iter = 0;
		gsl_root_fdfsolver_set (s, &FDF, 1);
		do{
			iter++;
			status = gsl_root_fdfsolver_iterate (s);
			x0 = x;
			x = gsl_root_fdfsolver_root (s);
			status = gsl_root_test_delta (x, x0, 0, 0.01);
			
			if (status == GSL_SUCCESS) fprintf (stderr,"Converged:\n");
			
			fprintf (stderr,"%5d %3.7f\n", iter, x);
			
		} while (status == GSL_CONTINUE && iter < max_iter);
		printf("%d\t%4.12f\n",Ti,x);
	}


	
	gsl_spline2d_free(spline_snB);
	gsl_interp_accel_free(xacc_snB);
	gsl_interp_accel_free(yacc_snB);
	free(za_snB);
	
	return 0;
}
