// Modified by Krisztian Sebestyen
// - retreive optimization information after the return of dtron()
// - added parameters: 
// int* convergence, int trace


// Package: kernlab
// Version: 0.9-14
// Date: 2011-11-01
// Title: Kernel-based Machine Learning Lab
// Author: Alexandros Karatzoglou, Alex Smola, Kurt Hornik
// Maintainer: Alexandros Karatzoglou <alexis@ci.tuwien.ac.at>
// Description: Kernel-based machine learning methods for classification,
        // regression, clustering, novelty detection, quantile regression
        // and dimensionality reduction.  Among other methods kernlab
        // includes Support Vector Machines, Spectral Clustering, Kernel
        // PCA and a QP solver.
// Depends: R (>= 2.10), methods
// LazyLoad: Yes
// License: GPL-2
// Packaged: 2011-11-02 18:09:44 UTC; hornik
// Repository: CRAN
// Date/Publication: 2011-11-02 18:36:27



#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <R_ext/BLAS.h> // includes <R_ext/RS.h> /* for F77_... */
#include <R.h>
#define PRINTF Rprintf

#ifndef FCONE
# define FCONE
#endif
extern void *xmalloc(size_t);
extern double mymin(double, double);
extern double mymax(double, double);
/* LEVEL 1 BLAS */
/*extern double dnrm2_(int *, double *, int *);*/
/*extern double ddot_(int *, double *, int *, double *, int *);*/
/* LEVEL 2 BLAS */
/*extern int dsymv_(char *, int *, double *, double *, int *, double *, int *, double *, double *, int *);*/
/* MINPACK 2 */
extern double dgpnrm(int, double *, double *, double *, double *);
extern void dspcg(int, double *, double *, double *, double *, double *, double, double, double *, int *);

// Krisz and Youyi
// In the kernlab package (and bsvm software) copy of dtron.c, there was an extra double* for 'wa' passed to dcauchy 
// but dcauchy() lacks it in its fn declaration. Apparently it does not affect compilation
// The extra double has no effect on numerical results, it is removed here
extern void dcauchy(int, double *, double *, double *, double *, double *, double, double *, double *, int);


double dtron(int n, double *x, double *xl, double *xu, double gtol, double frtol, double fatol, double fmin, int maxfev,double cgtol, 
     //new parameters
     int* convergence,int trace,
     int (*ugradhesptr)(int, double *, double *, double **), int (*ufvptr)(int, double *, double *)
) 
{
/*
c     *********
c
c     Subroutine dtron
c
c     The optimization problem of BSVM is a bound-constrained quadratic
c     optimization problem and its Hessian matrix is positive semidefinite. 
c     We modified the optimization solver TRON by Chih-Jen Lin and
c     Jorge More' into this version which is suitable for this
c     special case.
c
c     This subroutine implements a trust region Newton method for the
c     solution of large bound-constrained quadratic optimization problems
c
c           min { f(x)=0.5*x'*A*x + g0'*x : xl <= x <= xu }
c
c     where the Hessian matrix A is dense and positive semidefinite. The
c     user must define functions which evaluate the function, gradient, 
c     and the Hessian matrix.
c
c     The user must choose an initial approximation x to the minimizer,
c     lower bounds, upper bounds, quadratic terms, linear terms, and
c     constants about termination criterion.
c
c	parameters:
c
c       n is an integer variable.
c         On entry n is the number of variables.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x specifies the vector x.
c         On exit x is the final minimizer.
c
c       xl is a double precision array of dimension n.
c         On entry xl is the vector of lower bounds.
c         On exit xl is unchanged.
c
c       xu is a double precision array of dimension n.
c         On entry xu is the vector of upper bounds.
c         On exit xu is unchanged.
c
c       gtol is a double precision variable.
c         On entry gtol specifies the relative error of the projected 
c            gradient.
c         On exit gtol is unchanged. 
c
c       frtol is a double precision variable.
c         On entry frtol specifies the relative error desired in the
c            function. Convergence occurs if the estimate of the
c            relative error between f(x) and f(xsol), where xsol
c            is a local minimizer, is less than frtol.
c         On exit frtol is unchanged.
c
c       fatol is a double precision variable.
c         On entry fatol specifies the absolute error desired in the
c            function. Convergence occurs if the estimate of the
c            absolute error between f(x) and f(xsol), where xsol
c            is a local minimizer, is less than fatol.
c         On exit fatol is unchanged.
c
c       fmin is a double precision variable.
c         On entry fmin specifies a lower bound for the function.
c            The subroutine exits with a warning if f < fmin.
c         On exit fmin is unchanged.
c
c       maxfev is an integer variable.
c         On entry maxfev specifies the limit of function evaluations.
c         On exit maxfev is unchanged.
c
c       cgtol is a double precision variable.
c         On entry gqttol specifies the convergence criteria for
c            subproblems.
c         On exit gqttol is unchanged.
c
c     **********
*/

/*
// additional parameters
int* convergence  - 0 for failure, 1 if success
int trace
ugradhesptr: a function pointer to update gradient and hessian
ufvptr: a function pointer to update function value
*/


	/* Parameters for updating the iterates. */
	double eta0 = 1e-4, eta1 = 0.25, eta2 = 0.75;
	
	/* Parameters for updating the trust region size delta. */
	double sigma1 = 0.25, sigma2 = 0.5, sigma3 = 4;

	double p5 = 0.5, one = 1;
	double gnorm, gnorm0, delta, snorm;
	double alphac = 1.0, alpha, f, fc, prered, actred, gs;
	int search = 1, iter = 1, info, inc = 1;	
	double *xc = (double *) xmalloc(sizeof(double)*n);
	double *s = (double *) xmalloc(sizeof(double)*n);
	double *wa = (double *) xmalloc(sizeof(double)*n); 	
	double *g = (double *) xmalloc(sizeof(double)*n);
	double *A = NULL;

    *convergence = 0;
	(*ugradhesptr)(n, x, g, &A);
	(*ufvptr)(n, x, &f);	
	if(trace > 2) PRINTF(" dim: %d g: %f %f ... A: %f %f %f %f ... f: %f\n", n, g[0], g[1], A[0], A[1], A[2], A[3], f);
	gnorm0 = F77_CALL(dnrm2)(&n, g, &inc);
	delta = 1000.0*gnorm0;
	gnorm = dgpnrm(n, x, xl, xu, g);
	if (gnorm <= gtol*gnorm0)
	{
		
		if(trace>1)PRINTF("CONVERGENCE: GTOL TEST SATISFIED\n");
		
		search = 0;        
        *convergence = 1;
	}
    
    int cnt=0;
	while (search)
	{ 

		if(trace > 1) PRINTF("\n");
		
		/* Save the best function value and the best x. */   
		fc = f;
		memcpy(xc, x, sizeof(double)*n);
		

		/* Compute the Cauchy step and store in s. */	
		// Youyi: the Cauchy step has an undesirable property: when the dimensions are permutated, i.e. x1 x2 x3 -> x2 x1 x3, the resulting s differ
		dcauchy(n, x, xl, xu, A, g, delta, &alphac, s, 0 > trace - 1 ? 0 : trace - 1);
		if(trace >= 3) {
             PRINTF("x "); for (int i_s=0; i_s<n; i_s++) PRINTF("%f ", x[i_s]); PRINTF("\n");		
             PRINTF("s "); for (int i_s=0; i_s<n; i_s++) PRINTF("%f ", s[i_s]); PRINTF("\n");		
        }		
		
        /* Compute the projected Newton step. */		
		dspcg(n, x, xl, xu, A, g, delta, cgtol, s, &info);
		if(trace >= 3) {
             PRINTF("x "); for (int i_s=0; i_s<n; i_s++) PRINTF("%f ", x[i_s]); PRINTF("\n");		
             PRINTF("s "); for (int i_s=0; i_s<n; i_s++) PRINTF("%f ", s[i_s]); PRINTF("\n");		
        }		
		
        cnt=(*ufvptr)(n, x, &f);
		if(trace > 1) PRINTF("function evaluation # %i\n", cnt);
        if (cnt > maxfev)
		{
			
			PRINTF("ERROR: NFEV > MAXFEV\n");
			
			search = 0;
 
            *convergence = 0;
             
			continue;
		}

		/* Compute the predicted reduction. */
		memcpy(wa, g, sizeof(double)*n);
//		if(trace > 1) {
//             PRINTF("wa "); for (int i_s=0; i_s<n; i_s++) PRINTF("%f ", wa[i_s]); PRINTF("\n");		
//        }		
		F77_CALL(dsymv)("U", &n, &p5, A, &n, s, &inc, &one, wa, &inc FCONE);
//		if(trace > 1) {
//             PRINTF("x "); for (int i_s=0; i_s<n; i_s++) PRINTF("%f ", x[i_s]); PRINTF("\n");		
//             PRINTF("s "); for (int i_s=0; i_s<n; i_s++) PRINTF("%f ", s[i_s]); PRINTF("\n");		
//             PRINTF("wa "); for (int i_s=0; i_s<n; i_s++) PRINTF("%f ", wa[i_s]); PRINTF("\n");		
//        }		
		prered = -F77_CALL(ddot)(&n, s, &inc, wa, &inc);
		if(trace > 1) PRINTF("fc:%f f:%f prered:%f ratio:%f\n", fc, f, prered, (fc-f)/prered);
                        
		/* Compute the actual reduction. */
        actred = fc - f;                
                                    
		/* On the first iteration, adjust the initial step bound. */
		snorm = F77_CALL(dnrm2)(&n, s, &inc);
		if (iter == 1)
			delta = mymin(delta, snorm);

		/* Compute prediction alpha*snorm of the step. */
		gs = F77_CALL(ddot)(&n, g, &inc, s, &inc);

		if (f - fc - gs <= 0.0)
			alpha = sigma3;
		else
			alpha = mymax(sigma1, -0.5*(gs/(f - fc - gs)));

   		if(trace > 2)PRINTF("TRON: update trust region\n");
   		
		/* Update the trust region bound according to the ratio
		of actual to predicted reduction. */
		if (actred < eta0*prered){

            if(trace > 2)PRINTF("TRON: reduce delta for unsuccessfull step\n");

			/* Reduce delta. Step is not successful. */
			delta = mymin(mymax(alpha, sigma1)*snorm, sigma2*delta);
		}else 
		{
			if (actred < eta1*prered){
                if(trace > 2)PRINTF("TRON: reduce delta for not suff. successfull step\n");

				/* Reduce delta. Step is not sufficiently successful. */
				delta = mymax(sigma1*delta, mymin(alpha*snorm, sigma2*delta));
			} else if (actred < eta2*prered){
                if(trace > 2)PRINTF("TRON: The ratio of actual to predicted reduction is in the interval (eta1,eta2). We are allowed to either increase or decrease delta.\n");
				
					/* The ratio of actual to predicted reduction is in
					the interval (eta1,eta2). We are allowed to either
					increase or decrease delta. */
					delta = mymax(sigma1*delta, mymin(alpha*snorm, sigma3*delta));
            } else{
                if(trace > 2)PRINTF("TRON: The ratio of actual to predicted reduction exceeds eta2. Do not decrease delta.\n");
				
					/* The ratio of actual to predicted reduction exceeds eta2.
					Do not decrease delta. */
					delta = mymax(delta, mymin(alpha*snorm, sigma3*delta));
            }
		}

		/* Update the iterate. */
		if (actred > eta0*prered) 
		{
		
			/* Successful iterate. */
			iter++;

            if(trace > 2)PRINTF("TRON: Successful iterate.\n");
            
            
			(*ugradhesptr)(n, x, g, &A);
			gnorm = dgpnrm(n, x, xl, xu, g);	
            if(trace > 2) PRINTF("gnorm = %g, gtol*gnorm0 = %g \n", gnorm, gtol*gnorm0);
			if (gnorm <= gtol*gnorm0)
        		{
				
				if(trace>1) PRINTF("CONVERGENCE: GTOL TEST SATISFIED\n");
				
				search = 0;
                *convergence = 1;
                continue;
			} else if(trace>1) PRINTF("CONVERGENCE: GTOL TEST not SATISFIED");
		}
		else
		{
            if(trace > 2)PRINTF("TRON: Unsuccessful iterate.\n");

			/* Unsuccessful iterate. */
			memcpy(x, xc, sizeof(double)*n);
			f = fc;
		}
	
		/* Test for convergence */
		if (f < fmin)
		{
			if(trace)PRINTF("WARNING: F .LT. FMIN\n");
			search = 0; /* warning */
            *convergence = 1;
			continue;
		}
            if(trace > 2) PRINTF("\nfabs(actred) = %g, prered = %g, fatol = %g\n", fabs(actred), prered, fatol);
		if (fabs(actred) <= fatol && prered <= fatol)
		{
			if(trace>1)PRINTF(", FATOL TEST SATISFIED\n");
            *convergence = 1;
			search = 0;
			continue;
		} else if(trace>1) PRINTF(", FATOL TEST not SATISFIED");
        if(trace > 2) PRINTF("\n(actred) = %g, prered = %g, frtol = %g, fabs(f) = %g\n", (actred), prered, frtol, fabs(f));
		if (fabs(actred) <= frtol*fabs(f) && prered <= frtol*fabs(f))
		{
			
			if(trace>1) PRINTF(", FRTOL TEST SATISFIED\n");		
			*convergence = 1;
			search = 0;
			continue;
		} else if (trace>1) PRINTF(", FRTOL TEST not SATISFIED\n");		
	}	// end while loop

	free(g);
	free(xc);
	free(s);
	free(wa);

    if (trace) PRINTF("f: %.10f, iters %i\n", f, cnt);            

	return f;
}//end while
