#include <R.h>
#define PRINTF Rprintf

#include <errno.h>

// TRON declarations
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <R_ext/BLAS.h>
/* LEVEL 1 BLAS */
/*extern double ddot_(int *, double *, int *, double *, int *); */
/* LEVEL 2 BLAS */
/*extern int dsymv_(char *, int *, double *, double *, int *, double *, int *, double *, double *, int *);*/
/* MINPACK 2 */

#include "matrix.h" //crossprod() and tcrossprod() as R's source code
double dtron(int n, double *x, double *xl, double *xu, double gtol, double frtol, double fatol, double fmin, int maxfev, double cgtol,int* convergence, int verbose,
    int (*ugradhesptr)(int, double *, double *, double **), int (*ufvptr)(int, double *, double *) );


// MACROS for saud_dca
#define	S1(x)	log(1+exp(-b*x))/zeta
#define	S1_deriv(x) -b/(1+exp(b*x))/zeta
#define S1_dd(x) b*b/(1+exp(b*x))/(1+exp(-b*x))/zeta




static int nfev_sauc_dca; // counter for number of function evaluations
static double * eta; // hold intermediate variables used in gradient and hessian updates
static double * HESS_sauc_dca; // hessian

// data input
static int n1n2;
static double lambda;
static double * X_diff;
static double * S2_deriv_ij;
static double * K; //p x p, p = n for nonlinear kernel 
static double zeta;
static double b;


// x is read-only length p parameter vector
// g is gradient at x, H is hessian at x
int ugradhes_sauc_dca(int p, double *x, double *g, double **H)
{
    int ONE=1; 
    
    *H = HESS_sauc_dca;
    tcrossprod(x,&ONE,&p,X_diff,&n1n2,&p,eta);	

    /* Hessian */
	memset(HESS_sauc_dca, 0, sizeof(double)*p*p);	
    double t;
    for (int i=0; i<n1n2; i++) {
        t=S1_dd(eta[i]);
        //  H: = S1_dd(eta[i], zeta) * trans(X.diff[i,]) %*% X.diff[i,]) + H
    	F77_CALL(dger)(&p, &p, &t, X_diff+i, &n1n2, X_diff+i, &n1n2, HESS_sauc_dca, &p);
     }
     
     for (int i=0; i<p; i++) 
         for (int j=0; j<p; j++) 
             HESS_sauc_dca[i*p+j] += lambda * K[i*p + j]; 
    
    
    /* Gradient */
     
    // S1.deriv(eta, zeta) - S2.deriv.ij
    // use eta to store this value to save time/space
    for (int i=0; i<n1n2; i++) 
        eta[i] = S1_deriv(eta[i]) - S2_deriv_ij[i];
    	
    //  g := t(X.diff) %*% (S1.deriv(eta, zeta) - S2.deriv.ij) + 0
	crossprod(X_diff,&n1n2,&p,eta,&n1n2,&ONE,g);	

    //  g += lambda * beta or g += lambda * K %*% alpha; for linear kernel, K is pxp identity matrix, otherwise it is the p x p kernel matrix
	crossprod(x, &p,&ONE,K, &p, &p, eta);
	for (int i=0; i<p; i++) g[i] += lambda*eta[i];
	
	return 0;

}

int ufv_sauc_dca(int p, double *x, double *f)
{
    *f=0.0;
    int ONE=1;
    
    // loss    
    // eta := X.diff %*% x + 0
	tcrossprod(x,&ONE,&p,X_diff,&n1n2,&p,eta);  // t(x) %*% X, X _{n1n2 x p} stored column-major 
    //sum(S1(eta, zeta) - eta * S2.deriv.ij)
    for (int i=0; i<n1n2; i++) 
        *f += S1(eta[i]) - eta[i] * S2_deriv_ij[i];
    
    // penalty    
	//*f += lambda * (||beta||_2)^2 or *f += lambda * t(alpha) %*% K %*% alpha
	crossprod(x, &p, &ONE,K, &p, &p, eta); //eta = K %*% x	
	double xKx = 0.0;
	for (int i=0; i<p; i++) 
        xKx += eta[i] * x[i]; // += t(x) %*% K %*% x        
    *f += .5 * lambda * xKx;
	
    return ++nfev_sauc_dca;
}

void dcsauc_tron(
	double* _K, double* _X_diff, double* _S2_deriv_ij, 
	double* _lambda, int* _p, int* _n1n2, double * _zeta, double * _b,
	int* _maxfev, double* _gtol, double* _frtol, int* verbose, int *_q, //_q is not used, but is here to be consisent with sauc_dca_tron_decomp
	int * _pairs_to_exclude, int * _pairs_to_exclude_a, int * _pairs_to_exclude_b,
	//output variables
	double* solution, double* _f, int* convergence 
){

    n1n2 = *_n1n2;
    int p=*_p;
	K = _K;
    X_diff=_X_diff;
    S2_deriv_ij=_S2_deriv_ij;
    lambda=*_lambda;
    zeta=*_zeta;
    b=*_b;
    int frtol=*_frtol;
    
	eta = (double *) malloc(sizeof(double)*n1n2);
	HESS_sauc_dca = (double *) malloc(sizeof(double)*p*p);
	
	double* x_lo=(double *) malloc(sizeof(double)*p);
	double* x_up=(double *) malloc(sizeof(double)*p);
	for (int i=0; i<p; i++) {
        x_lo[i] = -DBL_MAX;
        x_up[i] = DBL_MAX;
    }	

	// copied from solvebqp.c of kernlab  
	double fatol = 0;
	double cgtol = 0.1;	
	// youyi
	double fmin = -DBL_MAX;

    //the settings below are relaxed if not convergent no more than 'MAX_ATTEMPT' number of times    
    //assure MAX_ATTEMPT is small enough so that SQRT_MACHINE_DOUBLE_EPS * 10^(MAX_ATTEMPT) <= 0.1 
    int MAX_ATTEMPT = (int)fabs(log(*_gtol) / log(10.0)); 
    int cnt = 1;

    while(!*convergence){
                         
        if(*verbose)PRINTF("Tron #%i, fatol(%.10f) frtol(%.10f) cgtol(%.10f) gtol(%.10f) maxfev(%i)\n",
                        cnt,fatol,frtol,cgtol,*_gtol,*_maxfev);
        nfev_sauc_dca = 0;         
        // p is number of covariates, it is called n in dtron
        *_f=dtron(p, solution, x_lo, x_up, *_gtol, frtol, fatol, fmin, *_maxfev, cgtol,convergence, *verbose, &ugradhes_sauc_dca, &ufv_sauc_dca);
        
        if(cnt++ >= MAX_ATTEMPT) break;
        fatol *= 10.0;
        frtol *= 10.0;
        cgtol *= 10.0;
        *_gtol *= 10.0;
        *_maxfev *= 2;    
    } 
    
	
	free(eta);
    free(HESS_sauc_dca);
    free(x_lo);
    free(x_up);
}


