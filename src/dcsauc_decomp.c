// with decomposition

// c++
//#include <sstream>
//#include "U.h"

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

//extern "C" {
    #include "matrix.h" //crossprod() and tcrossprod() as R's source code
    double dtron(int n, double *x, double *xl, double *xu, double gtol, double frtol, double fatol, double fmin, int maxfev, double cgtol,int* convergence, int verbose,
       int (*ugradhesptr)(int, double *, double *, double **), int (*ufvptr)(int, double *, double *));
//}

#include "working_set.h" // order functions from R

// MACROS for saud_dca
#define	S1(x)	(log(1+exp(-b*x))/zeta)
#define	S1_deriv(x)   (-b/(1+exp(b*x))/zeta)
#define S1_dd(x)    (b*b/(1+exp(b*x))/(1+exp(-b*x))/zeta)


//ofstream logFile;

// data input
static int n1n2;
static double lambda;
static double * X_diff;
static double * S2_deriv_ij;
static double * K; //p x p, p = n for nonlinear kernel 
static double zeta;
static double b;
static int p; // length of full parameter vector
static double * theta; // full parameter vector
static int verbose;


static double * eta; // hold intermediate variables used in gradient and hessian updates
static double * temp_v_n1n2; // a vector of length n1n2 for a temporary vector in the computation of gradient
static double * temp_v_p; // a vector of length p for a temporary vector 
static double * HESS_sauc_dca; // pointer to a hessian qxq matrix
static double * full_gradient_sauc_dca; // pointer to a gradient vector of length p
static int * full_gradient_sauc_dca_order;
static double * working_set_x; // pointer to a vector of length q
static int * working_set;
static int q; // size of working set


static int nfev_sauc_dca; // counter for number of function evaluations


// compute the full gradient vector at theta 
int get_full_gradient_sauc_dca ()
{
    int ONE=1; 
    
    //  g := t(X.diff) %*% (S1.deriv(eta, zeta) - S2.deriv.ij) + 0
    tcrossprod(theta,&ONE,&p,X_diff,&n1n2,&p,eta);	
    for (int i=0; i<n1n2; i++) temp_v_n1n2[i] = S1_deriv(eta[i]) - S2_deriv_ij[i];
	crossprod(X_diff,&n1n2,&p,temp_v_n1n2,&n1n2,&ONE,full_gradient_sauc_dca);	

    //  g += lambda * beta or g += lambda * K %*% alpha; for linear kernel, K is pxp identity matrix, otherwise it is the p x p kernel matrix
	crossprod(theta, &p,&ONE,K, &p, &p, temp_v_p);
	for (int i=0; i<p; i++) full_gradient_sauc_dca[i] += lambda*temp_v_p[i];
	
	return 0;
}


// update working_set and x
void sauc_dca_working_set (double *x, int pairs_to_exclude, int * pairs_to_exclude_a, int * pairs_to_exclude_b) {
     
     //if(verbose) logFile << "         inside sauc_dca_working_set\n";
     get_full_gradient_sauc_dca ();
    // take absolute value before sorting
    for (int i=0; i<p; i++) full_gradient_sauc_dca[i] = fabs(full_gradient_sauc_dca[i]);
     // initialize full_gradient_sauc_dca_order to 1:p
     for (int i=0; i<p; i++) full_gradient_sauc_dca_order[i]=i+1;
     // order the gradients
     revsort(full_gradient_sauc_dca, full_gradient_sauc_dca_order, p);
     
//     if (verbose) {
//        logFile << "full_gradient_sauc_dca_order: ";
//        U::printArray(full_gradient_sauc_dca_order, p, logFile);
//        logFile << "full_gradient_sauc_dca sorted abs: ";
//        U::printArray(full_gradient_sauc_dca, p, logFile);
//     }                           
     
     for (int i=0; i<q; i++) working_set[i]=full_gradient_sauc_dca_order[i]-1;     
//    logFile << "initial working set (0 indexed): ";
//    U::printArray(working_set, q, logFile);
        
     int n_considered=q;
     // if some values of K[working_set, working_set] are too close to 1, tron will behave badly, leading to too big a x
     for (int i=0; i<pairs_to_exclude; i++) {
         if (n_considered>=p) break;

         int match_a=-1;
         for (int j=0; j<q; j++) 
             if (working_set[j]==pairs_to_exclude_a[i]) {
                match_a=j;
                break;
             }

         int match_b=-1;
         for (int j=0; j<q; j++) 
             if (working_set[j]==pairs_to_exclude_b[i]) {
                match_b=j;
                break;
             }
         
         if (match_a!=-1 && match_b!=-1) {
            working_set[match_a]=full_gradient_sauc_dca_order[n_considered++]-1;
//             if (verbose) {
//                logFile << "changing working set (0 indexed): ";
//                U::printArray(working_set, q, logFile);
//             }                           
         }
     }         
     
     for (int i=0; i<q; i++) x[i] = theta[working_set[i]];

}


// q is the length of x
// x is read-only parameter vector
// g is gradient at x
// H is hessian at x
int ugradhes_sauc_dca_decomp(int q, double *x, double *g, double **H)
{
//     if(verbose) logFile << "       inside ugradhes_sauc_dca_decomp\n";
    int ONE=1; 
    
    // update theta at the working set    
    for (int i=0; i<q; i++) theta[working_set[i]] = x[i];     
     
    // eta = theta %*% t(X_diff)
    tcrossprod(theta,&ONE,&p,  X_diff,&n1n2,&p, eta);	

    /* Hessian */
    
    *H = HESS_sauc_dca;
    // initialize HESS_sauc_dca to 0    
	memset(HESS_sauc_dca, 0, sizeof(double)*q*q);	
    double t;
    for (int i=0; i<n1n2; i++) {
        t=S1_dd(eta[i]);
        for (int j=0; j<q; j++) working_set_x[j] = (X_diff+working_set[j]*n1n2)[i];
        //  H: = S1_dd(eta[i], zeta) * trans(working_set_x) %*% working_set_x) + H
    	F77_CALL(dger)(&q, &q, &t, working_set_x, &ONE, working_set_x, &ONE, HESS_sauc_dca, &q);
     }
     
     for (int i=0; i<q; i++) 
         for (int j=0; j<q; j++) 
             HESS_sauc_dca[i*q+j] += lambda * K[working_set[i]*p + working_set[j]]; 
    
//        logFile << "HESS_sauc_dca: " << endl;
//        U::print2Dmatrix(HESS_sauc_dca, p, p, logFile);
    
    /* Gradient */
     
    //  g := t(X.diff) %*% (S1.deriv(eta, zeta) - S2.deriv.ij)
    for (int i=0; i<n1n2; i++) 
        temp_v_n1n2[i] = S1_deriv(eta[i]) - S2_deriv_ij[i];    	
    for (int i=0; i<q; i++) {
        g[i]=0;
        for (int j=0; j<n1n2; j++) 
            g[i] += (X_diff+working_set[i]*n1n2)[j] * temp_v_n1n2[j];
    }

    //  g += lambda * beta or g += lambda * K %*% alpha; for linear kernel, K is pxp identity matrix, otherwise it is the p x p kernel matrix
    for (int i=0; i<q; i++) {
        for (int j=0; j<p; j++) 
            g[i] += lambda * K[working_set[i]*p + j] * theta[j];
    }
	
//        logFile << "gradient: \n";
//        U::printArray(g, q, logFile);
	return 0;

}

int ufv_sauc_dca_decomp(int q, double *x, double *f)
{
//     if(verbose) logFile << "         inside ufv_sauc_dca_decomp\n";
    *f=0.0;
    int ONE=1;
    
    // update theta at the working set    
    for (int i=0; i<q; i++) theta[working_set[i]] = x[i];     
    
    // loss    
    // eta := X.diff %*% theta + 0
	tcrossprod(theta,&ONE,&p,X_diff,&n1n2,&p,eta);  // t(x) %*% X, X _{n1n2 x p} stored column-major 
    //sum(S1(eta, zeta) - eta * S2.deriv.ij)
    for (int i=0; i<n1n2; i++) 
        *f += S1(eta[i]) - eta[i] * S2_deriv_ij[i];
        
//    if (verbose) logFile << "loss " << *f; 
    
    // penalty    
	//*f += lambda * (||beta||_2)^2 or *f += lambda * t(alpha) %*% K %*% alpha
	crossprod(theta, &p, &ONE,K, &p, &p, temp_v_p); //eta = K %*% x	
	double xKx = 0.0;
	for (int i=0; i<p; i++) 
        xKx += temp_v_p[i] * theta[i]; // += t(x) %*% K %*% x        
    *f += .5 * lambda * xKx;
//    if (verbose) logFile << ", penalty " << .5 * lambda * xKx << ", f: " << *f << endl; 
	
    return ++nfev_sauc_dca;
}


//extern "C" {  

void dcsauc_tron_decomposition (
	double* _K, double* _X_diff, double* _S2_deriv_ij, 
	double* _lambda, int* _p, int* _n1n2, double * _zeta, double * _b,
	int* _maxfev, double* _gtol, double *_frtol, int* _verbose, int * _q,
	int * _pairs_to_exclude, int * _pairs_to_exclude_a, int * _pairs_to_exclude_b,
	//output variables
	double* _theta, double* _f, int* convergence 
){

    n1n2 = *_n1n2;
    p=*_p;
	K = _K;
    X_diff=_X_diff;
    S2_deriv_ij=_S2_deriv_ij;
    lambda=*_lambda;
    zeta=*_zeta;
    b=*_b;
	theta = _theta;
	verbose = *_verbose;
	q = *_q;
    
//    if (verbose) U::openWrite (logFile, "sauc_dca_tron_decomp.log", cout, true);
    
	eta = (double *) malloc(sizeof(double)*n1n2);	
	HESS_sauc_dca = (double *) malloc(sizeof(double)*q*q);

	double* x_lo=(double *) malloc(sizeof(double)*q);
	double* x_up=(double *) malloc(sizeof(double)*q);
	for (int i=0; i<q; i++) {
        x_lo[i] = -DBL_MAX;
        x_up[i] = DBL_MAX;
    }	

    working_set = (int *) malloc(sizeof(int)*q);    
	working_set_x = (double *) malloc(sizeof(double)*q);	
	temp_v_n1n2 = (double *) malloc(sizeof(double)*n1n2);	
	temp_v_p = (double *) malloc(sizeof(double)*p);	

	full_gradient_sauc_dca = (double *) malloc(sizeof(double)*p);	
	full_gradient_sauc_dca_order = (int *) malloc(sizeof(int)*p);
	
	double* x=(double *) malloc(sizeof(double)*q);		

	// copied from solvebqp.c of kernlab  
	double fatol = 0;
	double cgtol = 0.1;	
	// youyi
	double fmin = -DBL_MAX;
	
	double f, fc=DBL_MAX; // fc hold previous function value, f current


    //the settings below are relaxed if not convergent no more than 'MAX_ATTEMPT' number of times    
    //assure MAX_ATTEMPT is small enough so that SQRT_MACHINE_DOUBLE_EPS * 10^(MAX_ATTEMPT) <= 0.1 
    //int MAX_ATTEMPT = (int)fabs(log(*_gtol) / log(10.0)); 
    int MAX_ATTEMPT = 100;
    int cnt = 1;
    
   while(1){
        
        if (verbose>1) PRINTF("\n");                 
        if(verbose) PRINTF("tron #%i, ", cnt);
        //if(verbose) PRINTF("tron #%i, fatol(%.10f) frtol(%.10f) cgtol(%.10f) gtol(%.10f) maxfev(%i) ", cnt,fatol,frtol,cgtol,*_gtol,*_maxfev);
        
//        double ** H=(double **) malloc(sizeof(double*)*1);
//
//     for (int i=0; i<q; i++) {
//         working_set[i]=i;
//         x[i] = theta[working_set[i]];
//     }
//        ugradhes_sauc_dca_decomp(q, x, full_gradient_sauc_dca, H);
//            logFile << "gr at theta: ";
//            U::printArray(full_gradient_sauc_dca, q, logFile);
        
        sauc_dca_working_set(x, *_pairs_to_exclude, _pairs_to_exclude_a, _pairs_to_exclude_b);
        
//        ugradhes_sauc_dca_decomp(q, x, full_gradient_sauc_dca, H);
//            logFile << "gr at x: ";
//            U::printArray(full_gradient_sauc_dca, q, logFile);

        
//        if (verbose) logFile << "call tron" << endl;
        nfev_sauc_dca = 0; // this needs to be reset for each dtron call              
        f=dtron(q, x, x_lo, x_up, *_gtol, *_frtol, fatol, fmin, *_maxfev, cgtol,convergence, verbose, &ugradhes_sauc_dca_decomp, &ufv_sauc_dca_decomp);
        
//        if (verbose) {
//            logFile << "tron returns" << endl;
//            logFile << "theta: ";
//            U::printArray(theta, p, logFile);
//            logFile << endl;
//        }
    
        //if (cnt>1 && (fc-f)/fc<frtol) { 
        if (cnt>1 && (fc-f)/fc<1e-6) { 
           break;
        } else fc=f;
        
        if(cnt++ >= MAX_ATTEMPT) break;
        
//        fatol *= 10.0;
//        frtol *= 10.0;
//        cgtol *= 10.0;
//        *_gtol *= 10.0;
//        *_maxfev *= 2;    
    } 
    
    *_f=f;
    	
	free(eta);
    free(HESS_sauc_dca);
    free(x_lo);
    free(x_up);
    
	free(temp_v_n1n2);
	free(temp_v_p);
    free(working_set);
    free(working_set_x);
    
    free(full_gradient_sauc_dca);
    free(full_gradient_sauc_dca_order);
    
    free(x);
    
//    if (verbose) logFile.close();
}


//}
