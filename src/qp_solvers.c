#include "loqo.h"
#include "minQuad_QP.h"
#include <R.h>
#define PRINTF Rprintf
#include <limits.h>

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
#include <errno.h>
#define	MAX(A,B)	((A) > (B) ? (A) : (B))

// TRON DECLARATIONS
//  min { f(x)=0.5*x'*A*x + g0'*x : xl <= x <= xu }
// notation below: A = HESS, g0 = LIN
extern double dtron(int n, double *x, double *xl, double *xu, double gtol, double frtol, double fatol, double fmin, int maxfev, double cgtol,int* convergence, int verbose,
       int (*ugradhesptr)(int, double *, double *, double **), int (*ufvptr)(int, double *, double *)
);
static int nfev;
static int inc = 1;
static double one = 1.0;
static double zero = 0.0;
static double *HESS = NULL; //hessian
static double *LIN = NULL; //linear
int ugradhes(int n, double *x, double *g, double **H)
{
	*H = HESS;

	/* evaluate the gradient g = A*x + g0 */
	memcpy(g, LIN, sizeof(double)*n);
	F77_CALL(dsymv)("U", &n, &one, HESS, &n, x, &inc, &one, g, &inc FCONE);

	return 0;
}

int ufv(int n, double *x, double *f)
{
	/* evaluate the function value f(x) = 0.5*x'*A*x + g0'*x */  
	double *t = (double *) malloc(sizeof(double)*n);
	if(!t){
	#ifdef DEBUG_QP_SOLVER
		FILE* _F = fopen("qp_solver.debug","w+");
		if(!_F){
			int errnum=errno;
			PRINTF("uvf() in tron could not open debug file, reason: %s\n",strerror(errnum));
			exit(1);
		}
		fprintf(_F,"Error: failure to allocate memory in function uvf() in tron, exiting with code 1.\n");
		fclose(_F);
		exit(1);
	#endif
	return INT_MAX;
	}
	F77_CALL(dsymv)("U", &n, &one, HESS, &n, x, &inc, &zero, t, &inc FCONE);
	*f = 0.5*F77_CALL(ddot)(&n, x, &inc, t, &inc) + F77_CALL(ddot)(&n, x, &inc, LIN, &inc);
	free(t);
	return ++nfev;
}


/* driver for positive semidefinite quadratic programing version of tron */
// dtron: min { f(x)=0.5*x'*H*x + b'*x : xl <= x <= xu }
void optimize_qp_tron(int* n,double* hessian, double* linear, double* x_init,double* x_lo,double* x_up,double* solution,
double* control,int* verbose, int* convergence)
{

    
    // dtron's global variables must be (re)set every time dtron is called
    HESS  = hessian;
    LIN = linear; 
    nfev = 0;

	// control parameters stored in alphabetical order
	double cgtol = control[0];
	double fatol = control[1];
	double fmin = control[2];//-DBL_MAX;
	double frtol = control[3];
	double gtol = control[4];
	int    maxfev = (int)control[5]; // 1e3-5e3     

    //the settings below are relaxed if not convergent no more than 'MAX_ATTEMPT' number of times    
    //assure MAX_ATTEMPT is small enough so that SQRT_MACHINE_DOUBLE_EPS * 10^(MAX_ATTEMPT) <= 0.1 
    int MAX_ATTEMPT = (int)fabs(log(gtol) / log(10.0)); 
    int cnt = 1;

    for(int q1 = 0; q1 < *n; q1++)solution[q1] =  x_init[q1];    
    while(!*convergence){
        if(*verbose)PRINTF("Run no. %i of %i of TRON optimizer with fatol(%.10f) frtol(%.10f) cgtol(%.10f) gtol(%.10f) maxfev(%i)\n",
                        cnt,MAX_ATTEMPT,fatol,frtol,cgtol,gtol,maxfev);
        nfev = 0; 
        dtron(*n , solution, x_lo, x_up, gtol, frtol, fatol, fmin, maxfev, cgtol,convergence, *verbose, &ugradhes, &ufv);
        if(cnt++ >= MAX_ATTEMPT) break;
        fatol *= 10.0;
        frtol *= 10.0;
        cgtol *= 10.0;
        gtol *= 10.0;
        maxfev *= 2;    
    } 
}


int optimize_qp_loqo(MNQD_QP* qp, double* solution,double* control,int verbose)               
{
   
    double primal_obj,dual_obj;
	int iterations,convergence,error_code;
	int maxit = (int)control[3];
	double* dual = (double*)malloc(MAX(qp->n_constr,1) * sizeof(double)); 
	if(!dual) return 0;
		
	double* v = NULL;
    double*	r = NULL;
	if(qp->n_constr){	
	    // lhs <= A*x <= rhs -> LOQO parametrization is: v <= A*x <= v + r
		v = qp->lhs_constr;
		r = (double*)malloc(qp->n_constr * sizeof(double)); 
		for(int i = 0;i < qp->n_constr;i++) r[i] = qp->rhs_constr[i] - qp->lhs_constr[i]; 
	}

	// for(int i = 0;i < qp->n;i++) solution[i] = 0.5;
	// print_QP(qp);
	//Rprintf("solution before ");for(int i = 0;i < qp->n;i++)Rprintf(" %f",solution[i]);Rprintf("\n");
    if(verbose)
		Rprintf("loqo controls bound(%f) inf(%f) margin(%f) sigfig(%f) maxit(%i) verbose(%i)\n",
			control[0],control[1],control[2],control[4],(int)control[3],verbose);
	
	loqo(&qp->n,qp->b,qp->H,qp->xl,qp->xu,
		&qp->n_constr,qp->mat_constr,v,r,
		NULL,//_lapack_solver
		control+4,&maxit,control+2,control,//sigfig_max,maxiter,margin,bound
		control+1,&verbose,NULL,NULL,//inf,verb,buf_size,buf
		solution,dual,&primal_obj,&dual_obj,
	    &iterations,&convergence,&error_code);

   //Rprintf("solution after  ");for(int i = 0;i < qp->n;i++)Rprintf(" %f",solution[i]);Rprintf("\n");
   free(dual);
   if(qp->n_constr) free(r);
   if((LOQO_STATUS)convergence == CONVERGED) return 1;
   if(verbose)Rprintf("Nonconvergence in LOQO: %s \n",LOQO_STATUS_NAMES[convergence]);
   return 0;
}
