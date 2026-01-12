// Author: Krisztian Sebestyen 05/2013
// C implementation of the R function ipop() in package 'kernlab'

//kernlab:
//ipop solves the quadratic programming problem
//minimize   c' * primal + 1/2 primal' * H * primal
//subject to b <= A*primal <= b + r
//           l <= x <= u
//returns primal and dual variables (i.e. x and the Lagrange
//multipliers for b <= A * primal <= b + r
//for additional documentation see
//     R. Vanderbei
//     LOQO: an Interior Point Code for Quadratic Programming, 1992
// Author:      R version Alexandros Karatzoglou, orig. matlab Alex J. Smola
// Created:     12/12/97
// R Version:   12/08/03
// Updated:     13/10/05
// This code is released under the GNU Public License

#include "loqo.h"
#include "matrix.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>


#define	MAX(A,B)	((A) > (B) ? (A) : (B))
#define	MIN(A,B)	((A) < (B) ? (A) : (B))
#define	ABS(A)  	((A) > 0 ? (A) : (-(A)))
#define EQ(x,y,eps)  ((ABS((x)-(y))) <= (eps) ? 1 : 0) 


// #define CONVERGED 0
// #define PRIMAL_DUAL_INFEASIBLE 1
// #define PRIMAL_INFEASIBLE 2
// #define DUAL_INFEASIBLE 3
// #define PRIMAL_INFEASIBLE_SOLUTION 4
// #define PRIMAL_UNBOUNDED 5
// #define DUAL_UNBOUNDED 6
// #define SLOW_CONVERGENCE 7 //change bound ?
// #define MEMORY_ALLOCATION_FAILURE 8
// #define LAPACK_ERROR 9 
// #define LOGIC_ERROR 10
const char* LOQO_STATUS_NAMES[] = {
		"converged","primal and dual infeasible",
		"primal infeasible","dual infeasible", 
		"primal infeasible solution",
		"primal unbounded","dual unbounded","slow convergence, change bound ?",
		"memory allocation failure","LAPACK error"};

//input: 			 
// m - no. of constraints
// n - sample size (no. of data points)
//output:
// double and int buffer sizes
void R_get_buffer_size_loqo(int* solver,int* update,int* _m,int* _n,int* dsize,int* isize){
get_buffer_size_loqo(*solver,*update,*_m,*_n,dsize,isize);
}

void get_buffer_size_loqo(int solver,int update_solver,int m,int n,int* dsize,int* isize){

	int lapack_solver = (LAPACK_SOLVER)(solver);
	m = MAX(1,m);
	n = MAX(1,n);
	
	int m_plus_n = m + n;
	*isize = m_plus_n;
	*dsize = m_plus_n * m_plus_n + 30 * m_plus_n;

	if(update_solver && (lapack_solver != LU)){
		Rprintf("Warning: updating is only supported for LU-factorization.");
		update_solver = 0;
	}
	
	
	if(lapack_solver == LU){
		if(update_solver){
			*dsize += 2 * m_plus_n * m_plus_n; //L,U
		}
	}else if(lapack_solver == QR){
		if(update_solver){
			*dsize += 2 * (m_plus_n + 1) * m_plus_n; //Q,R
			*dsize += 1 * m_plus_n * m_plus_n; // Rinv		
		}		
	}else if(lapack_solver == LDL){
		;
	}
	
	if(update_solver){
		*dsize += m_plus_n;//diag_KKT0_old
	}

	// add for (m = 1,n) in case NULL was passed by the user for the parameters belowd
	*dsize += 3*n + 2; 
	// double* lo  = _lo;//_{n x 1}
	// double* up  = _up;//_{n x 1}
	// double* A   = _A;//constr_mat;  // _{m x n}
	// double* b   = _b;//constr_vec1; // _{m x 1}
	// double* r   = _r;//constr_vec2; // _{m x 1}	
}
	

// KKT = [-H A' A I_m] + 'diag', a length 'm+n' vector to be added to the diagonal
// row_major: offset = i*ncx + j
// col_major: offset = i + j*nrx
// H_x_{n x n} in column-major order
// A_{m x n} in column-major order
// d: length 'n', subtracts from diag(KKT)[0..(n-1)]
// e: length 'm', sets diag(KKT)[n..(n + m)]
// output: KKT_{(m+n) x (m+n)}
void get_reduced_KKT_matrix(int n,int m,double* H,double* A,double* d,double* e,double* KKT){
	
	int m_plus_n = m + n;
    int i,j;
	
// Code below allows for 'updating' if any of the components are NULL
	
	// 'H' symmetric, fill out lower-triangular part including diagonal
	//-H_{n x n} in 2nd quadrant with diagonal
	if(H){
		for(i = 0;i < n;i++)
			for(j = 0;j <= i ;j++)
				KKT[i + j * m_plus_n] = -H[i + j * n];
	}
	// 'A' is not symmetric, fill out rectangular block of KKT below the diagonal
	// A_{m x n} in 3rd quadrant	
	if(A){
		for(i = 0;i < m;i++)
			for(j = 0;j < n;j++)
				KKT[(i + n) + j * m_plus_n] = A[i + j * m];
	}
	
	//I_m in 4th quadrant, zero-out below diagonal
	for(i = n;i < m_plus_n;i++)
		for(j = n;j < m_plus_n;j++)
			KKT[i + j * m_plus_n] = 0.0;
	
	// diagonal
	if(d){
		for(i = 0;i < n;i++) 
			KKT[i + i * m_plus_n] -= d[i];
	}
	
	if(e){
		for(i = n;i < m_plus_n;i++) 
			KKT[i + i * m_plus_n] = e[i-n];
	}
	
// fill out upper triangular part of KKT			
	for(i = 0;i < m_plus_n;i++)
		for(j = 0;j < i;j++)
			KKT[j + i * m_plus_n] = KKT[i + j * m_plus_n];

	
/* CODE BELOW WORKS	to return the initial(!) KKT matrix
	// fill out lower-triangular part
	for(i = 0;i < m_plus_n;i++){
		for(j = 0;j < i;j++){
			if(i < n){
					KKT[i + j * m_plus_n] = -H[i + j * n];//-H_{n x n} in 2nd quadrant
			}else{
				if(j < n){
					KKT[i + j * m_plus_n] = A[(i - n) + j * m];// A_{m x n} in 3rd quadrant
				}else{
					KKT[i + j * m_plus_n] = 0.0;//off-diagonal elements in 4th quadrant
				}
			}
		}
	}
	
	
	// // diagonal

	// subtract 'd' outside of this function
	for(i = 0;i < n;i++)
		KKT[i + i * m_plus_n] = -H[i + i * n] - 1.0
	
	// set 'e' instead of '1' outside of this function
	for(i = n;i < m_plus_n;i++)
		KKT[i + i * m_plus_n] = 1.0;
			
	
// fill out upper triangular part			
	for(i = 0;i < m_plus_n;i++)
		for(j = 0;j < i;j++)
			KKT[j + i * m_plus_n] = KKT[i + j * m_plus_n];
*/
}
void get_KKT(int* n,int* m,double* H,double* A,double* d,double* e,double* KKT){
get_reduced_KKT_matrix(*n,*m,H,A,d,e,KKT);
}



/*
The reduced KKT system 'KKT0' is SYMMETRIC QUASIDEFINITE

Cholesky positive definite decomposition fails :
SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )

Cholesky semi-positive definite decomposition fails :
DPSTRF()

LU factorization DGESV and in the loop split it into 2 calls: 
SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

Similarly for QR and LDL' decompositions
*/


/*
H_x:
We keep the alias 'H_x' for 'H' in case we use pass a decomposed 'H' that has dim (n x n0) 
if H is not decomposed:
in ipop(kernlab) H.x (initially the hessian 'H') is used to set the upper-left block of 
the reduced KKT matrix 'KKT0'. Only the diagonal of H.x changes from one iteration to the next: 
diag(H.x) <- diag(H) + d where d = rep(1,n) initially. 
Hence we do not create an extra copy of the matrix 'H': 'H_x <- H'.
Instead we keep a separate copy of the diagonal 'diag_KKT0_x' and update it accordingly, and pass
around the original copy of the hessian 'H' if needed to reset KKT0.

H_y:
H_y is the identity matrix and it is used to set the lower-right block of augmented matrix 'KKT0'.
Only its diagonal changes from one iterations to next: 
diag(H.y) <- e where e = rep(1,n) initially. Hence we only keep a copy of the diagonal diag_KKT0_y.
*/
// INPUT:
// _n           : number of data points, the 'sample size', number of primal variables
// linear       : 'c' in the optimization problem, length n
// hessian      : 'H' in the optimization problem, length n * n
// l            : lower bounds on primal variables, length n, l <= x
// u            : upper bounds on primal variables, length n  x <= u
// _m           : number of constraints of the form b <= Ax <= b + r
// A            : (m x n) constraint matrix
// b            : length-m vector in constraint
// r            : length-m vector in constraint
// lapack_solver   : 0 for LU-factorization, 1 for QR, 2 for LDL (symmetric using Bunch-Kaufman diagonal pivoting) 
// update_solver: 0/1, an attempt to rank-1 update the LU-factorization, it is slow hence not done for QR either
// sigfig       : Precision (hould default to 7 significant figures)
// maxiter      : Maximum number of iterations
// margin       : how close we get to the constrains
// bound        : Clipping bound for the variables}
// inf          : Large numeric value that is practically infinity for the objective functions (primal,dual)
//                should be less than DBL_MAX if loqo is called from R with .C(DUP = FALSE,NAOK = FALSE)                   
// verb         : Display convergence information during runtime,0,1,2,3,4}
// dbuffer_size : size of buffer passed externally, ignored if not greater than or eq. to one returned by get_buffer_size_loqo()
// _dbuffer     : pointer to external buffer

// OUTPUT:
// primal       : the solution, a length-n vector, 
// dual         : length-m vector
// primal_obj   : length-1, value of primal objective function at solution 'primal' 
// dual_obj     : length-1, value of dual objective function 
// counter      : number of iterations taken at end of alg.
// convergence  : convergence code
// error_code   : see loqo.h 

// Note: output variables cannot be NULL, the length of 'dual' must be at least '1', even if m = 0, see below. 
// Note: you can pass NULL to all input including the bounds and constraints on primal variables
// Note: you can also pass m = 0 or m = NULL for the number of constraints. In this case m = 1 is assumed
// and loqo expects dual to have length one. In this case LOQO internally does the following
// m = 0 is set to be m = 1, A = NULL is set to be A = (1,..,1) of length 'n' and the constraints amount to:
// lo <= x <= up with, i lo,up=NULL then lo = (0,..,0) and up=(1,..,1)
// sum(lo) <= sum(x) <= sum(up)

void loqo(int* _n,double* linear, double* hessian, double *_lo, double *_up,
int* _m, double* _A, double* _b,  double* _r, 
int* _lapack_solver,
double* _sigfig_max, int* _maxiter, double* _margin, double* _bound, 
double* _inf,int* _verb,
int* _dbuffer_size,double* _dbuffer,
double* primal,double* dual,
double* _primal_obj,double* _dual_obj,
int* counter,int* convergence,int* error_code)
{


///////////////////////////////////////////////////////////////////////	
	int i,n,m,m_plus_n,ione,info;
	double sigfig,alfa,primal_infeasibility,dual_infeasibility,primal_obj,dual_obj,mu,tmp,b_plus_1,c_plus_1;
	double* c   = linear;      // _{n x 1}
	double* H_x = hessian;     // _{n x n}
///////////////////////////////////////////////////////////////////////	

///////////////////////////////////////////////////////////////////////
// default parameters
	*convergence = LOGIC_ERROR;
	if(!_n) return;
	n = *_n; // no. of 'points'
	if(n < 1) return;
	m = 1;if(_m) m = *_m;// no. of constraints
	m = MAX(1,m); //if m < 1 then set m = 1 and pass 'b','r' so that b <= A*x <= b+r approximates unconstrained estimation
	m_plus_n = m + n;

	
	double sigfig_max = 7.0;if(_sigfig_max) sigfig_max = *_sigfig_max;
	int maxiter = 40;if(_maxiter) maxiter = *_maxiter;
	double margin = .05;if(_margin) margin = *_margin;
	double bound = 10.0;if(_bound) bound = *_bound;
	double inf = 1e6;if(_inf) inf = *_inf;//1e6 is default in kernlab's ipop
	int verb = 0;if(_verb) verb = *_verb;
	LAPACK_SOLVER lapack_solver = LU;if(_lapack_solver) lapack_solver = (LAPACK_SOLVER)*_lapack_solver;
///////////////////////////////////////////////////////////////////////
	
///////////////////////////////////////////////////////////////////////	
	int* ibuffer = NULL;
	double* dbuffer = NULL;
	double* _p = NULL; //alias to dbuffer, incremented accordingly
		
	int ibuffer_size = 0;
	int dbuffer_size = 0;
	get_buffer_size_loqo((int)lapack_solver,0,m,n,&dbuffer_size,&ibuffer_size);
	// allocate buffer and set all elements to 0 as required below by some variables
	
	ibuffer = (int *)malloc((size_t)ibuffer_size * sizeof(int));
	if(!ibuffer){
		//Rprintf("Unable to allcoate %i bytes in function %s\n",ibuffer_size * sizeof(int),"loqo");
		*convergence = MEMORY_ALLOCATION_FAILURE;		
		return;
	}	int* ipiv = ibuffer;

	
	int internal_buffer = 1;
	if(_dbuffer && _dbuffer_size)
		internal_buffer = dbuffer_size > (*_dbuffer_size);
	if(internal_buffer){
		dbuffer = (double *)calloc((size_t)dbuffer_size , sizeof(double));
		if(!dbuffer){
			//Rprintf("Unable to allcoate %i bytes in function %s\n",dbuffer_size * sizeof(double),"loqo");
			*convergence = MEMORY_ALLOCATION_FAILURE;
			return;
		}
	}else{
		dbuffer = _dbuffer;
		memset(dbuffer,0,(size_t)dbuffer_size * sizeof(double));
	}		
	int lwork = 10 * m_plus_n;


	_p = dbuffer;	
///////////////////////////////////////////////////////////////////////
	

	
///////////////////////////////////////////////////////////////////////	
// default parameters continued
	double* lo  = _lo;if(!lo){lo = _p;_p+=n;memset(lo,0,n*sizeof(double));}//_{n x 1}
	
	double* up  = _up;if(!up){up = _p;_p+=n;for(i = 0;i < n;i++)up[i] = 1.0;}//_{n x 1}
	
	// if _A == NULL then m = 1 and A is set to be _{m x n} with the constraint b <= sum(primal) <= b + r
	double* A   = _A;if(!A){A = _p;_p+=n;for(i = 0;i < n;i++)A[i] = 1.0;}//constr_mat;  // _{m x n}
	
	// if _b == NULL then m = 1 and b is set to be _{1 x 1}, b = -0.5*inf
	double* b   = _b;
	if(!b){
		b = _p;_p++;
		//b[0] = -0.5*inf;
		b[0] = 0.0;for(i = 0;i < n;i++)b[0] += lo[i];
	}//constr_vec1; // _{m x 1}
	
	// if _r == NULL then m = 1 and r is set to be _{1 x 1}, r = inf
	double* r   = _r;
	if(!r){
		r = _p;_p++;
		//r[0] = inf;
		r[0] = 0.0;for(i = 0;i < n;i++);r[0] += up[i];
		r[0] -= b[0]; // so that b+r <= sum(up)
	}//constr_vec2; // _{m x 1}		
///////////////////////////////////////////////////////////////////////
	
///////////////////////////////////////////////////////////////////////
	// declaration for KKT system KKT0 = [-H_x A' A H_y] [x, y] = [c_x c_y]
	// (x,y), (c_x,c_y) are contiguos !
	double* KKT0 = _p;_p+=m_plus_n * m_plus_n;
	double* x = _p;_p+=n;
	double* y = _p;_p+=m;
	double* c_x = _p;_p+=n;
	double* c_y = _p;_p+=m;
	double* d = _p;_p+=n;
	double* e = _p;_p+=m;
    double* diag_H = _p;_p+=n;	
	double* diag_KKT0 = _p;_p+=m_plus_n;
	double* diag_KKT0_x = diag_KKT0;
	double* diag_KKT0_y = diag_KKT0 + n;
	double* H_dot_x = _p;_p+=n;	
///////////////////////////////////////////////////////////////////////	
	
///////////////////////////////////////////////////////////////////////
// declare the following in contiguous memory for singular value decomposition
	// primal_infeasibility
    double* alpha = _p;_p+=m;
    double* rho = _p;_p+=m;
    double* nu = _p;_p+=n;
    double* tau = _p;_p+=n;
	// dual_infeasibility
	double* beta  = _p;_p+=m;
	double* sigma = _p;_p+=n;
///////////////////////////////////////////////////////////////////////
	
///////////////////////////////////////////////////////////////////////
    double* hat_alpha = _p;_p+=m;
    double* hat_beta = _p;_p+=m;
    double* hat_nu = _p;_p+=n;
    double* hat_tau = _p;_p+=n;	
///////////////////////////////////////////////////////////////////////	

///////////////////////////////////////////////////////////////////////
	double* g = _p;_p+=n;
	double* s = _p;_p+=n;
	double* t = _p;_p+=n;
	double* z = _p;_p+=n;

	double* p = _p;_p+=m;
	double* q = _p;_p+=m;
	double* v = _p;_p+=m;
	double* w = _p;_p+=m;
///////////////////////////////////////////////////////////////////////
	
///////////////////////////////////////////////////////////////////////	
    // (delta_x,delta_y) are contiguos to form solution in while-loop
	double* delta_x = _p;_p+=n;
	double* delta_y = _p;_p+=m;

	double* delta_g = _p;_p+=n;
	double* delta_s = _p;_p+=n;
	double* delta_t = _p;_p+=n;
	double* delta_z = _p;_p+=n;

	double* delta_p = _p;_p+=m;
	double* delta_q = _p;_p+=m;
	double* delta_v = _p;_p+=m;
	double* delta_w = _p;_p+=m;
		
	double* gamma_s = _p;_p+=n;
	double* gamma_z = _p;_p+=n;
	double* gamma_q = _p;_p+=m;
	double* gamma_w = _p;_p+=m;
///////////////////////////////////////////////////////////////////////
	double* work = _p;_p+=lwork;
///////////////////////////////////////////////////////////////////////	

///////////////////////////////////////////////////////////////////////
	ione = 1;
	info = 0;
	sigfig = 0.0;
	alfa = 1.0;
	primal_infeasibility = inf;
	dual_infeasibility = inf;
	dual_obj = inf;
	primal_obj = inf;
	mu = 0.0;
	tmp = 0.0;
	b_plus_1 = 1.0;
	c_plus_1 = 1.0;	
///////////////////////////////////////////////////////////////////////



    // b_plus_1 <- max(svd(b)$d) + 1
	C_singval_dgesvd(&ione,&m,b,&b_plus_1,&info);	
	b_plus_1++;
	if(info){
		*convergence = LAPACK_ERROR;	
		*error_code = info;
		free(ibuffer);
		if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
		return;
	}		

    // c_plus_1 <- max(svd(c)$d) + 1
	C_singval_dgesvd(&ione,&n,c,&c_plus_1,&info);	
	c_plus_1++;
	if(info){
		*convergence = LAPACK_ERROR;	
		*error_code = info;
		free(ibuffer);
		if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
		return;
	}		
///////////////////////////////////////////////////////////////////////	
	

///////////////////////////////////////////////////////////////////////
	
    // solve the system [-H_x A' A H_y] [x, y] = [c_x c_y]
	// AP <- matrix(0,m+n,m+n)
	// xp <- 1:(m+n) <= n
	// AP[xp,xp] <- -H_x
	// AP[xp == FALSE,xp] <- A
	// AP[xp,xp == FALSE] <- t(A)
	// AP[xp == FALSE, xp== FALSE] <- H_y
	// s.tmp <- solve(AP,c(c_x,c_y))
	// x <- s.tmp[1:n]
	// y <- s.tmp[-(1:n)]

/////////////////////////// starting point ///////////////////////////	
	// NOTE:
	// Performing the line below forces a copy of hessian H passed to ipop()
	// Instead we will only update AP[] instead by setting d = (1,1,..,1)
	// as it is done in while-loop
    // diag(H_x) <- H_diag + 1 <=> H_diag + d
    // H_y <- diag(1,m) initially, then 'e'
    // H_diag <- diag(H)  	
    // c_x <- c
    // c_y <- b	
	memcpy(c_x,c,n * sizeof(double));
	memcpy(c_y,b,m * sizeof(double));	
	
// set initial (d,e) so that diagonal of AP is not zero ! 	

	for(i = 0;i < n;i++)diag_H[i] = H_x[i*n + i];
	for(i = 0;i < n;i++){
		d[i] = 1.0;
		if(ABS(diag_H[i] - -1.0) < .1) d[i] = 2.0;	
	}
	for(i = 0;i < m;i++) e[i] = 1.0;
	
	
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////	
    // solve the system [-H_x A' A H_y] [x, y] = [c_x c_y]
	// KKT0 will be destroyed by LAPACK
	
	get_reduced_KKT_matrix(n,m,H_x,A,d,e,KKT0);	
	memcpy(x,c_x,m_plus_n * sizeof(double)); //x = [c_x,c_y] initially and will contain soln.
	
	
	// BOTH cholesky factorizations fail, the matrix is symmetric indefinite (or numerically close to it) 	
	// positive definite cholesky
	// F77_CALL(dpotrf)("U",&m_plus_n,KKT0,&m_plus_n,&info);
	// if(info)Rprintf("error code %d from Lapack routine '%s' %s\n", info, "dpotrf","matrix is not symmetric positive definite.");	

	// positive semidefinite cholesky
	// 	get_reduced_KKT_matrix(n,m,H_x,A,d,e,KKT0);
	// int RANK;
	// double TOL = 1e-8;
	// F77_CALL(dpstrf)("U", &m_plus_n,KKT0,&m_plus_n, ipiv, &RANK, &TOL, work, &info );
	// if(info)Rprintf("error code %d from Lapack routine '%s' %s\n", info, "dpstrf","matrix is not symmetric positive semi-definite.");	

	if(verb > 3){
		Rprintf("Augmented KKT matrix\n");
		//print_matrix(m+n,m+n,KKT0);
		for(i = 0;i < m_plus_n * m_plus_n;i++) Rprintf("%.20f ",KKT0[i]);
		Rprintf("\n");
	}


	if(lapack_solver == LDL){
		// LAPACK LDL' or UDU' factor of symmetric indefinite matrix
		
		int LWORK = -1;
		F77_CALL(dsytrf)("U",&m_plus_n,KKT0,&m_plus_n,ipiv,work,&LWORK,&info FCONE); 
		if(info){
			*convergence = LAPACK_ERROR;	
			*error_code = info;
			free(ibuffer);
			if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
			Rprintf("error code %d from Lapack routine '%s'\n", info, "dsytrf");
			return;
		}	
		
		double* WORK = work;
		LWORK = (int)work[0];
		if(LWORK > lwork) WORK = (double*)malloc((size_t)(LWORK * sizeof(double)));
		F77_CALL(dsytrf)("U",&m_plus_n,KKT0,&m_plus_n,ipiv,WORK,&LWORK,&info FCONE); 
		if(LWORK > lwork) free(WORK);
		if(info){
			*convergence = LAPACK_ERROR;	
			*error_code = info;
			free(ibuffer);
			if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
			Rprintf("error code %d from Lapack routine '%s'\n", info, "dsytrf");
			return;
		}
		
		// // fixed workspace
		// F77_CALL(dsytrf)("U",&m_plus_n,KKT0,&m_plus_n,ipiv,work,&lwork,&info); 
		// if(info){
			// *convergence = LAPACK_ERROR;	
			// *error_code = info;
			// free(ibuffer);
			// if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
			// Rprintf("error code %d from Lapack routine '%s'\n", info, "dsytrf");
			// return;
		// }	
		
		F77_CALL(dsytrs)("U",&m_plus_n,&ione, KKT0, &m_plus_n, ipiv, x,&m_plus_n,&info FCONE);
		if(info){
			Rprintf("error code %d from Lapack routine '%s'\n", info, "dsytrs");
			*convergence = LAPACK_ERROR;	
			*error_code = info;		
			free(ibuffer);
			if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
			return;
		}	
	}else if(lapack_solver == LU){
		// 1-step LU-solver
		// F77_CALL(dgesv)(&m_plus_n,&ione, KKT0, &m_plus_n, ipiv, delta_x,&m_plus_n,&info);		
		
		F77_CALL(dgetrf)(&m_plus_n,&m_plus_n,KKT0,&m_plus_n,ipiv,&info); 
		if(info){
			*convergence = LAPACK_ERROR;	
			*error_code = info;
			free(ibuffer);
			if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
			Rprintf("error code %d from Lapack routine '%s'\n", info, "dgetrf");
			return;
		}			
		

		F77_CALL(dgetrs)("N", &m_plus_n, &ione, KKT0, &m_plus_n,ipiv, x, &m_plus_n, &info);
		if(info){
			*convergence = LAPACK_ERROR;	
			*error_code = info;
			free(ibuffer);
			if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
			Rprintf("error code %d from Lapack routine '%s'\n", info, "dgetrs");			
			return;
		}	
	}else{
		*convergence = LOGIC_ERROR;	
		*error_code = 0;
		free(ibuffer);
		if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
		Rprintf("Unknown optimizer selected \n");
		return;
	}
///////////////////////////////////////////////////////////////////////

	for(i = 0;i < n;i++){
		if(EQ(lo[i],x[i],DBL_EPSILON)) 
			x[i] = lo[i];
		else if(EQ(up[i],x[i],DBL_EPSILON)) 
			x[i] = up[i];
	}
	
	//Rprintf("LAPACK initial solution=");for(i = 0;i < m_plus_n;i++) Rprintf("%.6e ",x[i]);Rprintf("\n");	
    
///////////////////////////////////////////////////////////////////////	
    // g <- pmax(abs(x - l), bound)
    // z <- pmax(abs(x), bound)
    // t <- pmax(abs(u - x), bound)
    // s <- pmax(abs(x), bound)
	for(i = 0;i < n;i++) g[i] = MAX(ABS(x[i] - 1.0),bound);
	for(i = 0;i < n;i++) z[i] = MAX(ABS(x[i]),bound);
	memcpy(s,z,n * sizeof(double));
	for(i = 0;i < n;i++) t[i] = MAX(ABS(up[i] - x[i]),bound);
	
    // v <- pmax(abs(y), bound)
    // w <- pmax(abs(y), bound)
    // p <- pmax(abs(r - w), bound)
    // q <- pmax(abs(y), bound)
	for(i = 0;i < m;i++) v[i] = MAX(ABS(y[i]) , bound);
	memcpy(w,v,m * sizeof(double));
	memcpy(q,v,m * sizeof(double));
	for(i = 0;i < m;i++) p[i] = MAX(ABS(r[i] - w[i]) , bound);


	// mu <- as.vector(crossprod(z,g) + crossprod(v,w) + crossprod(s,t) + crossprod(p,q))/(2 * (m + n))
	mu = 0.0;
	// for(i = 0;i < n;i++) mu += z[i]*g[i] + s[i]*t[i];
	// for(i = 0;i < m;i++) mu += v[i]*w[i] + p[i]*q[i];
	crossprod(z,&n,&ione,g,&n,&ione,&tmp);mu += tmp;
	crossprod(s,&n,&ione,t,&n,&ione,&tmp);mu += tmp;
	crossprod(v,&m,&ione,w,&m,&ione,&tmp);mu += tmp;
	crossprod(p,&m,&ione,q,&m,&ione,&tmp);mu += tmp;	
	mu /= (2.0 * (double)m_plus_n);
	
///////////////////////////////////////////////////////////////////////


	if((2 <= verb) && (verb <= 3)){
		Rprintf("Initial values:\n");
		Rprintf("mu=%.6e b+1=%.6e c+1=%.6e\n",mu,b_plus_1,c_plus_1);
		
		Rprintf("x=");for(i = 0;i < n;i++) Rprintf("%.6e ",x[i]);Rprintf("\n");		
		Rprintf("g=");for(i = 0;i < n;i++) Rprintf("%.6e ",g[i]);Rprintf("\n");
		Rprintf("s=");for(i = 0;i < n;i++) Rprintf("%.6e ",s[i]);Rprintf("\n");
		Rprintf("t=");for(i = 0;i < n;i++) Rprintf("%.6e ",t[i]);Rprintf("\n");
		Rprintf("z=");for(i = 0;i < n;i++) Rprintf("%.6e ",z[i]);Rprintf("\n");

		Rprintf("y=");for(i = 0;i < m;i++) Rprintf("%.6e ",y[i]);Rprintf("\n");		
		Rprintf("p=");for(i = 0;i < m;i++) Rprintf("%.6e ",p[i]);Rprintf("\n");
		Rprintf("q=");for(i = 0;i < m;i++) Rprintf("%.6e ",q[i]);Rprintf("\n");
		Rprintf("v=");for(i = 0;i < m;i++) Rprintf("%.6e ",v[i]);Rprintf("\n");
		Rprintf("w=");for(i = 0;i < m;i++) Rprintf("%.6e ",w[i]);Rprintf("\n");		
	}


	if (verb > 1)	                       // print at least one status report
	  Rprintf("Iter    PrimalInf  DualInf  SigFigs  Rescale  PrimalObj  DualObj\n");

	*counter = 0;	
    while ((*counter) < maxiter)
      {
        // update the iteration counter
        (*counter)++;	
	
        // central path (predictor)

		// H_dot_x <- H %*% x
		matprod(H_x,&n,&n,x,&n,&ione,H_dot_x);
		
        // rho <- b - A %*% x + w // m x 1
		matprod(A,&m,&n,x,&n,&ione,rho);
		for(i = 0;i < m;i++) rho[i] = b[i] - rho[i] + w[i];
		
        // nu <- l - x + g        // n x 1
		for(i = 0;i < n;i++) nu[i] = lo[i] - x[i] + g[i];

        // tau <- u - x - t       // n x 1
		for(i = 0;i < n;i++) tau[i] = up[i] - x[i] - t[i];

        // alpha <- r - w - p     // m x 1 
		for(i = 0;i < m;i++) alpha[i] = r[i] - w[i] - p[i];

		// CHECK LOGIC HERE: A FROM 'R' IS ACTUALLY tA
        // sigma <- c  - crossprod(A, y) - z + s + H_dot_x //n x 1
		crossprod(A,&m,&n,y,&m,&ione,sigma);
		for(i = 0;i < n;i++) sigma[i] = c[i] - sigma[i] - z[i] + s[i] + H_dot_x[i];
        
        // beta <- y + q - v
		for(i = 0;i < m;i++) beta[i] = y[i] + q[i] - v[i];

        // gamma_z <- - z
        // gamma_w <- - w
        // gamma_s <- - s
        // gamma_q <- - q
		for(i = 0;i < n;i++) gamma_z[i] = -z[i];
		for(i = 0;i < n;i++) gamma_s[i] = -s[i];
		for(i = 0;i < m;i++) gamma_w[i] = -w[i];
		for(i = 0;i < m;i++) gamma_q[i] = -q[i];
		
		if(verb > 3){
			Rprintf("Intermediate values:\n");
			Rprintf("mu=%.6e b+1=%.6e c+1=%.6e\n",mu,b_plus_1,c_plus_1);
			
			Rprintf("x=");for(i = 0;i < n;i++) Rprintf("%.6e ",x[i]);Rprintf("\n");		
			Rprintf("g=");for(i = 0;i < n;i++) Rprintf("%.6e ",g[i]);Rprintf("\n");
			Rprintf("s=");for(i = 0;i < n;i++) Rprintf("%.6e ",s[i]);Rprintf("\n");
			Rprintf("t=");for(i = 0;i < n;i++) Rprintf("%.6e ",t[i]);Rprintf("\n");
			Rprintf("z=");for(i = 0;i < n;i++) Rprintf("%.6e ",z[i]);Rprintf("\n");

			Rprintf("y=");for(i = 0;i < m;i++) Rprintf("%.6e ",y[i]);Rprintf("\n");
			Rprintf("p=");for(i = 0;i < m;i++) Rprintf("%.6e ",p[i]);Rprintf("\n");
			Rprintf("q=");for(i = 0;i < m;i++) Rprintf("%.6e ",q[i]);Rprintf("\n");
			Rprintf("v=");for(i = 0;i < m;i++) Rprintf("%.6e ",v[i]);Rprintf("\n");
			Rprintf("w=");for(i = 0;i < m;i++) Rprintf("%.6e ",w[i]);Rprintf("\n");		
		}		
	
        // instrumentation
        // x_dot_H_dot_x <-  crossprod(x, H_dot_x)
		double x_dot_H_dot_x;
		crossprod(x,&n,&ione,H_dot_x,&n,&ione,&x_dot_H_dot_x);
		
        // primal_infeasibility <- max(svd(rbind(rho, tau, matrix(alpha), nu))$d)/ b_plus_1
		int len = 2 * m_plus_n;
		double* alpha_rho_nu_tau = alpha;
		C_singval_dgesvd(&ione,&len,alpha_rho_nu_tau,&primal_infeasibility,&info);	
        primal_infeasibility /= b_plus_1;
		if(info){
			*convergence = LAPACK_ERROR;	
			*error_code = info;		
			free(ibuffer);
			if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
			return;
		}			
		// dual_infeasibility <- max(svd(rbind(sigma,t(t(beta))))$d) / c_plus_1
		C_singval_dgesvd(&ione,&m_plus_n,beta,&dual_infeasibility,&info);	
        dual_infeasibility /= c_plus_1;
		if(info){
			*convergence = LAPACK_ERROR;	
			*error_code = info;
			free(ibuffer);
			if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
			return;
		}			
        // primal_obj <- crossprod(c,x) + 0.5 * x_dot_H_dot_x
		primal_obj = 0.0;
		crossprod(c,&n,&ione,x,&n,&ione,&primal_obj);
		primal_obj += 0.5 * x_dot_H_dot_x;
		
		
        // dual_obj <- crossprod(b,y) - 0.5 * x_dot_H_dot_x + crossprod(l, z) - crossprod(u,s) - crossprod(r,q)
		dual_obj = -0.5 * x_dot_H_dot_x;
		crossprod(b,&m,&ione,y,&m,&ione,&tmp);dual_obj += tmp;
		crossprod(lo,&n,&ione,z,&n,&ione,&tmp);dual_obj += tmp;
		crossprod(up,&n,&ione,s,&n,&ione,&tmp);dual_obj -= tmp;
		crossprod(r,&m,&ione,q,&m,&ione,&tmp);dual_obj -= tmp;		
		// for(i = 0;i < n;i++) dual_obj += z[i]*l[i] - s[i]*u[i];
		// for(i = 0;i < m;i++) dual_obj += b[i]*y[i] - r[i]*q[i];
        
        // sigfig <- max(-log10(abs(primal_obj - dual_obj)/(abs(primal_obj) + 1)), 0)
        sigfig = MAX(-log10(ABS(primal_obj - dual_obj)/(ABS(primal_obj) + 1.0)), 0.0);
		
        // if (sigfig >= sigf) break
        // if (verb > 0)		      	# final report
          // cat( counter, "\t", signif(primal_infeasibility,6), signif(dual_infeasibility,6), sigfig, alfa, primal_obj, dual_obj,"\n")
		if ((1 < verb ) && (verb != 3)) Rprintf("%i    %.6e %.6e %.6e %.6e %.6e %.6e\n",(*counter), primal_infeasibility, dual_infeasibility, sigfig, alfa, primal_obj, dual_obj);
		if(verb == 3){ Rprintf("x=");for(i = 0;i < m_plus_n;i++) Rprintf("%.6f ",x[i]);Rprintf("\n");}		


		if(sigfig >= sigfig_max){
            *convergence = CONVERGED;
			break;
		}
///////////////////////////////////////////////////////////////////////		 
		// some more intermediate variables (the hat section)
        // hat_beta <- beta - v * gamma_w / w
        // hat_alpha <- alpha - p * gamma_q / q
        // hat_nu <- nu + g * gamma_z / z
        // hat_tau <- tau - t * gamma_s / s
		for(i = 0;i < m;i++) hat_alpha[i] = alpha[i] - p[i] * gamma_q[i] / q[i];				
		for(i = 0;i < m;i++) hat_beta[i] = beta[i] - v[i] * gamma_w[i] / w[i];
		for(i = 0;i < n;i++) hat_nu[i] = nu[i] + g[i] * gamma_z[i] / z[i];
		for(i = 0;i < n;i++) hat_tau[i] = tau[i] - t[i] * gamma_s[i] / s[i];
		
        // the diagonal terms
        // d <- z / g + s / t
        // e <- 1 / (v / w + q / p)
		for(i = 0;i < n;i++) d[i] = z[i]/g[i] + s[i] / t[i];
		for(i = 0;i < m;i++) e[i] = 1.0/(v[i]/w[i] + q[i] / p[i]);		

				
///////////////////////////////////////////////////////////////////////

		if(verb > 3){
			
			Rprintf("H***x=");for(i = 0;i < n;i++) Rprintf("%.6e ",H_dot_x[i]);Rprintf("\n");
			Rprintf("sigma=");for(i = 0;i < n;i++) Rprintf("%.6e ",sigma[i]);Rprintf("\n");
			Rprintf("d=");for(i = 0;i < n;i++) Rprintf("%.6e ",d[i]);Rprintf("\n");
			Rprintf("nu=");for(i = 0;i < n;i++) Rprintf("%.6e ",nu[i]);Rprintf("\n");
			Rprintf("tau=");for(i = 0;i < n;i++) Rprintf("%.6e ",tau[i]);Rprintf("\n");
			Rprintf("hat_nu=");for(i = 0;i < n;i++) Rprintf("%.6e ",hat_nu[i]);Rprintf("\n");
			Rprintf("hat_tau=");for(i = 0;i < n;i++) Rprintf("%.6e ",hat_tau[i]);Rprintf("\n");
			Rprintf("gamma_s=");for(i = 0;i < n;i++) Rprintf("%.6e ",gamma_s[i]);Rprintf("\n");
			Rprintf("gamma_z=");for(i = 0;i < n;i++) Rprintf("%.6e ",gamma_z[i]);Rprintf("\n");

			Rprintf("e=");for(i = 0;i < m;i++) Rprintf("%.6e ",e[i]);Rprintf("\n");
			Rprintf("gamma_w=");for(i = 0;i < m;i++) Rprintf("%.6e ",gamma_w[i]);Rprintf("\n");
			Rprintf("gamma_q=");for(i = 0;i < m;i++) Rprintf("%.6e ",gamma_q[i]);Rprintf("\n");
			Rprintf("rho=");for(i = 0;i < m;i++) Rprintf("%.6e ",rho[i]);Rprintf("\n");
			Rprintf("alpha=");for(i = 0;i < m;i++) Rprintf("%.6e ",alpha[i]);Rprintf("\n");
			Rprintf("beta=");for(i = 0;i < m;i++) Rprintf("%.6e ",beta[i]);Rprintf("\n");
			Rprintf("hat_alpha=");for(i = 0;i < m;i++) Rprintf("%.6e ",hat_alpha[i]);Rprintf("\n");
			Rprintf("hat_beta=");for(i = 0;i < m;i++) Rprintf("%.6e ",hat_beta[i]);Rprintf("\n");
		}	


///////////////////////////////////////////////////////////////////////
        // initialization before the big cholesky

        // diag(H_x) <- H_diag + d
        // diag(H_y) <- e
		for(i = 0;i < n;i++)diag_KKT0_x[i] = -(diag_H[i] + d[i]); // != 0.0
		memcpy(diag_KKT0_y,e,m * sizeof(double));	// != 0.0

        // c_x <- sigma - z * hat_nu / g - s * hat_tau / t
        // c_y <- rho - e * (hat_beta - q * hat_alpha / p)
		for(i = 0;i < m;i++)c_y[i] = rho[i] - e[i]*(hat_beta[i]-q[i]*hat_alpha[i]/p[i]);
		for(i = 0;i < n;i++)c_x[i] = sigma[i] - z[i]*hat_nu[i]/g[i]-s[i]*hat_tau[i]/t[i];
///////////////////////////////////////////////////////////////////////

// and solve the system [-H_x A' A H_y] [delta_x, delta_y] <- [c_x c_y]
		
///////////////////////////////////////////////////////////////////////		
		// AP[xp,xp] <- -H_x
		// AP[xp == FALSE, xp== FALSE] <- H_y
		// s1.tmp <- solve(AP,c(c_x,c_y))
		// delta_x <-s1.tmp[1:n]
		// delta_y <- s1.tmp[-(1:n)]
///////////////////////////////////////////////////////////////////////
		
///////////////////////////////////////////////////////////////////////	
		// and solve the system [-H_x A' A H_y] [delta_x, delta_y] = [c_x c_y]
		// split the 2 calls in dgesv in LAPACK to two separate calls dgesv = (dgetrf,dgetrs) to keep and reuse the factorization - twice as fast
				
		memcpy(delta_x,c_x,m_plus_n * sizeof(double)); //delta_x = [c_x,c_y] initially and will contain soln.	

		if(lapack_solver == LDL){
			// LDL' factorization for symmetric indefinite matrix
			get_reduced_KKT_matrix(n,m,H_x,A,d,e,KKT0);// KKT0 will be destroyed by LAPACK		
			
			int LWORK = -1;
			F77_CALL(dsytrf)("U",&m_plus_n,KKT0,&m_plus_n,ipiv,work,&LWORK,&info FCONE); 
			if(info){
				*convergence = LAPACK_ERROR;	
				*error_code = info;
				free(ibuffer);
				if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
				Rprintf("error code %d from Lapack routine '%s'\n", info, "dsytrf");
				return;
			}	
			double* WORK = work;
			LWORK = (int)work[0];
			if(LWORK > lwork) WORK = (double*)malloc((size_t)(LWORK * sizeof(double)));
			F77_CALL(dsytrf)("U",&m_plus_n,KKT0,&m_plus_n,ipiv,WORK,&LWORK,&info FCONE); 
			if(LWORK > lwork) free(WORK);
			if(info){
				*convergence = LAPACK_ERROR;	
				*error_code = info;
				free(ibuffer);
				if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
				Rprintf("error code %d from Lapack routine '%s'\n", info, "dsytrf");
				return;
			}			

			// // fixed workspace
			// F77_CALL(dsytrf)("U",&m_plus_n,KKT0,&m_plus_n,ipiv,work,&lwork,&info); 
			// if(info){
				// *convergence = LAPACK_ERROR;	
				// *error_code = info;
				// free(ibuffer);
				// if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
				// Rprintf("error code %d from Lapack routine '%s'\n", info, "dsytrf");
				// return;
			// }	

			F77_CALL(dsytrs)("U", &m_plus_n, &ione, KKT0, &m_plus_n,ipiv, delta_x, &m_plus_n, &info FCONE);
			if(info){
				*convergence = LAPACK_ERROR;	
				*error_code = info;
				free(ibuffer);
				if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
				Rprintf("error code %d from Lapack routine '%s'\n", info, "dsytrs");			
				return;
			}
		}else if(lapack_solver == LU){
		
				
			// LU-solve in 2 steps
			// store LU-factorization of KKT0 in 'KKT0' and pivot-permutation in 'ipiv'			
			get_reduced_KKT_matrix(n,m,H_x,A,d,e,KKT0);// KKT0 will be destroyed by LAPACK			
			F77_CALL(dgetrf)(&m_plus_n,&m_plus_n,KKT0,&m_plus_n,ipiv,&info); 
			if(info){
				*convergence = LAPACK_ERROR;	
				*error_code = info;
				free(ibuffer);
				if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
				Rprintf("error code %d from Lapack routine '%s'\n", info, "dgetrf");
				return;
			}					
			F77_CALL(dgetrs)("N", &m_plus_n, &ione, KKT0, &m_plus_n,ipiv, delta_x, &m_plus_n, &info);
			if(info){
				*convergence = LAPACK_ERROR;	
				*error_code = info;
				free(ibuffer);
				if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
				Rprintf("error code %d from Lapack routine '%s'\n", info, "dgetrs");			
				return;
			}
			
		}else{
			*convergence = LOGIC_ERROR;	
			*error_code = 0;
			free(ibuffer);
			if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
			Rprintf("Unknown optimizer selected \n");
			return;
		}
	  
		// (delta_x,delta_y) has solution      
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////        
        // backsubstitution
        // delta_w <- - e * (hat_beta - q * hat_alpha / p + delta_y)
        // delta_s <- s * (delta_x - hat_tau) / t
        // delta_z <- z * (hat_nu - delta_x) / g
        // delta_q <- q * (delta_w - hat_alpha) / p
        // delta_v <- v * (gamma_w - delta_w) / w
        // delta_p <- p * (gamma_q - delta_q) / q
        // delta_g <- g * (gamma_z - delta_z) / z
        // delta_t <- t * (gamma_s - delta_s) / s
		for(i = 0;i < m;i++) delta_w[i] = -e[i]*(hat_beta[i]-q[i]*hat_alpha[i]/p[i] + delta_y[i]);
		for(i = 0;i < n;i++) delta_s[i] = s[i]*(delta_x[i] - hat_tau[i]) / t[i];
		for(i = 0;i < n;i++) delta_z[i] = z[i]*(hat_nu[i] - delta_x[i]) / g[i];
		for(i = 0;i < m;i++) delta_q[i] = q[i]*(delta_w[i] - hat_alpha[i]) / p[i];

		for(i = 0;i < m;i++) delta_v[i] = v[i]*(gamma_w[i] - delta_w[i]) / w[i];		
		for(i = 0;i < m;i++) delta_p[i] = p[i]*(gamma_q[i] - delta_q[i]) / q[i];
		for(i = 0;i < n;i++) delta_g[i] = g[i]*(gamma_z[i] - delta_z[i]) / z[i];
		for(i = 0;i < n;i++) delta_t[i] = t[i]*(gamma_s[i] - delta_s[i]) / s[i];
///////////////////////////////////////////////////////////////////////
		if(verb > 3){

			Rprintf("delta_x=");for(i = 0;i < n;i++) Rprintf("%.6e ",delta_x[i]);Rprintf("\n");		
			Rprintf("delta_g=");for(i = 0;i < n;i++) Rprintf("%.6e ",delta_g[i]);Rprintf("\n");
			Rprintf("delta_s=");for(i = 0;i < n;i++) Rprintf("%.6e ",delta_s[i]);Rprintf("\n");
			Rprintf("delta_t=");for(i = 0;i < n;i++) Rprintf("%.6e ",delta_t[i]);Rprintf("\n");
			Rprintf("delta_z=");for(i = 0;i < n;i++) Rprintf("%.6e ",delta_z[i]);Rprintf("\n");

			Rprintf("delta_y=");for(i = 0;i < m;i++) Rprintf("%.6e ",delta_y[i]);Rprintf("\n");		
			Rprintf("delta_p=");for(i = 0;i < m;i++) Rprintf("%.6e ",delta_p[i]);Rprintf("\n");
			Rprintf("delta_q=");for(i = 0;i < m;i++) Rprintf("%.6e ",delta_q[i]);Rprintf("\n");
			Rprintf("delta_v=");for(i = 0;i < m;i++) Rprintf("%.6e ",delta_v[i]);Rprintf("\n");
			Rprintf("delta_w=");for(i = 0;i < m;i++) Rprintf("%.6e ",delta_w[i]);Rprintf("\n");		
		}


///////////////////////////////////////////////////////////////////////		
        // compute update step now (sebastian's trick)
        // alfa <- - (1 - margin) / min(c(delta_g / g, delta_w / w, delta_t / t, delta_p / p, delta_z / z, delta_v / v, delta_s / s, delta_q / q, -1))
		alfa = -1.0;
		for(i = 0;i < n;i++){
			alfa = MIN(alfa,delta_g[i]/g[i]);
			alfa = MIN(alfa,delta_s[i]/s[i]);
			alfa = MIN(alfa,delta_t[i]/t[i]);
			alfa = MIN(alfa,delta_z[i]/z[i]);
	    }
		for(i = 0;i < m;i++){
			alfa = MIN(alfa,delta_p[i]/p[i]);
			alfa = MIN(alfa,delta_q[i]/q[i]);
			alfa = MIN(alfa,delta_v[i]/v[i]);
			alfa = MIN(alfa,delta_w[i]/w[i]);
	    }	
		alfa = -(1.0 - margin) / alfa;
		
		// NOTE: check R-logic here (newmu)
        // newmu <- (crossprod(z,g) + crossprod(v,w) + crossprod(s,t) + crossprod(p,q))/(2 * (m + n))
        // newmu <- mu * ((alfa - 1) / (alfa + 10))^2
        // gamma_z <- mu / g - z - delta_z * delta_g / g
        // gamma_w <- mu / v - w - delta_w * delta_v / v
        // gamma_s <- mu / t - s - delta_s * delta_t / t
        // gamma_q <- mu / p - q - delta_q * delta_p / p		
		double newmu = ((alfa - 1.0) / (alfa + 10.0));
		newmu *= newmu;
		newmu *= mu;
	
		for(i = 0;i < n;i++) gamma_s[i] = mu / t[i] - s[i] - delta_s[i] * delta_t[i] / t[i];
		for(i = 0;i < n;i++) gamma_z[i] = mu / g[i] - z[i] - delta_z[i] * delta_g[i] / g[i];
		for(i = 0;i < m;i++) gamma_q[i] = mu / p[i] - q[i] - delta_q[i] * delta_p[i] / p[i];
		for(i = 0;i < m;i++) gamma_w[i] = mu / v[i] - w[i] - delta_w[i] * delta_v[i] / v[i];	
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////		
        // some more intermediate variables (the hat section)
        // hat_beta <- beta - v * gamma_w / w
        // hat_alpha <- alpha - p * gamma_q / q
        // hat_nu <- nu + g * gamma_z / z
        // hat_tau <- tau - t * gamma_s / s
		for(i = 0;i < m;i++) hat_alpha[i] = alpha[i] - p[i] * gamma_q[i] / q[i];
		for(i = 0;i < m;i++) hat_beta[i] = beta[i] - v[i] * gamma_w[i] / w[i];
		for(i = 0;i < n;i++) hat_nu[i] = nu[i] + g[i] * gamma_z[i] / z[i];
		for(i = 0;i < n;i++) hat_tau[i] = tau[i] - t[i] * gamma_s[i] / s[i];
///////////////////////////////////////////////////////////////////////		
		
		
///////////////////////////////////////////////////////////////////////		
		// NOTE: check R-logic here
		// diag(H_x) and H_y are not updated here
		
        // initialization before the big cholesky
        // //for (  i  in  1 : n H_x(i,i) <- H_diag(i) + d(i) ) {
        // //H_y <- diag(e)
        // c_x <- sigma - z * hat_nu / g - s * hat_tau / t
        // c_y <- rho - e * (hat_beta - q * hat_alpha / p)

		for(i = 0;i < m;i++) c_y[i] = rho[i] - e[i]*(hat_beta[i]-q[i]*hat_alpha[i]/p[i]);
		for(i = 0;i < n;i++) c_x[i] = sigma[i] - z[i]*hat_nu[i]/g[i]-s[i]*hat_tau[i]/t[i];
		
		// NOTE: check R-logic here
		// R-code unnecessarily resets KKT0 but C-code must

        // and solve the system [-H_x A' A H_y] [delta_x, delta_y] <- [c_x c_y]
		// AP[xp,xp] <- -H_x
		// AP[xp == FALSE, xp== FALSE] <- H_y
		// s1.tmp <- solve(AP,c(c_x,c_y))
		// delta_x<-s1.tmp[1:n] ; delta_y<-s1.tmp[-(1:n)]
///////////////////////////////////////////////////////////////////////	
				
		memcpy(delta_x,c_x,m_plus_n * sizeof(double)); //delta_x = [c_x,c_y] initially and will contain soln.	

		if(lapack_solver == LDL){
			// LDL-solve step reusing LDL-factorization in KKT0		
			F77_CALL(dsytrs)("U", &m_plus_n, &ione, KKT0, &m_plus_n,ipiv, delta_x, &m_plus_n, &info FCONE);
			if(info){
				*convergence = LAPACK_ERROR;	
				*error_code = info;
				free(ibuffer);
				if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
				Rprintf("error code %d from Lapack routine '%s'\n", info, "dsytrs");			
				return;
			}
		}else if(lapack_solver == LU){
			// LU-solve step reusing LU-factorization in KKT0
			F77_CALL(dgetrs)("N", &m_plus_n, &ione, KKT0, &m_plus_n,ipiv, delta_x, &m_plus_n, &info);
			if(info){
				*convergence = LAPACK_ERROR;	
				*error_code = info;
				free(ibuffer);
				if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
				Rprintf("error code %d from Lapack routine '%s'\n", info, "dgetrs");
				return;
			}
		}else{
			*convergence = LOGIC_ERROR;	
			*error_code = 0;
			free(ibuffer);
			if(internal_buffer)free(dbuffer);else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
			Rprintf("Unknown optimizer selected \n");
			return;
		}

		// (delta_x,delta_y) has soln. 		
		
		
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
     	  // backsubstitution
        // delta_w <- - e * (hat_beta - q * hat_alpha / p + delta_y)
        // delta_s <- s * (delta_x - hat_tau) / t
        // delta_z <- z * (hat_nu - delta_x) / g
        // delta_q <- q * (delta_w - hat_alpha) / p
        // delta_v <- v * (gamma_w - delta_w) / w
        // delta_p <- p * (gamma_q - delta_q) / q
        // delta_g <- g * (gamma_z - delta_z) / z
        // delta_t <- t * (gamma_s - delta_s) / s

		for(i = 0;i < m;i++) delta_w[i] = -e[i]*(hat_beta[i]-q[i]*hat_alpha[i]/p[i] + delta_y[i]);
		for(i = 0;i < n;i++) delta_s[i] = s[i]*(delta_x[i] - hat_tau[i]) / t[i];
		for(i = 0;i < n;i++) delta_z[i] = z[i]*(hat_nu[i] - delta_x[i]) / g[i];
		for(i = 0;i < m;i++) delta_q[i] = q[i]*(delta_w[i] - hat_alpha[i]) / p[i];
		for(i = 0;i < m;i++) delta_v[i] = v[i]*(gamma_w[i] - delta_w[i]) / w[i];		
		for(i = 0;i < m;i++) delta_p[i] = p[i]*(gamma_q[i] - delta_q[i]) / q[i];		
		for(i = 0;i < n;i++) delta_g[i] = g[i]*(gamma_z[i] - delta_z[i]) / z[i];
		for(i = 0;i < n;i++) delta_t[i] = t[i]*(gamma_s[i] - delta_s[i]) / s[i];
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////		
        // compute the updates
        // alfa <- - (1 - margin) / min(c(delta_g / g, delta_w / w, delta_t / t, delta_p / p, delta_z / z, delta_v / v, delta_s / s, delta_q / q, -1))
        // x <- x + delta_x * alfa
        // g <- g + delta_g * alfa
        // w <- w + delta_w * alfa
        // t <- t + delta_t * alfa
        // p <- p + delta_p * alfa
        // y <- y + delta_y * alfa
        // z <- z + delta_z * alfa
        // v <- v + delta_v * alfa
        // s <- s + delta_s * alfa
        // q <- q + delta_q * alfa

		alfa = -1.0;
		for(i = 0;i < n;i++){
			alfa = MIN(alfa,delta_g[i]/g[i]);
			alfa = MIN(alfa,delta_s[i]/s[i]);
			alfa = MIN(alfa,delta_t[i]/t[i]);
			alfa = MIN(alfa,delta_z[i]/z[i]);
	    }
		for(i = 0;i < m;i++){
			alfa = MIN(alfa,delta_p[i]/p[i]);
			alfa = MIN(alfa,delta_q[i]/q[i]);
			alfa = MIN(alfa,delta_v[i]/v[i]);
			alfa = MIN(alfa,delta_w[i]/w[i]);
	    }	
		alfa = -(1.0 - margin) / alfa;
		

		for(i = 0;i < n;i++) g[i] = g[i] + delta_g[i] * alfa;
		for(i = 0;i < n;i++) s[i] = s[i] + delta_s[i] * alfa;
		for(i = 0;i < n;i++) t[i] = t[i] + delta_t[i] * alfa;
		for(i = 0;i < n;i++) z[i] = z[i] + delta_z[i] * alfa;
		for(i = 0;i < m;i++) y[i] = y[i] + delta_y[i] * alfa;
		for(i = 0;i < m;i++) p[i] = p[i] + delta_p[i] * alfa;
		for(i = 0;i < m;i++) q[i] = q[i] + delta_q[i] * alfa;
		for(i = 0;i < m;i++) v[i] = v[i] + delta_v[i] * alfa;
		for(i = 0;i < m;i++) w[i] = w[i] + delta_w[i] * alfa;
		for(i = 0;i < n;i++) x[i] = x[i] + delta_x[i] * alfa;
		
		for(i = n;i < n;i++){
			if(EQ(lo[i],x[i],DBL_EPSILON)) 
				x[i] = lo[i];
			else if(EQ(up[i],x[i],DBL_EPSILON)) 
				x[i] = up[i];
		}		
///////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////				
        // // these two lines put back in ?
        // // mu <- (crossprod(z,g) + crossprod(v,w) + crossprod(s,t) + crossprod(p,q))/(2 * (m + n))
        // // mu <- mu * ((alfa - 1) / (alfa + 10))^2
        // mu <- newmu
		mu = newmu;
///////////////////////////////////////////////////////////////////////		

///////////////////////////////////////////////////////////////////////
// note: check for sigfig_max is done immediately after sigfig is computed 
		if((*counter) == maxiter){
			*convergence = SLOW_CONVERGENCE;		
			break;
		// }else if (sigfig >= sigfig_max){
            // *convergence = CONVERGED;
			// break;
		}else if ((primal_infeasibility > inf) && (dual_infeasibility > inf)){
			*convergence = PRIMAL_DUAL_INFEASIBLE;
			break;
		}else if (primal_infeasibility > inf){
			*convergence = PRIMAL_INFEASIBLE;
			break;
		}else if (dual_infeasibility > inf){
			*convergence = DUAL_INFEASIBLE;
			break;
		}else if (ABS(primal_obj) > inf){ 
			*convergence = PRIMAL_UNBOUNDED;
			break;
		}else if (ABS(dual_obj) > inf){ 
			*convergence = DUAL_UNBOUNDED;
			break;
		}
		// int infeasible_solution = 0;
		// for(i = 0;i < n;i++){
			// if(x[i] < lo[i]){ 
				// x[i] = lo[i];
				// infeasible_solution = 1;
			// }else if(x[i] > up[i]){
				// x[i] = up[i];
				// infeasible_solution = 1;				
			// }
		// }
		// if(infeasible_solution){
			// *convergence = PRIMAL_INFEASIBLE_SOLUTION;
			// break;
		// }		
///////////////////////////////////////////////////////////////////////
    }
    if (verb > 0){		      	// final report
        Rprintf("Iter    PrimalInf  DualInf  SigFigs  Rescale  PrimalObj  DualObj\n");
		Rprintf("%i    %.6e %.6e %.6e %.6e %.6e %.6e\n",
			(*counter), primal_infeasibility, dual_infeasibility, sigfig, alfa, primal_obj, dual_obj);
	}
	int infeasible_solution = 0;
	for(i = 0;i < n;i++){
		if(x[i] < lo[i]){ 
			infeasible_solution = 1;
			break;
		}else if(x[i] > up[i]){
			infeasible_solution = 1;
			break;
		}
	}
	if(infeasible_solution){
		*convergence = PRIMAL_INFEASIBLE_SOLUTION;
	}
		
	memcpy(primal,x,n * sizeof(double));
	memcpy(dual,y,m * sizeof(double));
	
	// Do not set return values to DBL_MAX.
	// DBL_MAX is R's 'Inf' and .C() call fails by default (unless NAOK = TRUE)
	if(inf == DBL_MAX) inf = 0.5*DBL_MAX;
	if(_primal_obj)
		*_primal_obj = primal_obj;
	if(_dual_obj)
		*_dual_obj = dual_obj;
		
	free(ibuffer);
	if(internal_buffer)free(dbuffer);//else memset(_dbuffer,0,*_dbuffer_size * sizeof(double));
}


/*
Buffer size  = (m+n)^2 + 10(m+n)(lwork) + 10(m+n)(as above) + 10(m+n) + (3n + 2)
In addition, updating requires more

lo      1  n
up      1  n
A       1  n
b       1
r       1
           3n + 2

KKT0       (m+n)^2
(x,y)      (m+n)
(cx,cy)    (m+n)
(d,e)      (m+n)
delta(x,y) (m+n)
diag_KKT0  (m+n)
lwork      10(m+n)
diag_H     n
H_dot_x    n
rho,tau    m+n
alpha,nu   m+n
beta,signa m+n

hat:
alpha,nu   m+n
beta,signa m+n

////////////////
g,s,t,z    4n
p,q,v,w    4m

delta:
g,s,t,z    4n
p,q,v,w    4m

gamma:
s,z        2n

gamma
q,w        2m
///////////////
10(m+n)


// p                        primal slack variable for inequality constraint
// q						dual   slack variable for inequality constraint
// t                       primal upper bound slack
// s                       dual   upper bound slack
// z                       primal lower bound slack
// g                       dual   lower bound slack

*/ 
// these below are not needed because we declare alpha-rho-nu-tau and beta-sigma contiguosly :
// buffer: primal_infeasibility m + n + m + n 
// buffer: dual_infeasibility   n + m 



// H_x_{n x n} = hessian, aside from negative sign, looks like
// (H_11 .. H_1n)
// ..
// (H_n1 .. H_nn)

// H_y_{m x m} = 
// (Y_11 .. Y_1m)
// ..
// (Y_m1 .. Y_mm)

// A_(m x n) consraint matrix
// (a_11 .. a_1n)
// ..
// (a_m1 .. a_mn)
// = (a_11..a_m1) .. (a_1n..a_mn) column-major order
// KKT0 = 
// (H_11 .. H_1n) | (a_11 .. a_m1)
// ..             ..
// (H_n1 .. H_nn) | (a_1n .. a_mn)
// -------------------------------
// (a_11 .. a_1n) | (Y_11 .. Y_1m)
// ..             ..
// (a_m1 .. a_mn) | (Y_m1 .. Y_mm)

// The matrix KKT0 needs to be filled-out in column-major order to pass to LAPACK
// But keep in mind that 'A' and 'H_y' are also stored in column-major order
// H_x = hessian, aside from negative sign, looks like
// (H_11 .. H_n1)(a_11 .. a_m1) 
// ..
// (H_1n .. H_nn)(a_1n .. a_mn) 
// (a_11 .. a_n1)(Y_11 .. a_1m) 
// ..
// (a_m1 .. a_mn)(Y_m1 .. Y_mm) 

// 'if-logic' references the conceptual matrix 'KKT0'
// on the other hand indexing is for column-major order

// KKT0 conceptually :
// KKT0 = [-H_x  A' ]
//      [ A    H_y]         
// in column-major order:
// KKT0 = [-H_x' A   ]
//      [ A'   H_y'] 

