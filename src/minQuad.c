// all headers common for C and/or C++ code
#include "minQuad_QP.h"
#include "minQuad.h"
#include "working_set.h"
#include "qp_exhaustive.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <R_ext/Lapack.h>
#include <R.h>
#include <Rinternals.h>

#include <math.h>
#include <float.h>
#define _Inf DBL_MAX
#define	MAX(A,B)	((A) > (B) ? (A) : (B))
#define	MIN(A,B)	((A) < (B) ? (A) : (B))
#define PROJ(L,X,U) MAX((L),MIN((X),(U))) // median(L,X,U)
#define intCEIL(X,Y) (1 + (((X) - 1) / (Y)))  
#define	ABS(A)  	((A) > 0 ? (A) : (-(A)))

// ceil(x/y) where 'x' and 'y' are integers
// equivalent to  (x + y - 1) / y where '/' stands for integer division
// it also avoids overflow in y
// http://stackoverflow.com/questions/2745074/fast-ceiling-of-an-integer-division-in-c-c 
#define PRINTF Rprintf

///////////////////////////////////////////////////////////////////////////////
/////////////////// Can define or comment out the following ///////////////////
// #define DO_VLA 
// a very bad idea to define this for any n = n1 + n2 s.t. the size of 'x' in 
// double x[n1*n2] is greater than a few kilobytes, unless your system has
// large amount of memory dedicated to the STACK

// Kahan-summing for updating 'Ha' 
//#define DO_KAHAN

///////////////////////////////////////////////////////////////////////////////

//#define CHECK_EXHAUSTIVE_SOLUTION
 

//#define _DEBUG
#ifdef _DEBUG
#include <errno.h>
FILE* _F = NULL;
#define _DEBUG_FILENAME "aucm.debug"
#define FPRINTF fprintf
#endif 


//#define MEASURE_ELAPSED_TIME
#ifdef MEASURE_ELAPSED_TIME
#include <time.h>
#endif

// quadratic program solver 'TRON'          
extern void optimize_qp_tron(int* n,double* hessian,double* linear,double* x_init,double* x_lo,double* x_up,double* solution,double* control_tron,int* verbose, int* convergence);   
extern int optimize_qp_loqo(MNQD_QP* qp, double* solution,double* control,int verbose);               


// GLOBAL VARIABLES to be freed if minQuad exits prematurely
typedef struct ptr_mnqd{
double** H;                  // The Hessian, a '**' for 'Q' or 'K' passed from R
double *solution;            // solution to quadratic problem of size of working set
double *a0;                  // vector of alpha's at previous iteration
double *Ha;                  // cached Qa
double *gradient;            // derivative of .5a'Qa + b'a
#ifdef DO_KAHAN
double *Ha_bits;             // Kahan sum correction bits for Ha
double *c_bits;              
#endif
double *c;                   
//double *hypercube;
int *working_set;           // Working set at each iteration if requested - may be memory intensive
int *violators;              // type of violator v[i] = 0,1,2 for i = 0,..,n1n2-1    
double* dbuf;               
int* ibuf;
int* ibuf2;
}PTR_MNQD;
PTR_MNQD state;

MNQD_QP main_qp;            // quadratic sub-problem corresponding to current working set
MNQD_QP sub_qp;             // 'a' quadratic problem that is a sub-problem of the quadratic sub-problem that corresponds to current working set
                            // sub_qp is only used by exhaustive method
							
//double machine epsilon, smallest number such that 1.0+x == x is true
void get_machine_double_eps(double* eps){            
    *eps = 1.0;
    do{*eps /= 2.0;}while ((1.0 + *eps) != 1.0);
}	
	
void free_minQuad(){
#ifdef _DEBUG
FPRINTF(_F,"free_minQuad\n");
#endif
//	PRINTF("free_minQuad\n");
    if(state.H){free(state.H);state.H = NULL;}
//    if(state.hypercube){free(state.hypercube);state.hypercube = NULL;}
    if(state.dbuf){free(state.dbuf);state.dbuf = NULL;}
    if(state.ibuf){free(state.ibuf);state.ibuf = NULL;}
    if(state.ibuf2){free(state.ibuf2);state.ibuf2 = NULL;}
   //PRINTF("free eps\n");
//  if(state.eps_history){free(state.eps_history);state.eps_history = NULL;}
//    PRINTF("free B\n");
    if(state.working_set){free(state.working_set);state.working_set = NULL;}
//    PRINTF("free aB\n");
    if(state.solution){free(state.solution);state.solution = NULL;}
//    PRINTF("free v\n");
    if(state.violators){free(state.violators);state.violators = NULL;}
//    PRINTF("free c\n");
    if(state.c){free(state.c);state.c = NULL;}
//    PRINTF("free a0\n");
    if(state.a0){free(state.a0);state.a0 = NULL;}
//    PRINTF("free Ha\n");
    if(state.Ha){free(state.Ha);state.Ha = NULL;}
//    PRINTF("free df\n");
    if(state.gradient){free(state.gradient);state.gradient = NULL;}
#ifdef DO_KAHAN	
    //PRINTF("free Ha_bits\n");
    if(state.Ha_bits){free(state.Ha_bits);state.Ha_bits = NULL;}
    //PRINTF("free cbits\n");
    if(state.c_bits){free(state.c_bits);state.c_bits = NULL;}
#endif
    memset(&state,0,sizeof(state));
//    PRINTF("free main and sub qp\n");
    free_QP(&main_qp);
    free_QP(&sub_qp);    
}


   
void dsignif(int n,double* x,int ndigits,double eps){
    double y = pow(10.0,(double)ndigits);
    for(int i =0;i < n;i++){x[i] = (x[i] * y + 0.5) / y;}
}   
   
// detect a repeating pattern (cycle) of length 'c' in 'x' of length 'n'
// pass 'eps = 0.0' for exact matching  
int IsCycle(int n, double* x, int c,double eps){
    int k = n / c; //integer division !
    int rem = n - k;
    x += rem;      // align x or expand x = c(x,x[1:(l-r)])
    n -= rem;
    for(int i = 0;i < k;i++)
        for(int j = 0;j < c;j++)
            if(ABS(x[i * c + j ] - x[j]) > eps) return 0;
    return 1;
}


// MAY BE INACCURATE for large 'n' - should add Kahan sum 
// DATA_PTR_TYPE type: unfortunately dynamic (run-time) typecasting does not exist in C
double objective(double machine_double_eps, int n1,int n2,int n1n2,DATA_PTR_TYPE type,void* _Q,double* a,double* b){

#ifdef _DEBUG
FPRINTF(_F,"objective()\n");
#endif

    double v = 0.0;  
    if(type == SINGLE_Q) 
    {
#ifdef DO_VLA
        double (*Q)[n1n2] = (double (*)[n1n2])(_Q);
#else		
        double* Q = (double*)_Q;
#endif		
        for(int p = 0; p < n1n2; p++) {
            if(ABS(a[p]) <= machine_double_eps) continue;
            for(int q = 0; q < n1n2; q++) {
                if(ABS(a[q]) <= machine_double_eps) continue;
#ifdef DO_VLA
                v   +=  Q[p][q] * a[p] * a[q];
#else				
				v   +=   Q[p*n1n2 + q]  * a[p] * a[q];
#endif			  
            }
        }        
    }else if(type == DOUBLE_Q)
    {
        double** Q = (double**) _Q;
        for(int p = 0; p < n1n2; p++) {
            if(ABS(a[p]) <= machine_double_eps) continue;
            for(int q = 0; q < n1n2; q++) {
				if(ABS(a[p]) <= machine_double_eps) continue;
                v   +=  Q[p][q] * a[p] * a[q];
            }
        }   
    }else if(type == SINGLE_K){
#ifdef DO_VLA		
        double (* K)[n1+n2] = (double (*)[n1+n2])_Q;
#else
        int n = n1+n2;
        double* K = (double*)_Q;
#endif        
        for(int p = 0; p < n1n2; p++) {
            if(ABS(a[p]) <= machine_double_eps) continue;
            int ip=intCEIL(p+1,n2);//int ip=(int)ceil((p+1.0)/n2); 
            int jp=(p+1)-(ip-1)*n2;
            for(int q = 0; q < n1n2; q++) {
				if(ABS(a[p]) <= machine_double_eps) continue;
                int iq=intCEIL(q+1,n2);//int iq=(int)ceil((q+1.0)/n2); 
                int jq=(q+1)-(iq-1)*n2;     
#ifdef DO_VLA		
                v += (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]) * a[p] * a[q];
#else
                v += a[p]*a[q]*(K[(ip-1)*n+(iq-1)]+K[(jp+n1-1)*n+(jq+n1-1)]-K[(ip-1)*n+(jq+n1-1)]-K[(jp+n1-1)*n+(iq-1)]);                                 
#endif				
            }
        }
    }else if(type == DOUBLE_K){
        double** K = (double**)_Q;
        for(int p = 0; p < n1n2; p++) {
			if(ABS(a[p]) <= machine_double_eps) continue;
            int ip=intCEIL(p+1,n2);//int ip=(int)ceil((p+1.0)/n2); 
            int jp=(p+1)-(ip-1)*n2;
            for(int q = 0; q < n1n2; q++) {
				if(ABS(a[p]) <= machine_double_eps) continue;
                int iq=intCEIL(q+1,n2);//int iq=(int)ceil((q+1.0)/n2); 
                int jq=(q+1)-(iq-1)*n2;                     
                v += (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]) * a[p] * a[q];
            }
        }
    }else{ PRINTF("Error in objective(): invalid type for 'Q'.");return(sqrt(-1.0));}   
    
    v *= 0.5;
    for(int p = 0; p < n1n2; p++) v += a[p] * b[p];
    return v;
}

// to be called from R
// result = 0.5*a'Ha + b'a
// 'H' is specified by 'type' 
// 'type' must match one of 'SINGLE_Q' or 'SINGLE_K' values defined in minQuad.h: typedef enum {SINGLE_Q = 0,SINGLE_K,DOUBLE_Q,DOUBLE_K}DATA_PTR_TYPE;
void R_objective(double* machine_double_eps,int* _n1,int* _n2,int* _n1n2,int* type,double* _Q,double* a,double* b,double* result){
	*result = objective(*machine_double_eps,*_n1,*_n2,*_n1n2,*type,(void*)_Q,a,b);
}


//// Q_pred is passed from R as {n1n2 x npred}
//void get_Q_pred(
//    //input
//    double* _K,
//    int* _n1,
//    int* _n2,
//    int* _n_pred,
//    double* _Q//output
//)
//{
//#ifdef _DEBUG
//FPRINTF(_F,"get_Q_pred()\n");
//#endif
//    int n1 = *_n1;
//    int n2 = *_n2;
//    int n_pred = *_n_pred;
//    int n1n2 = n1*n2;    
//    
//#ifdef DO_VLA        
//    double (*K)[n_pred] = (double (*)[n_pred])_K;	
//    double (*Q)[n_pred] = (double (*)[n_pred])_Q;	
//#endif
//    
//    int i, q;
//    for(i = 0; i < n_pred; i++) {
//        for(q = 0; q < n1n2; q++) {
//            int iq=intCEIL(q+1,n2);//int iq=(int)ceil((q+1.0)/n2);         
//            int jq=(q+1)-(iq-1)*n2;                     
//#ifdef DO_VLA        
//            Q[q][i] = K[iq-1][i] - K[jq+n1-1][i];
//#else			
//            _Q[q*n_pred + i] = _K[(iq-1)*n_pred + i] - _K[(jq+n1-1)*n_pred + i];
//#endif
//        }
//    }
//			PRINTF("%f\n",_Q[1]);
//}


// K is n.pred x (n1+n2)
// return a npred x n1n2 matrix
SEXP get_Q_pred(
    SEXP _K,
    SEXP _n1,
    SEXP _n2
)
{
 	int n_pred = nrows(_K);
	int n1 = *INTEGER(_n1);
	int n2 = *INTEGER(_n2);
    int n1n2 = n1*n2;    
	
    double *K=REAL(_K);
    
    SEXP _ans = PROTECT(allocMatrix(REALSXP, n_pred, n1n2));
    double *ans=REAL(_ans);
    
    int i, q;
    for(i = 0; i < n_pred; i++) {
        for(q = 0; q < n1n2; q++) {
            int iq=intCEIL(q+1,n2);//int iq=(int)ceil((q+1.0)/n2);         
            int jq=(q+1)-(iq-1)*n2;                     
            ans[q*n_pred + i] = K[(iq-1)*n_pred + i] - K[(jq+n1-1)*n_pred + i];
        }
    }

    UNPROTECT(1);
    return _ans;
}

// linear combination = Q_pred_{npred x n1n2} %*% x_{n1n2 x 1}
// entries in 'Q_pred' are implicitely computed via 'K', see function get_Q_pred() 
// order of nested loops switched since 'x' may have lots of zeros


SEXP get_Q_pred_x(
    SEXP _K,
    SEXP _n1,
    SEXP _n2,
	SEXP _x,
	SEXP _machine_double_eps
)
{
 	int n_pred = nrows(_K);
	int n1 = *INTEGER(_n1);
	int n2 = *INTEGER(_n2);
    int n1n2 = n1*n2;    
	double machine_double_eps = *REAL(_machine_double_eps);
	
    double *K=REAL(_K);
    double *x=REAL(_x);
    
    SEXP _ans = PROTECT(allocVector(REALSXP, n_pred));
    double *ans=REAL(_ans);
    
    int i, q;
    for(i = 0; i < n_pred; i++) {
    	ans[i] = 0; 
        for(q = 0; q < n1n2; q++) {
            if(ABS(x[q]) <= machine_double_eps) continue;
            int iq=intCEIL(q+1,n2);//int iq=(int)ceil((q+1.0)/n2);         
            int jq=(q+1)-(iq-1)*n2;                     
            ans[i] += x[q] *(K[(iq-1)*n_pred + i] - K[(jq+n1-1)*n_pred + i]);
        }
    }

    UNPROTECT(1);
    return _ans;
}





void get_Qx(
    //input
	double* _machine_double_eps,
    double* _K,
    int* _n1,
    int* _n2,
    double* x,
    double* Qx//output
)
{
    int n1   = *_n1;
    int n2   = *_n2;
    int n    = n1+n2;
    int n1n2 = n1*n2;
    int p, q;

	
#ifdef DO_VLA        
    double (*K)[n] = (double (*)[n])_K;	
#endif    
    double machine_double_eps = *_machine_double_eps;
	memset(Qx,0.0,n1n2 * sizeof(double));

	for(q = 0; q < n1n2; q++) {
		if(ABS(x[q]) <= machine_double_eps) continue; 
		int iq=intCEIL(q+1,n2);//int iq=(int)ceil(((double)q+1.0)/(double)n2);
		int jq=(q+1)-(iq-1)*n2; 
	    int jq_n1_1 = jq+n1-1;
		for(p = 0; p < n1n2; p++) {
			int ip=intCEIL(p+1,n2);
			int jp =(p+1)-(ip-1)*n2;  
			int jp_n1_1_n = (jp+n1-1)*n;
            int ip_1_n = (ip-1)*n;
#ifdef DO_VLA                
            Qx[p] += x[q] * (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]);
#else
//            Qx[p] += x[q] * (_K[(ip-1)*n + iq-1]+_K[(jp+n1-1)*n + jq+n1-1]-_K[(ip-1)*n + jq+n1-1]-_K[(jp+n1-1)*n + iq-1]);
            Qx[p] += x[q] * (_K[ip_1_n + iq-1]+_K[jp_n1_1_n + jq_n1_1]-_K[ip_1_n + jq_n1_1]-_K[jp_n1_1_n + iq-1]);
			
#endif
        }
    }  
}

// // original get_Q() - order of nested for-loops reversed
// void get_Qx1(
    // //input
    // double* machine_double_eps,
    // double* _K,
    // int* _n1,
    // int* _n2,
    // double* x,
    // //output
    // double* Qx,
// )
// {
    // int n1   = *_n1;
    // int n2   = *_n2;
    // int n    = n1+n2;
    // int n1n2 = n1*n2;
    // int p, q;

// #ifdef DO_VLA        
    // double (*K)[n] = (double (*)[n])_K;	
// #endif    
//    double machine_double_eps = *machine_double_eps;

    // for(p = 0; p < n1n2; p++) {
        // Qx[p] = 0.0;
        // int ip=intCEIL(p+1,n2);//int ip=(int)ceil(((double)p+1.0)/(double)n2)
        // int jp =(p+1)-(ip-1)*n2;  

        // for(q = 0; q < n1n2; q++) {
            // if(ABS(x[q]) <= machine_double_eps) continue; 
			// int iq=intCEIL(q+1,n2);//int iq=(int)ceil(((double)q+1.0)/(double)n2);
			// int jq=(q+1)-(iq-1)*n2;  

// #ifdef DO_VLA                
            // Qx[p] += x[q] * (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]);
// #else
            // Qx[p] += x[q] * (_K[(ip-1)*n + iq-1]+_K[(jp+n1-1)*n + jq+n1-1]-_K[(ip-1)*n + jq+n1-1]-_K[(jp+n1-1)*n + iq-1]);
// #endif
        // }
    // }  
// }

void get_Q(
    //input
    double* _K,
    int* _n1,
    int* _n2,
    double* _Q//output
)
{

#ifdef _DEBUG
FPRINTF(_F,"get_Q()\n");
#endif

    int n1   = *_n1;
    int n2   = *_n2;
    int n    = n1+n2;
    int n1n2 = n1*n2;
    

#ifdef DO_VLA    
    double (*K)[n]   = (double (*)[n])_K;	
    double (*Q)[n1n2]   = (double (*)[n1n2])_Q;	
#endif
    int p, q;
    for(p = 0; p < n1n2; p++) {
        int ip=intCEIL(p+1,n2);//int ip=(int)ceil((p+1.0)/n2);         
        int jp=(p+1)-(ip-1)*n2;
        for(q = 0; q < n1n2; q++) {
            int iq=intCEIL(q+1,n2);//int iq=(int)ceil((q+1.0)/n2);         
            int jq=(q+1)-(iq-1)*n2;                     
#ifdef DO_VLA    
            Q[p][q] = K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1];
#else
           _Q[p*n1n2+q]=_K[(ip-1)*n + iq-1]+_K[(jp+n1-1)*n + jq+n1-1]-_K[(ip-1)*n + jq+n1-1]-_K[(jp+n1-1)*n + iq-1];            
#endif
        }
    }
}


void optimize_1(double machine_double_eps,int ns,int* subset,double* C,int n1n2,int n1,int n2,DATA_PTR_TYPE type,void* H,double* b,double* c,double* kahan_bits,double* a)
{
    double Q_pp,Q_pq;
    int p,q;
    int N = n1n2;if(subset) N = ns;

    if(type == SINGLE_Q){
        double (*Q)[n1n2]   = (double (*)[n1n2])H;	
        for(int i = 0;i < N;i++){
            p = i;if(subset) p = subset[i];
            Q_pp = Q[p][p];
    ///////////////////////////////////////////////////////
    // new a[i] = median({0, C, (0.5b[i] - c[i]) / Kii} )
            double da_p = -a[p];
//            a[p] = (0.5*b[p] - c[p]) / Q_pp;
            a[p] = -(b[p] + c[p]) / Q_pp;
            a[p] = PROJ(0.0,a[p],C[p]);
    ///////////////////////////////////////////////////////
            da_p += a[p];
            if(ABS(da_p) <= machine_double_eps) continue;       
            
            for(q = 0; q < n1n2; q++){
                if(p == q) continue;
    #ifdef DO_KAHAN				
                double volatile kahan_summand = Q[p][q] * da_p;
                double volatile kahan_z = kahan_summand - kahan_bits[p];
                double volatile kahan_sum = c[p] + kahan_z;
                kahan_bits[p] = kahan_sum - c[p] - kahan_z;
                c[p] = kahan_sum;   
	#else
				c[p] += Q[p][q] * da_p;
    #endif            
            }
        }
    }else if(type == SINGLE_K){
        double (*K)[n1+n2]   = (double (*)[n1+n2])H;	
        for(int i = 0;i < N;i++){
            p = i;if(subset) p = subset[i];
            int ip=intCEIL(p+1,n2);//int ip=(int)ceil((p+1.0)/n2);         
            int jp=(p+1)-(ip-1)*n2;
            Q_pp = K[ip-1][ip-1]+K[jp+n1-1][jp+n1-1]-K[ip-1][jp+n1-1]-K[jp+n1-1][ip-1];
            
    ///////////////////////////////////////////////////////
    // new a[i] = median({0, C, (0.5b[i] - c[i]) / Kii} )
            double da_p = -a[p];
//          a[p] = (0.5*b[p] - c[p]) / Q_pp;
            a[p] = -(b[p] + c[p]) / Q_pp;
            a[p] = PROJ(0.0,a[p],C[p]);
    ///////////////////////////////////////////////////////
            da_p += a[p];
            if(ABS(da_p) <= machine_double_eps) continue;    
            
            for(q = 0; q < n1n2; q++) {
                if(p == q) continue;
                int iq=intCEIL(q+1,n2);//int iq=(int)ceil((q+1.0)/n2);         
                int jq =(q+1)-(iq-1)*n2;                     
                Q_pq = K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1];
    #ifdef DO_KAHAN				                
                double volatile kahan_summand = Q_pq * da_p;
                double volatile kahan_z = kahan_summand - kahan_bits[p];
                double volatile kahan_sum = c[p] + kahan_z;
                kahan_bits[p] = kahan_sum - c[p] - kahan_z;
                c[p] = kahan_sum;   
    #else
                c[p] += Q_pq * da_p;
	#endif
            }        
        }
    }
}

//c[i] = < Q[i,-i] , a[-i] > = < Q_i. , a > - Q_ii*a[i] 
void get_c_solveQuad(double machine_double_eps,int n1n2,int n1,int n2,DATA_PTR_TYPE type,void* _A,double* a,double* kahan_bits,double* c)
{
    if(type == SINGLE_Q){
        double (*Q)[n1n2] = (double (*)[n1n2])_A;	
        for(int p = 0;p < n1n2;p++){                                                             
            for(int q = 0; q < n1n2; q++){
                if(p == q) continue;
    #ifdef DO_KAHAN				                
                double volatile kahan_summand = Q[p][q] * a[q];
                double volatile kahan_z = kahan_summand - kahan_bits[p];
                double volatile kahan_sum = c[p] + kahan_z;
                kahan_bits[p] = kahan_sum - c[p] - kahan_z;
                c[p] = kahan_sum;   
	#else
                c[p] += Q[p][q] * a[q];
	#endif
            }
            
        } 
    }else if(type == SINGLE_K){
        double (*K)[n1+n2] = (double (*)[n1+n2])_A;	
        for(int p = 0;p < n1n2;p++){
            int ip=intCEIL(p+1,n2);//int ip=(int)ceil((p+1.0)/n2);         
            int jp=(p+1)-(ip-1)*n2;  
            for(int q = 0; q < n1n2; q++) {
                if(q == p) continue;
                int iq=intCEIL(q+1,n2);//int iq=(int)ceil((q+1.0)/n2);         
                int jq=(q+1)-(iq-1)*n2;                     
                double Q_pq = K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1];
    #ifdef DO_KAHAN				                
                double volatile kahan_summand = Q_pq * a[q];
                double volatile kahan_z = kahan_summand - kahan_bits[p];
                double volatile kahan_sum = c[p] + kahan_z;
                kahan_bits[p] = kahan_sum - c[p] - kahan_z;
                c[p] = kahan_sum;   
	#else
                c[p] += Q_pq * a[q];
	#endif
            }
        }    
    }
}

//c[i] = < Q[i,-i] , a[-i] > = < Q_i. , a > - Q_ii*a[i] 
void update_c_solveQuad(double machine_double_eps,int ns,int* subset,int n1n2,int n1,int n2,DATA_PTR_TYPE type,void* _A,double* a0,double* a,double* kahan_bits,double* c)
{
    int p,q;
    int N = n1n2;if(subset) N = ns;
        
    if(type == SINGLE_Q){    
        double (*Q)[n1n2] = (double (*)[n1n2])_A;	
        for(int i = 0;i < N;i++){
            if(subset) p = subset[i];else p = i;
            double da_p = a[p] - a0[p];
            if(ABS(da_p) <= machine_double_eps) continue;              
            for(q = 0; q < n1n2; q++){
                if(q == p) continue;
    #ifdef DO_KAHAN				                
                double volatile kahan_summand = Q[p][q] * da_p;
                double volatile kahan_z = kahan_summand - kahan_bits[p];
                double volatile kahan_sum = c[p] + kahan_z;
                kahan_bits[q] = kahan_sum - c[q] - kahan_z;
                c[p] = kahan_sum;                                            
	#else
                c[p] += Q[p][q] * da_p;
	#endif
            }
        }
    }else if(type == SINGLE_K){
        double (*K)[n1+n2]   = (double (*)[n1+n2])_A;	    
        for(int i = 0;i < N;i++){
            if(subset) p = subset[i];else p = i;
            double da_p = a[p] - a0[p];
            if(ABS(da_p) <= machine_double_eps) continue;
            int ip=intCEIL(p+1,n2);//int ip=(int)ceil((p+1.0)/n2);         
            int jp=(p+1)-(ip-1)*n2; 
            for(q = 0; q < n1n2; q++) {
                if(p == q) continue;
                int iq=intCEIL(q+1,n2);//int iq=(int)ceil((q+1.0)/n2);         
                int jq=(q+1)-(iq-1)*n2;                     
                double Q_pq = K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1];       
    #ifdef DO_KAHAN				                
                double volatile kahan_summand = Q_pq * da_p;
                double volatile kahan_z = kahan_summand - kahan_bits[p];
                double volatile kahan_sum = c[p] + kahan_z;
                kahan_bits[p] = kahan_sum - c[p] - kahan_z;
                c[p] = kahan_sum;               
    #else
                c[q] += Q_pq * da_p;
    #endif	
            }
        }
    }    
}

void minQuad(
    //input
    int* _n1,
    int* _n2,
    int* _n1n2,        // n1n2 = n_diseased * n_non_diseased from R
    int* mem_efficient,//  1/0 <-> _H = K_(n1+n2 x n1+n2) or _H = Q_n1n2xn1n2 
    double* _H,        
    double* b,         // length n1n2
    double* C,         // length n1n2
    double* _a,        // initial vector of alpha's, length n1n2
	int* n_constr,  // whether to do the constraint below, to be satisfied for each sub-problem, pass 0 if no constraint
    double* mat_constr, // (n_constr x n1n2): lhs <= M*x <= rhs
	double* lhs_constr, // length n_constr
	double* rhs_constr, // length n_constr
    int* _rank,
    int* _s,           // do a sweep via solveQuad every s-th iterations  
    int* _q,           // size of working set
    int* working_set_method,   // heuristics for working set selection 0,1,2,3 for fixed,random 1,2,3
    int* return_working_set, //save working-set at each iteration  up to 'maxit'?
    int* optim,              // which optimizer
	double* optim_control,
    int* maxit,        // max no. of iterations, >= 1
    double* tol,       // tolerance criterion
	double* _machine_double_eps,
    int* _trace,       // 0,1,2: print at each iteration (2), before returning to caller (1) no print(0)   
    //output
    double* a,         // return final alpha of length 'n1n2' to caller
    double* obj,       // return value of objective function to caller
    int* iter,         // return the final number of iteration to caller
    int* convergence,  // did alg. converge (Y/N) <-> (0/1)
    double*  epsilon,  // return value in to termination condition to caller
    int* working_set   // save ((B_1),(B_2),..,(B_maxit)) B_k = (p,q), length = 2 * maxit
)
{

#ifdef _DEBUG
_F = fopen(_DEBUG_FILENAME,"w+");
if(!_F){
	int errnum=errno;
	PRINTF("minQuad() could not open debug file, reason: %s\n",strerror(errnum));
	return;
}
FPRINTF(_F,"minQuad() initialization\n");
#endif

// #ifdef _DEBUG
// PRINTF("minQuad quadratic optimizer\n");
// #endif

    #ifdef MEASURE_ELAPSED_TIME
    clock_t t0_init = clock();
    clock_t t1_init,t0_loop,t1_loop,t0_end,t1_end;
    clock_t t0_init_cache,t1_init_cache;
    #endif
    
    int trace = *_trace;
    int p,q,q1,ip,jp,iq,jq;
    int n1n2 = *_n1n2;
    int n1 = *_n1;
    int n2 = *_n2;
    int done = 0;
    int UPDATE_C = 0; //disables certain logic 

	double machine_double_eps = *_machine_double_eps;
	if(machine_double_eps + 1.0 != 1.0)get_machine_double_eps(&machine_double_eps);
    
	
    char* QP_METHOD_NAMES[] = {"tron","loqo","exhaustive","exhaustive2"};
    typedef enum {TRON = 0,LOQO = 1,EXHAUSTIVE = 2,EXHAUSTIVE2 = 3}QP_METHOD_TYPE;
    
    char* WORKINGSET_METHOD_NAMES[]       = {"rvwg","rv2wg","rv","rv2","v","v2","greedy"}; 
    
    typedef enum {RANDOM_VL_WEIGHTED_GRAD = 0,RANDOM_VL2_WEIGHTED_GRAD = 1,RANDOM_VL = 2,
				  RANDOM_VL2 = 3, VL = 4,VL2 = 5,GREEDY = 6}WS_METHOD_TYPE;
    

    if(trace > 2){
        PRINTF("quadratic optimizer(%s) working set heuristics(%s)\n",
        QP_METHOD_NAMES[*optim],WORKINGSET_METHOD_NAMES[*working_set_method]);
        PRINTF("maxit(%i) n1n2(%i), n1(%i), n2(%i), dim(working set(%i)), bounds(%f,%f)\n",
                *maxit,*_n1n2,*_n1,*_n2,*_q,0.0,C[0]);
        PRINTF("Sweep every %ith iterations\n",*_s);
        PRINTF("machine_double_eps %.20e\n",machine_double_eps);
    }
    
//  optimizing (a subset of) all variables, one at a time     
    int s = *_s;if(s == 0) s = *maxit + 1;//never occurs
     
// working set variables    
    int nV  = 0;      // total number of a's which violate the KKT conditions, nV = nV0 + nVC
    int nV0 = 0;      // total number of a's which violate the KKT conditions w.r.t. '0'
    int nVC = 0;      // total number of a's which violate the KKT conditions w.r.t. 'C'
    int dimB = *_q;
	int dynamic = 0;
	if((QP_METHOD_TYPE)(*optim) == TRON)
		if(*_q == 0)
			dynamic = 1;
	if(dynamic) dimB = MAX(dimB,(int)ceil(sqrt((double)n1n2)));
	int nB = dimB;    // size of actual working set that can be selected at kth iteration, MIN(nV,dimB)


	
	
    memset(&state,0,sizeof(state));// store pointers
//    double* hypercube = NULL;//get_hypercube(dimB,0.0,C[0]);state.hypercube = hypercube;
    state.dbuf = (double*)calloc(n1n2 , sizeof(double));if(!state.dbuf){free_minQuad();return;}
    state.ibuf = (int*)calloc(n1n2 , sizeof(int));if(!state.ibuf){free_minQuad();return;}
    state.ibuf2 = (int*)calloc(n1n2 , sizeof(int));if(!state.ibuf2){free_minQuad();return;}
    // working set    
    // int _B[(unsigned int)dimB];
    // int* B = &_B[0];      
    state.working_set = (int*)calloc(dimB , sizeof(int));if(!state.working_set){free_minQuad();return;}
    int* B = state.working_set ;
    // double _aB[(unsigned int)dimB];
    // double* aB = &_aB[0];
    state.solution = (double*)calloc(dimB , sizeof(double));if(!state.solution){free_minQuad();return;}
    double* aB = state.solution;
    // vector of violators {0,1,2} <-> {not, 0, C}    
    // int _v[(unsigned int)n1n2];  
    // int* v = &_v[0];    
    state.violators = (int*)calloc(n1n2 , sizeof(int));if(!state.violators){free_minQuad();return;}
    int* v = state.violators;
    // Caching c_i = < K_iB, a_B > in 'solveQuad' for 1 variable working set
    // double _c[(unsigned int)n1n2];    // double* c = &_v[0];
    double* c = NULL;
    if(UPDATE_C){
        state.c = (double*)calloc(n1n2 , sizeof(double));if(!state.c){free_minQuad();return;}
        c = state.c;
    }
    // alpha's from previous iteration    
    // double _a0[(unsigned int)n1n2]; 
    // double* a0 = &_a0[0];
    state.a0 = (double*)calloc(n1n2 , sizeof(double));if(!state.a0){free_minQuad();return;}
    double* a0 = state.a0;
    // double _df[(unsigned int)n1n2]; // derivative Ha + b
    // double* df = &_df[0];   
    state.gradient = (double*)calloc(n1n2 , sizeof(double));if(!state.gradient){free_minQuad();return;}
    double* df = state.gradient;
    // caching Ha = Q %*% alpha
    // double _Ha[(unsigned int)n1n2]; 
    // double* Ha = &_Ha[0];
    state.Ha = (double*)calloc(n1n2 , sizeof(double));if(!state.Ha){free_minQuad();return;}
    double* Ha = state.Ha;
    // Kahan Sum bits for avoiding accumulation of error in Ha    
    // double _kahan_bits[(unsigned int)n1n2];
    // double* kahan_bits = &_kahan_bits[0];
    // memset(kahan_bits,0,n1n2*sizeof(double)); //kahan sum bits for Ha
#ifdef DO_KAHAN			
    state.Ha_bits = (double*)calloc(n1n2 , sizeof(double));if(!state.Ha_bits){free_minQuad();return;}
    double* Ha_bits = state.Ha_bits;
    if(UPDATE_C){state.c_bits = (double*)calloc(n1n2 , sizeof(double));if(!state.c_bits){free_minQuad();return;}}
    double* c_bits = state.c_bits;
#else
	double* c_bits = NULL;
#endif
    
    void* H = NULL;
    double** Q = NULL; // used in main loop
    DATA_PTR_TYPE type_H;    
    if(!*mem_efficient) {
        if(n1n2 <= 0){PRINTF("Cannot allocate size n1n2 = 0 memory.\n");done = 1;}
        state.H = (double**) calloc(n1n2,sizeof(double*));
        if(!state.H){PRINTF("minQuad: Unable to allocate %lld(bytes)\n",n1n2*sizeof(double*));free_minQuad();return;}
        Q = state.H;
        for(p = 0; p < n1n2; p++)Q[p] = _H + p*n1n2;
        H = (void*)_H;
        type_H = SINGLE_Q;
    }else{
        state.H = (double**) calloc(n1+n2,sizeof(double*));if(!state.H){free_minQuad();return;}
        Q = state.H;        
        for(p = 0; p < (n1+n2); p++)Q[p] = _H + p*(n1+n2);
        H = (void*)_H;
        type_H = SINGLE_K;
    }
#ifdef CHECK_EXHAUSTIVE_SOLUTION
     double aB_check[dimB];
#endif

//  initialize working set        
    for(p = 0; p < dimB; p++) B[p] = p;    
//  initialize 'quadratic problem' structures        
    if(!init_QP(dimB,*n_constr,machine_double_eps,&main_qp)){free_minQuad();return;}
    if(!init_QP(dimB,*n_constr,machine_double_eps,&sub_qp)){free_minQuad();return;}


// Initialization in C instead of R if .C(...,DUP = FALSE)
    *obj = 0.0;
    *convergence = 1;
    *iter = 0;
    *epsilon = 2.0 * (*tol);    

// DUP = FALSE: problem here using memcpy() if 'alpha' from a previous call to minQuad in R 
// is passed as initial value ('_a'), unless local variables are cleared in caller (in R)
    
    for(p = 0; p < n1n2; p++) a[p] = _a[p];
    for(p = 0; p < n1n2; p++)a0[p] = _a[p];
    
///////////////////////////// Caching Initialization //////////////////////////////////////    
#ifdef MEASURE_ELAPSED_TIME    
    t0_init_cache = clock();
#endif  
// Caching Ha
    if(!*mem_efficient) {
    // LAPACK recompute Ha        
//           memset(&Ha[0], 0, sizeof(double)*n1n2);
//           double one = 1.0;int inc = 1;F77_CALL(dsymv)("U", &n1n2, &one, _H, &n1n2, a, &inc, &one, &Ha[0], &inc);


    // Kahan sum       
        for(p = 0; p < n1n2; p++) {
            Ha[p] = 0.0;
            for(q = 0; q < n1n2; q++){
#ifdef DO_KAHAN			
                double volatile y = Q[p][q] * a[q] - Ha_bits[p];
                double volatile t = Ha[p] + y;
                Ha_bits[p] = t - Ha[p] - y;
                Ha[p] = t;                    
#else
                Ha[p] += Q[p][q] * a[q];
#endif				
            }
        }
       
    } else {
//        double (*K)[(unsigned int)(n1+n2)] = (double (*)[(unsigned int)(n1+n2)])_H;	
        double** K = Q;
        for(p = 0; p < n1n2; p++) {
            ip=intCEIL(p+1,n2);//ip=(int)ceil((p+1.0)/n2);         
            jp=(p+1)-(ip-1)*n2; 
            for(q = 0; q < n1n2; q++) {
                iq=intCEIL(q+1,n2);//iq=(int)ceil((q+1.0)/n2);         
                jq=(q+1)-(iq-1)*n2;   
#ifdef DO_KAHAN			              
//			   Ha[p] = KAHAN_ADD(Ha[p],a[q] * (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]),Ha_bits[p])
               double volatile y = a[q] * (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]) - Ha_bits[p];
               double volatile t = Ha[p] + y;                
               Ha_bits[p] = t - Ha[p] - y;
               Ha[p] = t;                    
#else
               Ha[p] += a[q] * (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]);
#endif				
            }
        }
    }    
// initialize c_i = < Q_[i,-i], a[-i] >       
    if(UPDATE_C) get_c_solveQuad(machine_double_eps,n1n2,n1,n2,type_H,H,a,c_bits,(double*)&c[0]); 
#ifdef MEASURE_ELAPSED_TIME
    t1_init_cache = clock();
#endif

//  *epsilon = getEpsilonKKTV(n1n2,a,(double*)&df[0],C);  
//  if(*epsilon < *tol){PRINTF("Warning: initial epsilon (%f) < tol(%f), resetting epsilon\n",*epsilon,*tol);*epsilon = 2.0 * (*tol);}
             
#ifdef MEASURE_ELAPSED_TIME
    t1_init = clock();
    t0_loop = t1_init;
#endif

#ifdef _DEBUG
FPRINTF(_F,"minQuad() main loop\n");
#endif

    int ncache = 0;
    int* cache = NULL;	
    while(!done) {
        (*iter)++;
        
        memcpy(a0,a,n1n2 * sizeof(double)); // save alpha


        // gradient with rounding
        for(p = 0; p < n1n2; p++){
            df[p] = Ha[p] + b[p];
            // if(ABS(df[p]) <= machine_double_eps){
                // df[p] = 0.0;
                // Ha[p] = -b[p];
                // Ha_bits[p] = 0.0;
            // }
        }


#ifdef _DEBUG
FPRINTF(_F,"Get violators and working set\n");
//PRINTF("Get violators and working set\n");
#endif         
    
        
        nV = getKKTViolators_0C(n1n2,a,df,C,&nV0,&nVC,v);        
        nB = MIN(dimB,nV); 			
        if(dynamic)nB = (int)ceil(pow((double)nV,.5)); 
				
        switch((WS_METHOD_TYPE)(*working_set_method)){
            // case RANDOM_V_WEIGHTED_RANK:
                //nB = getWorkingSet_rvwr(nB,B,nV0,nVC,n1n2,df,v,state.ibuf,state.ibuf2,state.ibuf3,state.dbuf);
                // break;
            case GREEDY: 
                nB = getWorkingSet_greedy_0C(nB,B,n1n2,a,df,C,v,state.ibuf,state.ibuf2,state.dbuf);                    
                break;
            case VL:
                nB = getWorkingSet_v(nB,B,nV0,nVC,n1n2,a,df,v,state.ibuf,state.ibuf2,state.dbuf);            
                break;            
            case VL2:
                nB = getWorkingSet_v2(nB,B,nV0,nVC,n1n2,a,df,v,state.ibuf,state.ibuf2,state.dbuf);
                break;
            case RANDOM_VL:
                nB = getWorkingSet_rv(nB,B,nV0,nVC,n1n2,v,state.ibuf,state.ibuf2);
                break;
            case RANDOM_VL2:
                nB =  getWorkingSet_rv2(nB,B,nV0,nVC,n1n2,v,state.ibuf,state.ibuf2);
                break;
            case RANDOM_VL2_WEIGHTED_GRAD:
                nB = getWorkingSet_rv2wg(nB,B,nV0,nVC,n1n2,df,v,state.ibuf,state.ibuf2,state.dbuf); 
                break;
            case RANDOM_VL_WEIGHTED_GRAD:
                nB  = getWorkingSet_rvwg(nB,B,nV0,nVC,n1n2,df,v,state.ibuf,state.ibuf2,state.dbuf);
                break;
            default: 
                 PRINTF("\nInvalid working set method.\n");
                 *convergence = 1;
                 done = 1;
                 break;
        }                
//        if(trace > 3)if(nB < dimB)PRINTF("Unable to select %i sized working set from %i violators, selected %i only.\n",dimB,nV,nB);
        
        if(*return_working_set) {          
            for(q1 = 0; q1 < nB; q1++) working_set[(*iter - 1)*dimB + q1] = B[q1] + 1; //'+1' for 1-based indexing
            for(q1 = nB; q1 < dimB; q1++) working_set[(*iter - 1)*dimB + q1] = -1;
        }
#ifdef _DEBUG
//PRINTF("only min(20,%i) shown:\n",nB);for(q1 = 0;q1 < MIN(20,nB);q1++)PRINTF("%i ",B[q1]);PRINTF("\n");
#endif 

////////////////////////////////////////////////////////////////////////////////////////////        
#ifdef _DEBUG
FPRINTF(_F,"Get epsilon\n");        
//PRINTF("Get epsilon\n");        
#endif  

        *epsilon = getEpsilonKKTV(n1n2,a,df,C);         
        if(*epsilon <= *tol) {
            done = 1;
            *convergence = 0;
            break;
        }
////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
FPRINTF(_F,"Optimization: \n");        
//PRINTF("Optimization: \n");        
#endif  

// UPDATE and CACHING
// if cache is NULL then all alphas are updated else cache points to working set. 
// ncache is either 1 (if working set collapsed) or 2 <= . 
// sweep - optimize (a subset of) all variables 1 variable at-a-time
// sweep - all alpha's and c's may change hence we must update at more alpha's than dimension of working set

        int success = 0;
	    int verbose = MAX(0,trace - 2);
        if(nB == 0) break;
        if((((*iter) % s) == 0)){
            ncache = n1n2;cache = NULL; 
            optimize_1(machine_double_eps,n1n2,cache,C,n1n2,n1,n2,type_H,H,b,c_bits,&c[0],a);//c,a updated   
            success = 1;
        }else{
            switch((QP_METHOD_TYPE)*optim){
                case EXHAUSTIVE2:
                    ncache = nB;cache = &B[0]; 
                    if(!get_QP(n1,n2,(int)type_H,H,&Ha[0],a,b,nB,&B[0],C,0,NULL,NULL,NULL,1,0,*_rank,&main_qp)) 
						break;                     
                    success = optimize_qp_exhaustive_fast(&main_qp,&sub_qp,&aB[0],0,verbose);
                    if(!success) break;                    
                    for(q1 = 0; q1 < nB; q1++)a[B[q1]] = aB[q1];
                    if(UPDATE_C)update_c_solveQuad(machine_double_eps,nB,&B[0],n1n2,n1,n2,type_H,H,a0,a,c_bits,c);                    
                    break;  
                case EXHAUSTIVE:         
                    ncache = nB;cache = &B[0];                    
                    // if(nB == 1){
                        // optimize_1(machine_double_eps,nB,&B[0],C,n1n2,n1,n2,type_H,H,b,c_bits,&c[0],a);//c,a updated   
                        // success = 1;
                        // break;
                    // }
                    
// single exhaustive 'run'                    
                    if(trace > 5) PRINTF("Get 'qp' structure:\n");
                    if(!get_QP(n1,n2,(int)type_H,H,&Ha[0],a,b,nB,&B[0],C,0,NULL,NULL,NULL,1,0,*_rank,&main_qp)) 
						break; 
                    if(trace > 5) PRINTF("Optimize Exhaustive\n");
                    success = optimize_qp_exhaustive(&main_qp,&sub_qp,&aB[0],0,verbose);
                    if(!success) break;                    
                
#ifdef CHECK_EXHAUSTIVE_SOLUTION                    
// detect mismatch w.r.t tron or exahustive_fast
					optimize_qp_tron(&main_qp.n,main_qp.H,main_qp.b,main_qp.x0,main_qp.xl,main_qp.xu,&aB_check[0],
									 optim_control,&verbose,&success); 
                    int mismatch = 0;                    
                    for(q1 = 0; q1 < nB; q1++)if(ABS(aB[q1] - aB_check[q1]) > .1){
                        mismatch = 1;
                        break;
                    };
                    if(mismatch){
                        *maxit = *iter;
                         optimize_qp_exhaustive(&main_qp,&sub_qp,&aB[0],0,10);
//                         optimize_qp_exhaustive_fast(&main_qp,&sub_qp,&aB_check[0],0,10);
                         
                        PRINTF("exhaustive f(x)=%+.5e x=",objective_QP(&main_qp,aB));for(q1 = 0; q1 < nB; q1++)PRINTF("%+.3e ",aB[q1]);PRINTF("\n"); 
                        PRINTF("check      f(x)=%+.5e x=",objective_QP(&main_qp,&aB_check[0]));for(q1 = 0; q1 < nB; q1++)PRINTF("%+.3e ",aB_check[q1]);PRINTF("\n"); 
                    }
#endif                        
                    for(q1 = 0; q1 < nB; q1++)a[B[q1]] = aB[q1];
                    if(UPDATE_C)update_c_solveQuad(machine_double_eps,nB,&B[0],n1n2,n1,n2,type_H,H,a0,a,c_bits,c);                    
                break;
                case LOQO:
                    ncache = nB;cache = &B[0];         
                    if(!get_QP(n1,n2,(int)type_H,H,&Ha[0],a,b,nB,&B[0],C,
							*n_constr,mat_constr,lhs_constr,rhs_constr,0,0,*_rank,&main_qp)) 
								break; 
                    if(trace > 5){PRINTF("\nloqo optimizer\n");print_QP(&main_qp);}                    
				
					success = optimize_qp_loqo(&main_qp, &aB[0],optim_control,verbose);               

					
                    //PRINTF("success %i ",success);for(q1 = 0; q1 < nB; q1++)PRINTF(" %f ",aB[q1]);PRINTF("\n");
                    if(!success) break;
                    for(q1 = 0; q1 < nB; q1++)a[B[q1]] = aB[q1];
                    if(UPDATE_C)update_c_solveQuad(machine_double_eps,nB,&B[0],n1n2,n1,n2,type_H,H,a0,a,c_bits,c);                                
                break;
                case TRON:
                    ncache = nB;cache = &B[0]; 
                    if(!get_QP(n1,n2,(int)type_H,H,&Ha[0],a,b,nB,&B[0],C,0,NULL,NULL,NULL,0,0,*_rank,&main_qp)) 
								break;
                    if(trace > 5){PRINTF("\ntron optimizer\n"); print_QP(&main_qp);}
					optimize_qp_tron(&main_qp.n,main_qp.H,main_qp.b,main_qp.x0,main_qp.xl,main_qp.xu,&aB[0],
									 optim_control,&verbose,&success); 
                    if(!success) break;
                    for(q1 = 0; q1 < nB; q1++)a[B[q1]] = aB[q1];
                    if(UPDATE_C)update_c_solveQuad(machine_double_eps,nB,&B[0],n1n2,n1,n2,type_H,H,a0,a,c_bits,c);    
                    break;
                default:        
                    PRINTF("Invalid optimizer selected\n");
                    break;                
            }
        }        
        
#ifdef _DEBUG
        FPRINTF(_F,"Optimization done.\n");
//        PRINTF("Optimization done.\n");
#endif       

        for(p = 0;p < nB;p++){if(ABS(a[B[p]] - C[B[p]]) <= machine_double_eps) a[B[p]] = C[B[p]];}
        for(p = 0;p < nB;p++)if(ABS(a[B[p]]) <= machine_double_eps) a[B[p]] = 0.0;      
        if(!success) break;
      
//////////////////////////////////////////////////////////////////////////////////

        
// at this point we can return 'final' alpha 'a' to R
        if(*iter == *maxit) {
            done = 1;
            *convergence = 1;
            if(trace > 1)trace = 1;
            break;
        }        
   
#ifdef _DEBUG
//PRINTF("Caching (%i) entries.",ncache);        
//if(cache)PRINTF("cache[0](%i) cache[1](%i),...",cache[0],cache[1]);PRINTF("\n");    
//PRINTF("da != 0:  ");for(p = 0; p < n1n2; p++)if(a[p]!=a0[p]) PRINTF("(p a[p] a0[p])) = (%i %+.3f %+.3f) ",p,a[p],a0[p]); PRINTF(" \n");        
FPRINTF(_F,"Caching H * alpha.\n");
//PRINTF("Caching H * alpha.\n");
#endif


    
// Computing Ha = H %*% a
// Use 'update' to reduce complexity and/or low-order representation Q = KK' 
// if recompute then probably best to use LAPACK 
        if(!*mem_efficient) {

//    Naive update of Ha - NOT PRECISE ENOUGH FOR SOME PROBLEMS    
            for(q1 = 0; q1 < ncache; q1++) {
                if(cache) q = cache[q1];else q = q1;
                double da = a[q] - a0[q];if(ABS(da) <= machine_double_eps) continue;
                for(p = 0; p < n1n2; p++){
#ifdef DO_KAHAN			
                    double y = Q[p][q] * da - Ha_bits[p];
                    double t = Ha[p] + y;
                    Ha_bits[p] = t - Ha[p] - y;
                    Ha[p] = t;                    
#else
                    Ha[p] += Q[p][q] * da; // this 'naive' sum will lose digits and alg. will diverge for some cases, e.g. data with outliers
#endif					
                }
            }    
            
        } else {
//            double (*K)[(unsigned int)(n1+n2)] = (double (*)[(unsigned int)(n1+n2)])_H;	
            double** K = Q;
            for(q1 = 0; q1 < ncache; q1++) {
                if(cache) q = cache[q1];else q = q1;
                double da = a[q] - a0[q];if(ABS(da) <= machine_double_eps) continue;
                iq=intCEIL(q+1,n2);//iq=(int)ceil((q+1.0)/n2);         
                jq=(q+1)-(iq-1)*n2;                     
                for(p = 0; p < n1n2; p++) {
                    ip=intCEIL(p+1,n2);//ip=(int)ceil((p+1.0)/n2);         
                    jp=(p+1)-(ip-1)*n2; 
#ifdef DO_KAHAN			                    
                    double y = -Ha_bits[p] + da * (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]);
                    double t = Ha[p] + y;
                    Ha_bits[p] = t - Ha[p] - y;
                    Ha[p] = t;
#else                    
                    Ha[p] += da * (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]);
#endif					
                }
            }
        }
 

// #ifdef _DEBUG
// PRINTF("Compute objective funcion\n");   
// #endif 
// compute objective value based on update
          *obj = 0.0;
		  if(trace > 1)
			for(q = 0; q < n1n2; q++)
				*obj += a[q] * (0.5*Ha[q] + b[q]);

//        *obj = objective(machine_double_eps,n1,n2,n1n2,type_H,H,a,b);
                
        if(trace > 1){
            PRINTF("m(%s) ws(%s) ",QP_METHOD_NAMES[*optim],WORKINGSET_METHOD_NAMES[*working_set_method]);
            PRINTF("i(%i) f(%+.5e) e(%+.5e) ", *iter,*obj,*epsilon);
            PRINTF("q(%i) nB(%i) nV(%i) = nV0(%i) + nVC(%i) ",dimB,nB,nV,nV0,nVC);
			PRINTF("B(");for(q = 0;q < nB;q++)PRINTF("%i ",B[q]+1);PRINTF(") ");
            //  int nb = MIN(nB,10);
            // PRINTF("B(");for(q = 0;q < nb;q++)PRINTF("%i ",B[q]+1);PRINTF(") ");
            // PRINTF("-df[B](");for(q = 0;q < nb;q++)PRINTF("%+.3e ",-df[B[q]]);PRINTF(") ");
            // PRINTF("a[B](");for(q = 0;q < nb;q++)PRINTF("%+.3e ",a[B[q]]);PRINTF(") ");
			// PRINTF("singular values ");
			// for(q = 0; q < nb/2; q++)PRINTF("%.3f ",main_qp.singular_values[q]);
			// for(q = 0; q < main_qp.n-nb/2; q++)PRINTF("%.4f ",main_qp.singular_values[q]);
            PRINTF("\n");
        }
           
    } // end while loop	
#ifdef _DEBUG
FPRINTF(_F,"minQuad main loop done\n");
#endif

	
#ifdef MEASURE_ELAPSED_TIME
    t1_loop = clock();
    t0_end = t1_loop;
#endif
    
// in case maxit = 1 and a sweep occured    
    *epsilon = getEpsilonKKTV(n1n2,a,df,C);  
//    *obj = objective(n1,n2,n1n2,type_H,H,a,b);
    *obj = 0.0;for(q = 0; q < n1n2; q++)*obj += a[q]*(0.5*Ha[q] + b[q]);

/* 
//  for debugging in R
    memcpy(__df,&df[0],n1n2 * sizeof(double));
    memcpy(__Ha,&Ha[0],n1n2 * sizeof(double));
    memcpy(__v,&v[0],n1n2 * sizeof(int));
*/    
    if(trace){
        PRINTF("m(%s) ws(%s) ",QP_METHOD_NAMES[*optim],WORKINGSET_METHOD_NAMES[*working_set_method]);
        PRINTF("i(%i) f(%+.5e) e(%+.5e) ", *iter,*obj,*epsilon);
        PRINTF("q(%i) nB(%i) nV(%i) = nV0(%i) + nVC(%i) ",dimB,nB,nV,nV0,nVC);
        PRINTF("B(");for(q = 0;q < nB;q++)PRINTF("%i ",B[q]+1);PRINTF(") ");
//        int nb = MIN(nB,10);
//        PRINTF("B(");for(q = 0;q < nb/2;q++)PRINTF("%i ",B[q]+1);for(q = nB - nb/2;q < nB;q++)PRINTF("%i ",B[q]+1);PRINTF(") ");
//        PRINTF("-df[B](");for(q = 0;q < nb/2;q++)PRINTF("%+.3e ",-df[B[q]]);for(q = nB - nb/2;q < nB;q++)PRINTF("%+.3e ",-df[B[q]]);PRINTF(") ");
//        PRINTF("a[B](");for(q = 0;q < nb/2;q++)PRINTF("%+.3e ",a[B[q]]);for(q = nB - nb/2;q < nB;q++)PRINTF("%+.3e ",a[B[q]]);PRINTF(") ");
//        PRINTF("singular values ");
//		for(q = 0; q < nb/2; q++)PRINTF("%.3f ",main_qp.singular_values[q]);
//		for(q = 0; q < main_qp.n-nb/2; q++)PRINTF("%.4f ",main_qp.singular_values[q]);
        PRINTF("\n");
    }

    free_minQuad();      

    
#ifndef SCYTHE_COMPILE_DIRECT    
    PutRNGstate();    
#endif

#ifdef MEASURE_ELAPSED_TIME
    t1_end = clock();
    PRINTF("minQuad elapsed initialization : %f seconds\n", (double)(t1_init - t0_init) / CLOCKS_PER_SEC);
    PRINTF("minQuad elapsed init cache     : %f seconds\n", (double)(t1_init_cache - t0_init_cache) / CLOCKS_PER_SEC);
    PRINTF("minQuad elapsed loop           : %f seconds\n", (double)(t1_loop - t0_loop) / CLOCKS_PER_SEC);
    PRINTF("minQuad elapsed exit           : %f seconds\n", (double)(t1_end - t0_end) / CLOCKS_PER_SEC);    
#endif

#ifdef _DEBUG
FPRINTF(_F,"minQuad() closing %s and returning\n",_DEBUG_FILENAME);
fclose(_F);
#endif
}

