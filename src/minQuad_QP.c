#include "minQuad_QP.h"
#include "minQuad.h"
#include "matrix.h"
#include <R_ext/Lapack.h>
#include <R.h>
#define PRINTF Rprintf
#include <float.h>
#define _Inf DBL_MAX
#define CEIL(X,Y) (1 + (((X) - 1) / (Y)))  
#define CX(i,j,ld) (((j)*(ld))+(i)) //indexing matrices in column-major order, 'ld' = 'leading dimension of matrix'

// ceil(x/y) where 'x' and 'y' are integers
// equivalent to  (x + y - 1) / y where '/' stands for integer division
// it also avoids overflow in y
// http://stackoverflow.com/questions/2745074/fast-ceiling-of-an-integer-division-in-c-c 


void free_QP(MNQD_QP* qp){
    if(qp->H) free(qp->H);
    if(qp->iH) free(qp->iH);
    if(qp->b)  free(qp->b);
    if(qp->x0) free(qp->x0);
    if(qp->xl) free(qp->xl);
    if(qp->xu) free(qp->xu);
    if(qp->dbuf) free(qp->dbuf);
    if(qp->ibuf) free(qp->ibuf);
    if(qp->singular_values) free(qp->singular_values);
	if(qp->n_constr){
		if(qp->mat_constr)free(qp->mat_constr);			
		if(qp->lhs_constr)free(qp->lhs_constr);			
		if(qp->rhs_constr)free(qp->rhs_constr);			

	}

    memset(qp,0,sizeof(MNQD_QP));
}
// double objective_QP(MNQD_QP* qp,double* x){
    
    // double (*H)[qp->n] = (double (*)[qp->n])qp->H;
    // double v = 0.0;
    // for(int p = 0; p < qp->n; p++) {
        // for(int q = 0; q < qp->n; q++) {
            // v   +=   H[p][q] * x[p] * x[q];
        // }
    // }
    // v *= 0.5;
    // for(int p = 0; p < qp->n; p++) v += x[p] * qp->b[p];    
    // return v;
// }

double objective_QP(MNQD_QP* qp,double* x){    
	// /* evaluate the function value f(x) = 0.5x'*H*x + b'*x */  
    int n = qp->n;double one = 1.0,zero = 0.0;int inc = 1;
	F77_CALL(dsymv)("U", &n, &one, qp->H, &n, x, &inc, &zero, qp->dbuf, &inc FCONE);
	return 0.5*F77_CALL(ddot)(&n, x, &inc, qp->dbuf, &inc) + F77_CALL(ddot)(&n, x, &inc, qp->b, &inc);
}

// return norm of gradient
// double objective_QP(MNQD_QP* qp,double* x){   
    // double norm = 0.0;
    // double one = 1.0;int inc = 1;
    // double* grad = qp->dbuf;
	// memcpy(grad,qp->b, sizeof(double)*qp->n);
	// F77_CALL(dsymv)("U", &(qp->n), &one, qp->H, &(qp->n), x, &inc, &one, grad, &inc);    
    // for(int p = 0; p < qp->n; p++) norm += grad[p] * grad[p];
    // return norm;
// }


void print_QP(MNQD_QP* qp){
    PRINTF("\nMinQuad QP-structure:\nQuadratic Minimization Problem f(x) = 0.5x'Hx + b'x subject to xl <= x <= xu\n");
    PRINTF("Current number variables (%i) and allocated space (%i)\n",(int)qp->n,(int)qp->N);
    PRINTF("initial values 'x0' ");for(int i = 0;i < qp->n;i++)PRINTF("%+.3e ",qp->x0[i]);PRINTF("\n");
    PRINTF("linear term 'b'     ");for(int i = 0;i < qp->n;i++)PRINTF("%+.3e ",qp->b[i]);PRINTF("\n");
    PRINTF("lower bounds 'xl'   ");for(int i = 0;i < qp->n;i++)PRINTF("%+.3e ",qp->xl[i]);PRINTF("\n");
    PRINTF("upper bounds 'xu'   ");for(int i = 0;i < qp->n;i++)PRINTF("%+.3e ",qp->xu[i]);PRINTF("\n");
    PRINTF("# of equality constraints");PRINTF(" %i\n",qp->n_constr);
    PRINTF("Matrix H\n");
//    double (*H)[qp->n] = (double (*)[qp->n])qp->H;
    for(int i = 0;i < qp->n;i++){
//        for(int j = 0;j < qp->n;j++)PRINTF("%+.3e ",H[i][j]);
        for(int j = 0;j < qp->n;j++)PRINTF("%+.3e ",qp->H[i*qp->n+j]);
        PRINTF("\n");
    } 
    if(qp->iH){
        PRINTF("Inverse Matrix of H\n");
//        double (*H)[qp->n] = (double (*)[qp->n])qp->iH;
        for(int i = 0;i < qp->n;i++){
//            for(int j = 0;j < qp->n;j++)PRINTF("%+.3e ",H[i][j]);
			for(int j = 0;j < qp->n;j++)PRINTF("%+.3e ",qp->iH[i*qp->n+j]);
            PRINTF("\n");
        } 
    }
    if(qp->singular_values){
        PRINTF("singular values     ");
        for(int i = 0;i < qp->n;i++)PRINTF("%+.3e ",qp->singular_values[i]);
        PRINTF("\n");
    }
}

// Q is symmetric 
//get and store submatrix Q[B,B] specified by index set B of length 'nB' in _sub_Q
void get_sub_Q(int n1,int n2,int n1n2,int type_H,void* H,int nB,int* B,double* sub_Q){

    DATA_PTR_TYPE type = type_H;
    double (*Q_BB)[nB] = (double (*)[nB])(sub_Q);
    if(type == SINGLE_Q) 
    {
//        double (*Q)[n1n2] = (double (*)[n1n2])(H);
        double* Q = (double*)H;
        for(int p = 0; p < nB; p++) {
            for(int q = 0; q < nB; q++) {
//                Q_BB[p][q] = Q[B[p]][B[q]];
				Q_BB[p][q] = Q[B[p] * n1n2 + B[q]];
				
            }
        }        
    }else if(type == DOUBLE_Q)
    {
        double** Q = (double**) H;
        for(int p = 0; p < nB; p++) {
            for(int q = 0; q < nB; q++) {
                Q_BB[p][q] = Q[B[p]][B[q]];
            }
        }   
    }else if(type == SINGLE_K){
//      double (* K)[n1+n2] = (double (*)[n1+n2])H;
        double* K = (double*)H;
		int n = n1 + n2;
        for(int q1 = 0; q1 < nB; q1++) {
            int p = B[q1];
            int ip=CEIL(p+1,n2);//int ip=(int)ceil((p+1.0)/n2); 
            int jp=(p+1)-(ip-1)*n2;
            for(int q2 = 0; q2 < nB; q2++) {
                int q = B[q2];
                int iq=CEIL(q+1,n2);//int iq=(int)ceil((q+1.0)/n2); 
                int jq=(q+1)-(iq-1)*n2;                     
//                Q_BB[q1][q2] =  K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1];
				Q_BB[q1][q2] = K[(ip-1)*n + iq-1]+K[(jp+n1-1)*n + jq+n1-1]-K[(ip-1)*n + jq+n1-1]-K[(jp+n1-1)*n + iq-1];
            }
        }
    }else if(type == DOUBLE_K){
        double** K = (double**)H;
        for(int q1 = 0; q1 < nB; q1++) {
            int p = B[q1];
            int ip=CEIL(p+1,n2);//int ip=(int)ceil((p+1.0)/n2); 
            int jp=(p+1)-(ip-1)*n2;
            for(int q2 = 0; q2 < nB; q2++) {
                int q = B[q2];
                int iq=CEIL(q+1,n2);//int iq=(int)ceil((q+1.0)/n2); 
                int jq=(q+1)-(iq-1)*n2;                     
                Q_BB[q1][q2] = (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]);
            }
        }
    }          
}


// output: 'sub' contains quadratic optimization sub-problem
// it is a 'sub' problem of the above sub-problem 
// 'restriced' further to a sub-combination of {B[1]..B[q]}
// its pointers should be freed at program termination
// n = Qp->n if ix == NULL else n = length(ix) 
int get_sub_QP(int n,int* ix,int do_inverse,int check_min_existence,MNQD_QP* Qp,MNQD_QP* sub){
    if(n < 1) return 0;
    MNQD_QP* qp = sub;   
    qp->n = n;
    if(qp->N < n){    
        qp->H  = (double*) realloc(qp->H, n * n * sizeof(double));if(!qp->H){qp->H = NULL;return 0;}
        qp->iH = (double*) realloc(qp->iH, n * n * sizeof(double));if(!qp->iH){qp->iH = NULL;return 0;}
        qp->b  = (double*) realloc(qp->b, n * sizeof(double));if(!qp->b){qp->b = NULL;return 0;}
        qp->x0 = (double*) realloc(qp->x0,n * sizeof(double));if(!qp->x0){qp->x0 = NULL;return 0;}
        qp->xl = (double*) realloc(qp->xl,n * sizeof(double));if(!qp->xl){qp->xl = NULL;return 0;}
        qp->xu = (double*) realloc(qp->xu,n * sizeof(double));if(!qp->xu){qp->xu = NULL;return 0;}
		if(qp->n_constr){
			qp->mat_constr = (double*) realloc(qp->mat_constr,qp->n_constr * n * sizeof(double));if(!qp->mat_constr){qp->mat_constr = NULL;return 0;}			
			qp->lhs_constr = (double*) realloc(qp->lhs_constr,qp->n_constr * sizeof(double));if(!qp->lhs_constr){qp->lhs_constr = NULL;return 0;}
			qp->rhs_constr= (double*) realloc(qp->rhs_constr,qp->n_constr * sizeof(double));if(!qp->rhs_constr){qp->rhs_constr = NULL;return 0;}
        }
		qp->dbuf = (double*) realloc(qp->dbuf,n * sizeof(double));if(!qp->dbuf){qp->dbuf = NULL;return 0;}
        qp->ibuf = (int*) realloc(qp->ibuf,n * sizeof(int));if(!qp->ibuf){qp->ibuf = NULL;return 0;}
        qp->N  = n;
        qp->singular_values = (double*) realloc(qp->singular_values, n * sizeof(double));
        if(!qp->singular_values){qp->singular_values = NULL;return 0;}        
    }
    qp->machine_double_eps = Qp->machine_double_eps;
    // double (*H)[Qp->n] = (double (*)[Qp->n])Qp->H;
    // double (*h)[qp->n] = (double (*)[qp->n])qp->H;
    // for(int i = 0;i < n;i++){
        // int p = i;if(ix) p = ix[i];
        // for(int j = 0;j < n;j++){
            // int q = j;if(ix) q = ix[j];
            // h[i][j] = H[p][q];
        // }
    // }

    for(int i = 0;i < n;i++){
        int p = i;if(ix) p = ix[i];
        for(int j = 0;j < n;j++){
            int q = j;if(ix) q = ix[j];
			qp->H[i*qp->n + j] = Qp->H[p*Qp->n + q];
        }
    }	
   
    for(int i = 0;i < n;i++){
        int ki = i;if(ix) ki = ix[i];
        qp->xl[i] = Qp->xl[ki]; 
        qp->xu[i] = Qp->xu[ki]; 
        qp->x0[i] = Qp->x0[ki];
// lin term in restricted obj. function is computed directly here:        
        qp->b[i] = Qp->b[ki]; 
        for(int j = 0; j < Qp->n - n; j++){
            int kj = j;if(ix){if(j == ix[j])continue;}
            qp->b[i] += Qp->H[ki*Qp->n + kj] * Qp->x0[kj];
//            qp->b[i] += H[ki][kj] * Qp->x0[kj];
        }
    }
	
// constraints matrix is in COLUMN-major order and not symmetric in general
	if(qp->n_constr){
		for(int j = 0;j < n;j++){
			int kj = j;if(ix) kj = ix[j];
			qp->lhs_constr[j] = Qp->lhs_constr[kj];
			qp->rhs_constr[j] = Qp->rhs_constr[kj];
			// for(int i = 0;i < qp->n_constr;i++){
				// qp->mat_constr[CX(i,j,qp->n_constr)] = Qp->mat_constr[CX(i,kj,Qp->n_constr)];
			// }
		}			
		get_sub_matrix(NULL,ix,"c",Qp->n_constr,Qp->n,Qp->mat_constr,"c",qp->n_constr,qp->n,qp->mat_constr);		
	}
	

    for(int i = 0;i < n;i++)qp->singular_values[i] = -1.0;     
    if(do_inverse){
        if(luinv(qp->n,qp->ibuf,qp->dbuf,qp->H,qp->iH))
            ginv(sqrt(qp->machine_double_eps),qp->n,qp->n,qp->H,qp->iH,qp->singular_values);        
    }

    if(check_min_existence){
// 0.5x'Hx + b'x has min iff H >= 0 (p.s.d) and (I-HH*)b = 0 for H* = gen. inv(H)
// (I-HH*)b <=> b - (HH*)b == 0 <=> b - H(Hinv b) == 0
//F77_CALL(dsymm)("n", "n", &N, &M, &M, &ALPHA, qp->H, &N, &u[0],&M, &BETA,qp->iH,&N);
            double b[qp->n];
            for(int j = 0;j < qp->n;j++)b[j] = qp->b[j];  
            double buf[qp->n];
            double zero = 0.0,one = 1.0;int inc = 1;
            F77_CALL(dsymv)("U", &n, &one, qp->iH, &n, &b[0], &inc, &zero, &buf[0], &inc FCONE);
            F77_CALL(dsymv)("U", &n, &one, qp->H, &n, &buf[0], &inc, &zero, &b[0], &inc FCONE);
            zero = 0.0; for(int j = 0;j < qp->n;j++)zero += (b[j] - qp->b[j]) * (b[j] - qp->b[j]);
            if(zero > .1){
//                print_QP(qp);
//                PRINTF("Warning: minimum of sub-sub-problem may not be found by Moore-Penrose inverse, 0 != ||(I - HH*)b|| = %f for problem\n",zero);
            }
    }
    return 1;            
}

// output: 'qp' contains quadratic optimization problem
// it is a 'sub' problem as it is based on working set B
// it can be passed to external optimizers
// its pointers should be freed at program termination
// working set 'B' subsets Q/K,a,b 

// Q[B,-B] %*% a[-B] = Q[B,B]%*%a[B] - Ka[B]
// bB = b[B] + Q[B,-B] %*% a[-B]

int get_QP(int n1,int n2,int type_H,void* H,double* Ha,double* a,double* b,int nB,int* B,double* C,
int n_constr,double* mat_constr,double* lhs_constr,double* rhs_constr,
int do_inverse,int check_min_existence,int rank,MNQD_QP* qp){
    int i,n = nB;
    if(n < 1) return 0;
    
    // assume all dynamic pointers in qp were initialized to NULL!
    // in minQuad() this should execute one time only
    if(qp->N < n){    
        qp->H  = (double*) realloc(qp->H, n * n * sizeof(double));if(!qp->H){qp->H = NULL;return 0;}
        qp->iH = (double*) realloc(qp->iH,n * n * sizeof(double));if(!qp->iH){qp->iH = NULL;return 0;}
        qp->b  = (double*) realloc(qp->b, n * sizeof(double));if(!qp->b){qp->b = NULL;return 0;}
        qp->x0 = (double*) realloc(qp->x0,n * sizeof(double));if(!qp->x0){qp->x0 = NULL;return 0;}
        qp->xl = (double*) realloc(qp->xl,n * sizeof(double));if(!qp->xl){qp->xl = NULL;return 0;}
        qp->xu = (double*) realloc(qp->xu,n * sizeof(double));if(!qp->xu){qp->xu = NULL;return 0;}

		if(n_constr){
			qp->n_constr = n_constr; 
			qp->mat_constr = (double*) realloc(qp->mat_constr,n_constr * n * sizeof(double));if(!qp->mat_constr){qp->mat_constr = NULL;return 0;}			
			qp->lhs_constr = (double*) realloc(qp->lhs_constr,n_constr * sizeof(double));if(!qp->lhs_constr){qp->lhs_constr = NULL;return 0;}
			qp->rhs_constr = (double*) realloc(qp->rhs_constr,n_constr * sizeof(double));if(!qp->rhs_constr){qp->rhs_constr = NULL;return 0;}
		}
		
        qp->dbuf = (double*) realloc(qp->dbuf,n * sizeof(double));if(!qp->dbuf){qp->dbuf = NULL;return 0;}
        qp->ibuf = (int*) realloc(qp->ibuf,n * sizeof(int));if(!qp->ibuf){qp->ibuf = NULL;return 0;}
        qp->N  = n;
        qp->singular_values = (double*) realloc(qp->singular_values, n * sizeof(double));
        if(!qp->singular_values){qp->singular_values = NULL;return 0;}        
    }
    // qp->_x = a;
    // qp->_H = (double*)H;
    // qp->_b = b;
    // qp->_Hx = Ka;
    // qp->_n = n1*n2;
    // qp->working_set = B;
    qp->n = n;   

	qp->n_constr = n_constr;
    for(i = 0;i < n;i++){
        qp->xl[i] = 0.0; 
        qp->xu[i] = C[B[i]];
        qp->x0[i] = a[B[i]];
		if(qp->n_constr){
			qp->lhs_constr[i] = lhs_constr[B[i]];
			qp->rhs_constr[i] = rhs_constr[B[i]];
		}
    }
// constraints matrix is in COLUMN-major order
	if(n_constr){
		for(int j = 0;j < n;j++){
			qp->lhs_constr[j] = lhs_constr[B[j]];
			qp->rhs_constr[j] = rhs_constr[B[j]];
			// for(int i = 0;i < qp->n_constr;i++){
				// qp->mat_constr[CX(i,j,qp->n_constr)] = mat_constr[CX(i,B[j],qp->n_constr)];
			// }
		}	
		get_sub_matrix(NULL,B,"c",n_constr,n1*n2,mat_constr,"c",qp->n_constr,qp->n,qp->mat_constr);				
	}	
	
    get_sub_Q(n1,n2,n1*n2,type_H,H,nB,B,qp->H);

    // linear term in restricted objective function
//    double (*H)[n] = (double (*)[n])qp->H;
    for(i = 0; i < n; i++) {
        qp->b[i] = b[B[i]] + Ha[B[i]];
        for(int j = 0; j < n; j++) qp->b[i] -= qp->H[i*n + j] * a[B[j]];
//        for(int j = 0; j < n; j++) qp->b[i] -= H[i][j] * a[B[j]];
    }   
    qp->rank = rank;
    for(int i = 0;i < n;i++)qp->singular_values[i] = -1.0; 
    if(do_inverse){
        if(qp->rank < 1){
            qp->rank = 0;
			double ginv_tol = sqrt(qp->machine_double_eps);
			ginv(ginv_tol,qp->n,qp->n,qp->H,qp->iH,qp->singular_values); // to get an approximation for rank(H)
            for(int i = 0;i < n;i++)if(qp->singular_values[i] > ginv_tol)qp->rank++; 
        }else {
			if(luinv(qp->n,qp->ibuf,qp->dbuf,qp->H,qp->iH))
				ginv(sqrt(qp->machine_double_eps),qp->n,qp->n,qp->H,qp->iH,qp->singular_values); 
		}
    }
    qp->max_order = qp->rank;
    if(check_min_existence){
// 0.5x'Hx + b'x has min iff H >= 0 (p.s.d) and (I-HH*)b = 0 for H* = gen. inv(H)
// (I-HH*)b <=> b - (HH*)b == 0 <=> b - H(Hinv b) == 0
//F77_CALL(dsymm)("n", "n", &N, &M, &M, &ALPHA, qp->H, &N, &u[0],&M, &BETA,qp->iH,&N);
            double b[qp->n];
            for(int j = 0;j < qp->n;j++)b[j] = qp->b[j];  
            double zero = 0.0,one = 1.0;int inc = 1;
            F77_CALL(dsymv)("U", &n, &one, qp->iH, &n, &b[0], &inc, &zero, qp->dbuf, &inc FCONE);//dbuf <- (H*)b 
            F77_CALL(dsymv)("U", &n, &one, qp->H, &n, qp->dbuf, &inc, &zero, &b[0], &inc FCONE);//b <- H(H*)b
            zero = 0.0; for(int j = 0;j < qp->n;j++)zero += (b[j] - qp->b[j]) * (b[j] - qp->b[j]);
            if(zero > .1){
//                print_QP(qp);
//                PRINTF("Warning: minimum of sub-problem may not be found by Moore-Penrose inverse, 0 != ||(I - HH*)b|| = %f for problem\n",zero);
            }else qp->max_order = n;
    }
    return 1;
}


int init_QP(int n,int n_constr,double machine_double_eps,MNQD_QP* qp){
    if(n < 1) return 0;
    memset(qp,0,sizeof(MNQD_QP));
    if(qp->N < n){    
        qp->H  = (double*) realloc(qp->H, n * n * sizeof(double));if(!qp->H){qp->H = NULL;return 0;}
        qp->iH = (double*) realloc(qp->iH, n * n * sizeof(double));if(!qp->iH){qp->iH = NULL;return 0;}
        qp->b  = (double*) realloc(qp->b, n * sizeof(double));if(!qp->b){qp->b = NULL;return 0;}
        qp->x0 = (double*) realloc(qp->x0,n * sizeof(double));if(!qp->x0){qp->x0 = NULL;return 0;}
        qp->xl = (double*) realloc(qp->xl,n * sizeof(double));if(!qp->xl){qp->xl = NULL;return 0;}
        qp->xu = (double*) realloc(qp->xu,n * sizeof(double));if(!qp->xu){qp->xu = NULL;return 0;}
		if(n_constr){
			qp->n_constr = n_constr;
			qp->mat_constr = (double*) realloc(qp->mat_constr,n * n_constr * sizeof(double));if(!qp->mat_constr){qp->mat_constr = NULL;return 0;}
			qp->lhs_constr = (double*) realloc(qp->lhs_constr,n_constr * sizeof(double));if(!qp->lhs_constr){qp->lhs_constr = NULL;return 0;}
			qp->rhs_constr = (double*) realloc(qp->rhs_constr,n_constr * sizeof(double));if(!qp->rhs_constr){qp->rhs_constr = NULL;return 0;}
		}
        qp->dbuf = (double*) realloc(qp->dbuf,n * sizeof(double));if(!qp->dbuf){qp->dbuf = NULL;return 0;}
        qp->ibuf = (int*) realloc(qp->ibuf,n * sizeof(int));if(!qp->ibuf){qp->ibuf = NULL;return 0;}
        qp->N  = n;
        qp->singular_values = (double*) realloc(qp->singular_values, n * sizeof(double));
        if(!qp->singular_values){qp->singular_values = NULL;return 0;}        
    }
    qp->n = n;
    for(int i = 0;i < n;i++)qp->singular_values[i] = -1.0;     
    for(int i = 0;i < n;i++){
        qp->xl[i] = 0.0; 
        qp->xu[i] = 1.0; 
        qp->x0[i] = 0.5;
    }
	qp->machine_double_eps = machine_double_eps;
    return 1;            
}
