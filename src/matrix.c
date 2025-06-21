// Author: Krisztian Sebestyen ksebestyen@gmail.com

#include "matrix.h"
#include <float.h> //DBL_EPSILON
#include <R_ext/Lapack.h>
#include <R.h>
#include <Rinternals.h>
#ifndef FCONE
# define FCONE
#endif
//#include <R_ext/Applic.h>
//#include <R_ext/BLAS.h>
//#include <R_ext/RS.h> //R definitions for 'extending' R, registering functions,...


#define PRINTF Rprintf
#define	MAX(A,B)	((A) > (B) ? (A) : (B))
#define	MIN(A,B)	((A) < (B) ? (A) : (B))


// LEADING DIMENSION 'LDA' or 'LD', comments collected from the 'net'
// 'ld' : 'leading dimension of matrix' is equal to the number of elements in its major dimension.
//      : the distance between two neighboring elements along minor dimension.
// It is used to access a sub-matrix of a larger matrix
// The leading dimension for a two-dimensional array is an increment that is used to find the 
// starting point for the matrix elements in each successive column of the array. 
// The leading dimension must always be positive. It must always be greater than or equal to the
// minor dimension of the matrix.
#define CX(i,j,ld) (((j)*(ld))+(i)) //indexing matrices in column-major order
#define RX(i,j,ld) (((i)*(ld))+(j)) //indexing matrices in row-major order
typedef enum {COLUMN_MAJOR = 0,ROW_MAJOR = 1}MAJOR_ORDER;




// 'Increment': increment of a vector is the number of storage elements from one element to the next.
// increment = 1 <=> vector is stored contiguously in memory. It can be used to :
// - manipulate rows of a matrix stored in column-major order  
// - manipulate columns of a matrix stored in row-major order  

// input: LAPACK's ipiv
// output: row-permutation 'p'
// note: LAPACK's ipiv definition is: "row i of the matrix was interchanged with row ipiv(i)"
// hence ipiv is NOT a permutation in the regular sense
// ipiv represents a permutation in a 'sequenced' form and this functions 'unsequences' it
void rowperm_ipiv(int* _n,int* ipiv,int* perm){
	int n = *_n;
	int i,tmp;
	 
	for(i = 0;i < n;i++) perm[i] = i + 1;
    for(i = 0;i < n;i++){
        tmp = perm[ipiv[i] - 1];
        perm[ipiv[i] - 1] = perm[i];
        perm[i] = tmp;
    }
}

// generate left permutation matrix corresponding to 'rowperm'
void row_PL(int* _n, int* rowperm,int* iP,double* dP){
	int n = *_n;
	int* perm = rowperm;
	int i;
    if(iP){
		for (i = 0; i < n; i++) 
			iP[i*n + perm[i] - 1] = 1;
    }else if(dP){
		for (i = 0; i < n; i++) 
			dP[i*n + perm[i] - 1] = 1.0;	
	}
}
void invperm(int* n,int* p,int* ip){ for(int i = 0;i < *n;i++)ip[p[i]-1] = i + 1;}

// row_major: offset = i*ncx + j
// col_major: offset = i + j*nrx
// // x_{nrx x ncx}
// // L_{nrx x ncx}
// // extract L_{nrx x k} lower triangular matrix
// // diag: value the diagonal should be set to 
// void lower_tri(int nrx,int ncx,double* x,int unit,int k,double* L){
	// memset(L,0,(size_t)(nrx * ncx * sizeof(double)));
	// for(int i = k;i < nrx;i++){
		// if(i < 0) continue;
		// for(int j = 0;j < i;j++){
			// L[i + j * nrx] = x[i + j * nrx];//column major order
		// }
	// }
	// if(unit)
		// for(int i = MAX(0,k-1);i < nrx;i++) 
			// L[i + i * nrx] = 1.0;
// }

// // x_{nrx x ncx}
// // U_{nrx x ncx}
// // extract U_{k x ncx} upper triangular matrix
// // diag: value the diagonal should be set to 
// void upper_tri(int nrx,int ncx,double* x,int unit,int k,double* U){
	// memset(U,0,(size_t)(nrx * ncx * sizeof(double)));
	// for(int j = k;j < ncx;j++){
		// if(j < 0) continue;
		// for(int i = 0;i < j;i++){
			// U[i + j * nrx] = x[i + j * nrx];//column major order
		// }
	// }
	// if(unit)
		// for(int j = MAX(0,k-1);j < ncx;j++) 
			// U[j + j * nrx] = 1.0;	
// }
// // R interface
// void tri_upper(int* nrx,int* ncx,double* x,int* unit,int* k,double* U){upper_tri(*nrx,*ncx,x,*unit,*k,U);}
// void tri_lower(int* nrx,int* ncx,double* x,int* unit,int* k,double* L){lower_tri(*nrx,*ncx,x,*unit,*k,L);}




// input:
// x_{nrx x ncx} matrix in column-major order
// k the position of the 'diagonal' of the trapezoid 
// diag: values the diagonal should be set to, max length min(nrx,ncx). if NULL, diagonal = diag(x).
// output:
// L_{nrx x ncx} in column-major order

// put the matrix in the 4th quadrant with (0,0) entry at the origin
// k = {-(K - 1), .. ,0, .. ,(K - 1)}, K = MAX(nrx,ncx)
// if k = 0 we are at the longest diagonal and extract the lower triangular part (including diagonal)
// if k > 0 the diagonal lines move up
// if k < 0 the diagonal lines move down
void lower_trap(int nrx,int ncx,double* x,double* diag,int k,double* L){
	int i,j,d;
	int K = MAX(nrx,ncx);
	memset(L,0,(size_t)(nrx * ncx * sizeof(double)));

	k = -k;//k > 0 <=> the diagonal lines move down
	
	// diagonal
	d = k;
	int dx = 0;
	for(j = 0;j < ncx;j++){
		i = j + d;
		if((i < 0) || (nrx <= i)) continue;
		L[i + j * nrx] = (diag == NULL) ? x[dx++] : diag[dx++];//column major order
		
		// Rprintf("L(%i,%i)=%f d=%f\n",i,j,L[i + j * nrx],diag);
	}
	// Rprintf("length of diagonal %i\n",dx);
	
	for(d = k + 1;d < K;d++){ // length of 'current' diagonal
		for(j = 0;j < ncx;j++){ // move along diagonal 
			i = j + d;// 'y-intercept' of diagonal
			if((i < 0) || (nrx <= i)) continue;
			L[i + j * nrx] = x[i + j * nrx];//column major order
		}
	}
}

// input:
// x_{nrx x ncx} matrix in column-major order
// k the position of the 'diagonal' of the trapezoid 
// diag: values the diagonal should be set to, max length min(nrx,ncx). if NULL, diagonal  = diag(x)
// output:
// U_{nrx x ncx} in column-major order

// put the matrix in the 4th quadrant with (0,0) entry at the origin
// k = {-(K - 1), .. ,0, .. ,(K - 1)}, K = MAX(nrx,ncx)
// if k = 0 we are at the longest diagonal and extract the lower triangular part (including diagonal)
// if k > 0 the diagonal lines move down  
// if k < 0 the diagonal lines move up
void upper_trap(int nrx,int ncx,double* x,double* diag,int k,double* U){
	int i,j,d;
	int K = MAX(nrx,ncx);
	memset(U,0,(size_t)(nrx * ncx * sizeof(double)));

	// set diagonal
	d = k;
	int dx = 0;
	for(j = 0;j < ncx;j++){
		i = j - d;
		if((i < 0) || (nrx <= i)) continue;
		U[i + j * nrx] = (diag == NULL) ? x[dx++] : diag[dx++];//column major order
	}
	for(d = k + 1;d < K;d++){ // length of 'current' diagonal
		for(j = 0;j < ncx;j++){ // move along diagonal 
			i = j - d;// 'y-intercept' of diagonal
			if((i < 0) || (nrx <= i)) continue;
			U[i + j * nrx] = x[i + j * nrx];//column major order
		}
	}
}



/////////////////////////// rep-rbind() and rep-cbind() //////////////////////////

// x <- matrix(c(11,12,21,22,31,32),nrow = 3,ncol = 2,byrow = TRUE);x
     // [,1] [,2]
// [1,]   11   12
// [2,]   21   22
// [3,]   31   32
// c(x)
// [1] 11 21 31 12 22 32

// c(reprbind(x,times = 2))
// 11 21 31 11 21 31 12 22 32 12 22 32

// c(reprbind(x,each = 2))
// 11 11 21 21 31 31 12 12 22 22 32 32

// 0->0,1
// 1->2,3
// 2->4,5
// i -> i*each,i*each+1

// x = _{nrx,ncx} in column-major order
// y = rep-cbind(x)
// exactly one of times,each > 0, the other is 0
// times is integer >= 0
// each is integer >= 0
// each takes precedence
// vec_each is length nrx, takes precedence over each 
// y = rep-rbind(x) 
void rrbind(double* _x,int nrx,int ncx,int times,int _each,int* vec_each,double* _y){

	double* x = _x;
	double* y = _y;
	
	if((_each > 0) || vec_each){
		for(int j = 0;j < ncx;j++){
			for(int i = 0;i < nrx;i++){
				int each = vec_each ? vec_each[i] : _each;
				int ix = CX(i,j,nrx);
				for(int k = 0;k < each;k++){
					*y = x[ix]; // with nesting order i(j(k)) we fill out y in column order automatically 
					y++;
				}
			}
		}
		return;
	}
	
	if(times > 0){
		for(int j = 0;j < ncx;j++){
			for(int k = 0;k < times;k++){
				memcpy(y,x,nrx*sizeof(double));
				y+=nrx;
			}
			x += nrx;
		}
	}
	return;		
}

// example for 'x' above 
// c(ccbind(x,times = 2))
// 11 21 31 12 22 32 11 21 31 12 22 32
// c(ccbind(x,each = 2))
// 11 21 31 11 21 31 12 22 32 12 22 32

// x = _{nrx,ncx} in column-major order
// y = rep-cbind(x)
// exactly one of times,each > 0, the other is 0
// times is integer >= 0
// each is integer >= 0
// each takes precedence
// vec_each is length ncx takes precedence over each
void rcbind(double* _x,int nrx,int ncx,int times,int each,int* vec_each,double* _y){

	double* x = _x;
	double* y = _y;
		
	if(vec_each){
		for(int j = 0;j < ncx;j++){
			int each_j = vec_each[j];
			for(int i = 0;i < nrx;i++){
				for(int k = 0;k < each_j;k++){
					*y = x[CX(i,j,nrx)];
					y++;
				}
			}
		}
		return;
	}else if(each > 0){
		for(int j = 0;j < ncx;j++){
			for(int k = 0;k < each;k++){
				memcpy(y,x,nrx*sizeof(double));
				y+=nrx;
			}
			x += nrx;
		}
		return;
	}
	if(times > 0){
		int size_x = nrx * ncx;
		for(int k = 0;k < times;k++){
			memcpy(y,x,size_x*sizeof(double));
			y+=size_x;
		}
	}
	return;
}


///////////////////////////////////////////////////////////////////////////////////

// functions crossprod() and tcrossprod() are directly from R

// x,y are passed in column-major order from R 
//x %*% y
void matprod(double *x, int* nrx, int* ncx,
		      double *y, int* nry, int* ncy, double *z)
{
    char *transa = "N", *transb = "N";
    double one = 1.0, zero = 0.0;
    if (*nrx > 0 && *ncx > 0 && *nry > 0 && *ncy > 0) {
	F77_CALL(dgemm)(transa, transb, nrx, ncy, ncx, &one,
			x, nrx, y, nry, &zero, z, nrx FCONE FCONE);
    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < (*nrx)*(*ncy); i++) z[i] = 0;
    }
}

//t(x) %*% y
void crossprod(double *x, int* nrx, int* ncx,
		      double *y, int* nry, int* ncy, double *z)
{
    char *transa = "T", *transb = "N";
    double one = 1.0, zero = 0.0;
    if (*nrx > 0 && *ncx > 0 && *nry > 0 && *ncy > 0) {
	F77_CALL(dgemm)(transa, transb, ncx, ncy, nrx, &one,
			x, nrx, y, nry, &zero, z, ncx FCONE FCONE);
    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < (*ncx)*(*ncy); i++) z[i] = 0;
    }
}

//x %*% t(y)
void tcrossprod(double *x, int* nrx, int* ncx,
		      double *y, int* nry, int* ncy, double *z)
{
    char *transa = "N", *transb = "T";
    double one = 1.0, zero = 0.0;
    if (*nrx > 0 && *ncx > 0 && *nry > 0 && *ncy > 0) {
	F77_CALL(dgemm)(transa, transb, nrx, nry, ncx, &one,
			x, nrx, y, nry, &zero, z, nrx FCONE FCONE);
    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < (*nrx)*(*nry); i++) z[i] = 0;
    }
}

// matrix multiplication D %*% X %*% D2 for D,D2 diagonal matrices and X square matrix in column-major order
// x,y are square matrices in column order, d1 and d2 represent the diagonals of diagonal matrices
void dxd_(int* _n, double* d1,double* x,double* d2,double* y){
	int i,j;
	int n = *_n;
	for(j = 0;j < n;j++)
		for(i = 0;i < n;i++)
			y[j*n + i] = d1[i] * x[j*n + i] * d2[j];
}

void dxd(int* _n, double* d1,double* x,double* d2,double* _ans){
	int n = *_n;
	double* ans = _ans;

	int i,j;
	for(j = 0;j < n;j++){
		for(i = 0;i < n;i++){
			*ans = d1[i] * (*x) * d2[j];
			ans++; x++;
		}
	}
}



// z = c1*x op c2*y, op = {+,-,*,/}<=>{0,1,2,3} 
void vec_op(double* z,double c1,double* x,int op,double c2,double* y,int n){
	int i;
	switch((op_type)op){
		case ADD: 
			for(i = 0;i < n;i++)z[i] = c1*x[i] + c2*y[i];
			break;
		case SUB: 
			for(i = 0;i < n;i++)z[i] = c1*x[i] - c2*y[i];
			break;
		case MULT: 
			for(i = 0;i < n;i++)z[i] = c1*x[i] * c2*y[i];
			break;
		case DIV: 
			for(i = 0;i < n;i++)z[i] = c1*x[i] / (c2*y[i]);
			break;			
	}
}

// x in column major order
void print_matrix(int m,int n,double* x){
    if(!x) return;
//    double (*x)[n] = (double (*)[n])_x;
    for(int i = 0;i < m;i++){
        for(int j = 0;j < n;j++)
			Rprintf("%+.4e ",x[CX(i,j,m)]);
        Rprintf("\n");
    } 
}


void get_sub_matrix(int* ix,int* jx,const char* _major_x,int nrx,int ncx,double* x,const char* _major_y,int nry,int ncy,double* y)
{
	int row_major_x = ((int)*_major_x == (int)'r') || ((int)*_major_x == (int)'R');	
	int row_major_y = ((int)*_major_y == (int)'r') || ((int)*_major_y == (int)'R');	
		
	if(!row_major_x && !row_major_y){  
		for(int i = 0;i < nry;i++){
			int p = (ix ? ix[i] : i);
			for(int j = 0;j < ncy;j++){
				int q = (jx ? jx[j] : j);
				y[CX(i,j,nry)] = x[CX(p,q,nrx)]; 
			}
		}
		return;
	}
	
	if(!row_major_x && row_major_y){  
		for(int i = 0;i < nry;i++){
			int p = (ix ? ix[i] : i);
			for(int j = 0;j < ncy;j++){
				int q = (jx ? jx[j] : j);
				y[RX(i,j,ncy)] = x[CX(p,q,nrx)]; 
			}
		}	
		return;
	}
	
	if(row_major_x && row_major_y){  
		for(int i = 0;i < nry;i++){
			int p = (ix ? ix[i] : i);
			for(int j = 0;j < ncy;j++){
				int q = (jx ? jx[j] : j);
				y[RX(i,j,ncy)] = x[RX(p,q,ncx)]; 
			}
		}	
		return;
	}

	if(row_major_x && !row_major_y){  
		for(int i = 0;i < nry;i++){
			int p = (ix ? ix[i] : i);
			for(int j = 0;j < ncy;j++){
				int q = (jx ? jx[j] : j);
				y[CX(i,j,nry)] = x[RX(p,q,ncx)]; 
			}
		}	
		return;
	}	
}
void R_get_sub_matrix(int* major_x,int* nrx,int* ncx,double* x,int* major_y,int* nry,int* ix,int* ncy,int* jx,double* y){		
	const char* MAJOR_X = (*major_x == ROW_MAJOR) ? "r" : "c";	
	const char* MAJOR_Y = (*major_y == ROW_MAJOR) ? "r" : "c";	
	get_sub_matrix(ix,jx,MAJOR_X,*nrx,*ncx,x,MAJOR_Y,*nry,*ncy,y);	
}



// obsolete
// // 'x' in column-major order
// // sub = x[ix,jx] assuming row-major order but for us x is symmetric
// // (ix,jx) = (NULL,NULL) implies no restriction
// void get_sub_matrix(int nrow,int ncol,double* x,int nix,int* ix,int njx,int* jx,double* sub)
// {    
// // if no restriction by ix or jx and 0 was passed for nix or njx then reset     
    // nix = nix ? nix : nrow;
    // njx = njx ? njx : ncol;
    // for(int i = 0;i < nix;i++){
        // int p = (ix ? ix[i] : i);
        // for(int j = 0;j < njx;j++){
            // int q = (jx ? jx[j] : j);
            // sub[i*nix + j] = x[p*nrow + q]; 
    // }}
// }

// obsolete
// void get_sub_matrix(int nrow,int ncol,double* _x,int nix,int* ix,int njx,int* jx,double* _sub)
// {    
// // if no restriction by ix or jx and 0 was passed for nix or njx then reset     
    // nix = nix ? nix : nrow;
    // njx = njx ? njx : ncol;
    
    // double (*sub)[njx] = (double (*)[njx])(_sub);
    // double (*x)[ncol] = (double (*)[ncol])(_x);
    
    
    // for(int i = 0;i < nix;i++){
        // int p = (ix ? ix[i] : i);
        // for(int j = 0;j < njx;j++){
            // int q = (jx ? jx[j] : j);
            // sub[i][j] = x[p][q]; 
    // }}
// }



// invert matrices by LU-factorization
// ipiv,work,x,are length-n vectors
// x,xinv are length n*n vectors
// return 0(success) or error code, side effect: 'x' is overwritten  
int luinv(int n,int* ipiv,double* work,double* x,double* xinv){
    int info = 0;
    memcpy(xinv,x,n * n * sizeof(double)); //dgetrf changes its argument
    F77_CALL(dgetrf)(&n,&n,xinv,&n,ipiv,&info); 
    if(info)return info;
    F77_CALL(dgetri)(&n,xinv,&n,ipiv,work,&n,&info);    
    return info;
}

// invert matrices by LU-factorization
// ipiv,work,x,are length-n vectors
// x is overwritten by xinv, length n*n vectors
// return 0(success) or error code, 'side effect': 'x' is overwritten  
int luinv0(int n,int* ipiv,double* work,double* x){
    int info = 0;
    F77_CALL(dgetrf)(&n,&n,x,&n,ipiv,&info); 
    if(info)return info;
    F77_CALL(dgetri)(&n,x,&n,ipiv,work,&n,&info);    
    return info;
}

// invert matrix by LU-factorization
void R_inv(int* _n,double* x,double* xinv,int* info){
    int n = *_n;
    int* ipiv = (int*) calloc(n, sizeof(int));if(!ipiv)return;
    double* work = (double*) calloc(n, sizeof(double));
	if(!work){free(ipiv);return;}
    memcpy(xinv,x,n * n * sizeof(double)); //dgetrf changes 'x'
    F77_CALL(dgetrf)(&n,&n,xinv,&n,ipiv,info); 
    if(!*info)
		F77_CALL(dgetri)(&n,xinv,&n,ipiv,work,&n,info);  
	free(ipiv);
	free(work);
}


//// use dposv to invert positive symmetric matrix, still won't work
//char uplo='U';
//Matrix<> tmp=VV - VB * t(VB);
//Matrix <double,Row,Concrete> B(2, 1); B=Vr;
//F77_CALL(dposv)(&uplo, &_n, &_nrhs, tmp.getArray(), &_n, B.getArray(), &_n, &info); 
//if (info==0) {
//    //success
//    crit= (t(Vr) * B) (0,0);
//} else {
//    crit=R_NegInf;
//}           
//
//Matrix<> tmp=VV - VB * t(VB);
//Matrix <double,Row,Concrete> tmpinv(2, 2);
//R_inv(&_n, tmp.getArray(), tmpinv.getArray(), &info);
//crit = (t(Vr) * (tmpinv * Vr))(0,0);   
//


// return generalized inverse of 'x' in 'xinv' and the singular values in 'sval' if not NULL
// x = U_(m x m) * S_(m x n) * V'_(n x n)  
// S = nonnegative at 1..min(M,N) diagonal
// D = MP(S), the Moore-Penrose inverse of S
// D = (d_ii) = 1/s_ii if s_ii!= 0, 
//                0    if s_ii = 0
//                0    elsewhere    
//MP(x) = VDU' = t(V') D t(U) _(n x m) is Moore-Penrose inverse
int ginv(double tol,int M, int N,double* x,double* xinv,double* sval){

    int i,info = 0;
    int q = MIN(M,N);
    
// buffer    
    int lwork = MAX(3*q + MAX(M,N),5*q);
    lwork = MAX(1,lwork);
	
	// u,vt,s,D,work
	int buffer_size = M*M + N*N + q + N * M + lwork; 
    double *_p = (double *) malloc(buffer_size * sizeof(double));
	if(!_p){Rprintf("Memory allocation error\n");return 1;}
	double* p = _p;
	
	double* u = p;p+=M*M;
	double* vt = p;p+=N*N;
	double* D = p;p+=N*M;
	double* s = p;p+=q;
	double* work = p;
	
// Note that dgesvd destroys its input matrix !    	
    memcpy(xinv,x,M*N*sizeof(double));
        
 // SVD decomposition   
    F77_CALL(dgesvd)("A","A",&M,&N,xinv,&M,&s[0],&u[0],&M,&vt[0],&N,&work[0],&lwork,&info FCONE FCONE);               
    if(info){
		free(_p);
		return(info);
	}
        
// limiting the condition number of the matrix:
// "small" singular values less than a certain factor times the highest singular value
// s[0] - highest singular value from LAPACK

// http://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse: cutoff = MACHINE_EPS * max(M,N) * max(singual val)
// implementation in R in ginv(): cutoff = MACHINE_EPS^.5 * max(singular value)

    if(tol < 0.0) tol = sqrt(DBL_EPSILON);

    double cutoff = MAX(tol * s[0],0.0);
//    double cutoff = ((double)MAX(M,N)) * MAX(tol * s[0],0.0);
    
    for(i = 0;i < q;i++)if(s[i] < cutoff) s[i] = 0.0;
    if(sval){
		for(i = 0;i < q;i++)
			sval[i] = s[i];
    }
 // D is in COLUMN_MAJOR order for FORTRAN   
 //   double D[N*M];// D_(N x M)
    memset(&D[0],0,N*M*sizeof(double));
    for(i = 0;i < q;i++)
        if(s[i] > 0.0)
            D[i*N + i] = 1.0/s[i];
            
// compute DU' and store in xinv on first call to dgemm
// then copy xinv -> D so that xinv can be overwitten on 2nd call to dgemm with the final result 
    double ALPHA = 1.0;
    double BETA = 0.0;
    F77_CALL(dgemm)("n", "t", &N, &M, &M, &ALPHA, &D[0], &N, &u[0],&M, &BETA,xinv,&N FCONE FCONE);
    memcpy(&D[0],xinv,N*M*sizeof(double));
    F77_CALL(dgemm)("t", "n", &N, &M, &N, &ALPHA, &vt[0],&N, &D[0], &N, &BETA,xinv,&N FCONE FCONE);
    
	free(_p);	
    return 0;
}



// test ginv from R
void R_ginv(double* _tol,int* _M, int* _N, double* _s,double* _u0,double* _u,double* _v,double* x,double* xinv){
    int i,j,info = 0;
    int M = *_M;
    int N = *_N;
    double tol = *_tol;
    int lda = MAX(1,M); // needed if subsetting 'x'
    int q = MIN(M,N);
    int ldu = M;//q ?
    int du[2] = {ldu,M};
    int ldvt = N;//q?
    int dv[2] = {ldvt,N};
    //    double *u = (double*) calloc(du[0] * du[1], sizeof(double)); // zeros
    //    double *vt = (double*) calloc(dv[0] * dv[1], sizeof(double)); // zeros
    //    double* s = (double*) calloc(q, sizeof(double)); // singular valued in desc. order
    int lwork = MAX(3*q + MAX(M,N),5*q);
    lwork = MAX(1,lwork);
    //    double* work = (double*) calloc(lwork, sizeof(double));
    double u[du[0]*du[1]];
    double vt[dv[0]*dv[1]];
    double s[N];
    double work[lwork];  

// Note that dgesvd destroys its input matrix !    	
    memcpy(xinv,x,M*N*sizeof(double));

    // A = U * SIGMA * t(V)  
    // U (MxM) if job 'A'
    // SIGMA min(M,N)
    // t(V) LDVTxN, LDVT >= N if job 'A'

    PRINTF("x contiguous in memory\n");
    for(i = 0;i < M*N;i++) PRINTF("%f ",xinv[i]);

    PRINTF("\n\n");
    
    PRINTF("x by mod-col\n");
    for(i = 0;i < M;i++){
        for(j = 0;j < N;j++){
          PRINTF("%f ",xinv[i*N + j]);
        }
    PRINTF("\n");
    } 
    PRINTF("\n\n");
    
    PRINTF("x by mod-row\n");
    for(j = 0;j < N;j++){
        for(i = 0;i < M;i++){
          PRINTF("%f ",xinv[j*M + i]);
        }
    PRINTF("\n");
    } 
    PRINTF("\n\n");
    
    PRINTF("M(%i) N(%i) du(%i %i) dv(%i %i) min(M,N)(%i) lda(%i) lwork(%i)\n",M,N,du[0],du[1],dv[0],dv[1],q,lda,lwork);
    F77_CALL(dgesvd)("A","A",(int*)&M,(int*)&N,xinv,(int*)&lda,(double*)&s[0],
    (double*)&u[0],(int*)&ldu,(double*)&vt[0],(int*)&ldvt,(double*)&work[0],
    (int*)&lwork,(int*)&info FCONE FCONE);

    memcpy(_u0,&u[0],M*M*sizeof(double));
    memcpy(_v,&vt[0],N*N*sizeof(double));
    memcpy(_s,&s[0],N*sizeof(double));

         
// X (mxn)          
// svd(X) = USV' is the singular value decomposition of X
// U(mxm) has orthonormal columns
// S(mxn) has first min(n,p) diagonals as singular values, 0 ow.
// V(nxn) orthogonal       
// D = S* = ginv(S)
// d_ii = 1/s_ii if s_ii!=0, 0 if s_ii = 0             
// ginv = VDU'

// http://mathforum.org/kb/thread.jspa?threadID=529857&messageID=1606150 ??
// limiting the condition number of the matrix:
// "small" singular values less than a certain factor times the highest singular value

    if(tol < 0.0) tol = sqrt(DBL_EPSILON);
    double cutoff = MAX(tol * s[0],0.0);
    PRINTF("singular values ");for(i = 0;i < q;i++)PRINTF("%f ",s[i]);PRINTF("\n\n");
    int nzero = 0;
    for(i = 0;i < q;i++){
        if(s[i] > cutoff){
            s[i] = 1.0/s[i];
        }else{
            s[i] = 0.0;
            nzero++;
        }
    }
 
    PRINTF("u_(Mxq) stored in column order\n");
    PRINTF("u_(Mxq) printed by mod row-dim\n");
    for(i = 0;i < M;i++){
        for(j = 0;j < q;j++){
            PRINTF("%f ",u[i*q + j]);  
        }
        PRINTF("\n");
    }
    PRINTF("\n");
    
    PRINTF("u_(Mxq) printed by mod col-dim\n");
    for(i = 0;i < M;i++){
        for(j = 0;j < q;j++){
            PRINTF("%f ",u[j*M + i]);  
        }
        PRINTF("\n");
    }   
    PRINTF("\n");

    PRINTF("vt_(qxN) stored in column order\n");
    PRINTF("vt_(qxN) printed by mod row-dim\n");
    for(i = 0;i < N;i++){
        for(j = 0;j < q;j++){
            PRINTF("%f ",vt[i*q + j]);  
        }
        PRINTF("\n");
    }
    PRINTF("\n");
    PRINTF("vt_(qxN) printed by mod col-dim\n");
    for(i = 0;i < N;i++){
        for(j = 0;j < q;j++){
            PRINTF("%f ",vt[j*N + i]);  
        }
        PRINTF("\n");
    } 
    PRINTF("\n");
    
   
    for(i = 0;i < M;i++){
        for(j = 0;j < q;j++){
            u[i*q + j]*=s[i];  
        }
        PRINTF("\n");
    } 
    
   // q = q
   // u  Mxq
   // vt qxN   
// VDU' = Nxq qxq 
// VDU' = V(UD')' = V(UD)' = t(V') 
// VDU' = (UD'V')'
   double ALPHA = 1.0;
   double BETA = 0.0;

  
   F77_CALL(dgemm)("t", "t", &N, &M, &q, &ALPHA, &vt[0], &dv[0], &u[0],&du[1] , &BETA,xinv,&lda FCONE FCONE);
   PRINTF("x contiguous in memory\n");for(i = 0;i < M*N;i++) PRINTF("%f ",xinv[i]);PRINTF("\n\n");
   F77_CALL(dgemm)("n", "n", &N, &M, &q, &ALPHA, &u[0],&du[1] , &vt[0], &dv[0],&BETA,xinv,&lda FCONE FCONE);
   PRINTF("x contiguous in memory\n");for(i = 0;i < M*N;i++) PRINTF("%f ",xinv[i]);PRINTF("\n\n");

   memcpy(_u,&u[0],M*M*sizeof(double));
}

// Interface to LAPACK's dgesvd svd routine
// see LAPACK documentation for more details: http://www.netlib.org/lapack/double/dgesvd.f
// to be called from C or R via .C() hence the int*- and double*-only arguments
// matrices are stored in column-major order (transpose in row-major)
// x (m x n) - DESTROYED !
// s length min(m,n) for singular values sorted
// u (m x m) at most
// vt(n x n) at most

void C_dgesvd(int* jobu,int* jobv,int* nrx,int* ncx,
double* x,double* s,double* u,double* vt,int* info)
{
	char const jobs[] = "NOSA";
	char JOBU[2];JOBU[0] = jobs[*jobu];JOBU[1] = '\0';
	char JOBV[2];JOBV[0] = jobs[*jobv];JOBV[1] = '\0';
	// Rprintf("jobi(%i %i) jobs(%s,%s)\n",*jobu,*jobv,&JOBU[0],&JOBV[0]);
	
	// set leading dimensions to default values no matrices are submatrices here
	int ldx = MAX(1,*nrx); 
	int ldu = 1;
	if((JOBU[0] == 'S') || (JOBU[0] == 'A'))
		ldu = *nrx;
	int ldvt = 1;
	if(JOBV[0] == 'S')
		ldvt = MIN(*nrx,*ncx); 
	else if(JOBV[0] == 'A')		
		ldvt = *ncx;

    // Rprintf("n=%i p=%i ldx=%i ldu=%i ldvt=%i\n",*nrx,*ncx,ldx,ldu,ldvt);
	// dgesvd
    int lwork = -1;
	double _work;
    F77_CALL(dgesvd)(JOBU, JOBV, nrx,ncx,x,&ldx,s,u,&ldu, vt,&ldvt,&_work,&lwork,info FCONE FCONE);
	if(*info){
		Rprintf("Illegal arguments to Lapack routine '%s' returning error code %d", "dgesvd" ,*info);
		return;
	}
	lwork = (int)_work;
    double *work = (double *) malloc(lwork * sizeof(double));
    F77_CALL(dgesvd)(JOBU, JOBV, nrx,ncx,x,&ldx,s,u,&ldu, vt,&ldvt,work,&lwork,info FCONE FCONE);
	free(work);
	if(*info){
		Rprintf("error code %d from Lapack routine '%s'", *info, "dgesvd");
		//return;
	}		
}



// Interface to LAPACK's dgesdd()
// to be called from C or R via .C() hence the int*- and double*-only arguments
// matrices are stored in column-major order (transpose in row-major)
// x (m x n) - DESTROYED !
// s length min(m,n) for singular values sorted
// u (m x m) at most
// vt(n x n) at most
// see LAPACK documentation for more details: http://www.netlib.org/lapack/double/dgesdd.f

void C_dgesdd(int* jobu,int* nrx,int* ncx,
double* x,double* s,double* u,double* vt,int* info)
{
	char const jobs[] = "NOSA";
	char JOBU[2];JOBU[0] = jobs[*jobu];JOBU[1] = '\0';
	//Rprintf("jobi(%i) jobs(%s)\n",*jobu,&JOBU[0]);
	
	// set leading dimensions to default values no matrices are submatrices here
	int ldx = MAX(1,*nrx); 
	int ldu = 1;
	int ldvt = 1;
	
	if(JOBU[0] == 'S'){ 
		ldu = *nrx;
		ldvt = MIN(*nrx,*ncx); 
	}else if(JOBU[0] == 'A'){ 
		ldu = *nrx;
		ldvt = *ncx;
	}
    // Rprintf("n=%i p=%i ldx=%i ldu=%i ldvt=%i\n",*nrx,*ncx,ldx,ldu,ldvt);
		
	// dgesdd
    int lwork = -1;
	double _work;
	int *iwork = (int*) malloc(8*(size_t)(MIN(*nrx,*ncx) * sizeof(int)));
    F77_CALL(dgesdd)(JOBU,nrx,ncx,x,&ldx,s,u,&ldu,vt, &ldvt,&_work, &lwork, iwork, info FCONE);
	if(*info){
		Rprintf("Illegal arguments to Lapack routine '%s' returning error code %d", "dgesdd" , *info);
		free(iwork);
		return;
	}
	lwork = (int)_work;
    double *work = (double *) malloc(lwork * sizeof(double));
    F77_CALL(dgesdd)(JOBU,nrx, ncx, x, &ldx,s,u, &ldu,vt, &ldvt,work, &lwork, iwork, info FCONE);
	free(work);
	free(iwork);
	if(*info){
		Rprintf("error code %d from Lapack routine '%s'", *info, "dgesdd");
	}	
}


// Interface to LAPACK's dgesvd svd routine to get singular values only
// see LAPACK documentation for more details: http://www.netlib.org/lapack/double/dgesvd.f
// to be called from C or R via .C() hence the int*- and double*-only arguments
// _x is stored in column-major order (transpose in row-major)
// _x (m x n) - NOT DESTROYED, a copy is made internally !
// s length min(m,n) for singular values returned
// singular values are sorted in decreasing order
void C_singval_dgesvd(int* nrx,int* ncx,double* _x,double* s,int* info)
{

	double* u = NULL;
	double* vt = NULL;
	
	// set leading dimensions to default values no matrices are submatrices here
	int ldx = MAX(1,*nrx); 
	int ldu = 1;
	int ldvt = 1;

	
	// dgesvd
    int lwork = -1;
	double _work;
    F77_CALL(dgesvd)("N", "N", nrx,ncx,NULL,&ldx,s,NULL,&ldu,NULL,&ldvt,&_work,&lwork,info FCONE FCONE);
	if(*info){
		Rprintf("Illegal arguments to Lapack routine '%s' returning error code %d", "dgesvd" ,*info);
		return;
	}
	lwork = (int)_work;
    double *work = (double *) malloc(lwork * sizeof(double));
	double* x = (double *) malloc((size_t)(*nrx * *ncx) * sizeof(double));
	memcpy(x,_x,(size_t)(*nrx * *ncx) * sizeof(double));
    F77_CALL(dgesvd)("N", "N", nrx,ncx,x,&ldx,s,u,&ldu, vt,&ldvt,work,&lwork,info FCONE FCONE);
	free(work);
	free(x);
	if(*info){
		Rprintf("error code %d from Lapack routine '%s'", *info, "dgesvd");
		//return;
	}		
}

// x    : input symmetric pos. def. matrix, will be overwritten
// xinv : output inv(x)
int ldlinv(int n,double* x,double* xinv){

		int* ipiv = (int *)malloc((size_t)n * sizeof(int));
		if(!ipiv){
			Rprintf("Unable to allcoate %i bytes in function %s\n",n * sizeof(int),"ldlSolve");
			return 1;
		}	
		int LWORK = -1;
		double WORK;
		int info = 0;		
		F77_CALL(dsytrf)("U",&n,x,&n,ipiv,&WORK,&LWORK,&info FCONE); 
		if(info){
			free(ipiv);
			Rprintf("error code %d from Lapack routine '%s'\n", info, "dsytrf");
			return 1;
		}	
		
		LWORK = (int)WORK;
		double* work = (double*)malloc((size_t)(LWORK * sizeof(double)));
		if(!work){
			free(ipiv);
			Rprintf("Unable to allcoate %i bytes in function %s\n",LWORK * sizeof(double),"ldl_inv");
			return 1;		
		}
		F77_CALL(dsytrf)("U",&n,x,&n,ipiv,work,&LWORK,&info FCONE); 
		if(info){
			free(ipiv);
			free(work);
			Rprintf("error code %d from Lapack routine '%s'\n", info, "dsytrf");
			return 1;
		}
		// set xinv to identity and solve for it
		memset(xinv,0,n*n*sizeof(double));
		for(int i = 0;i < n;i++)xinv[i*n+i] = 1.0;
		F77_CALL(dsytrs)("U",&n,&n, x, &n, ipiv, xinv,&n,&info FCONE);
		if(info){
			Rprintf("error code %d from Lapack routine '%s'\n", info, "dsytrs");	
			free(ipiv);
			free(work);
			return 1;
		}	
		free(ipiv);
		free(work);		
		return 0;
}
