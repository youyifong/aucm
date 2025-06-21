// Author: Krisztian Sebestyen <ksebestyen@gmail.com>

#include <R.h>


#ifndef MATRIX_UTILS
#define MATRIX_UTILS
// note: do not declare integers as character constants if longer than 1 character as it is dependent
// on machine byte order 'Little-' vs 'Big-Endian'
typedef enum {ADD = '+',SUB = '-',MULT = '*',DIV = '/'}op_type;
//typedef enum {COLUMN_MAJOR = 0,ROW_MAJOR = 1}MAJOR_ORDER;
#endif

// z = c1*x 'op' c2*y
// 'int op' will be silently typecasted to 'op_type'
void vec_op(double* z,double c1,double* x,int op,double c2,double* y,int n);
void matprod(double *x, int* nrx, int* ncx,double *y, int* nry, int* ncy, double *z);
void crossprod(double *x, int* nrx, int* ncx,double *y, int* nry, int* ncy, double *z);
void tcrossprod(double *x, int* nrx, int* ncx,double *y, int* nry, int* ncy, double *z);			  

// y = x[ix,jx], major_ = "r" or "R" for row-major order , "c" or "C" for column-major order  
void get_sub_matrix(int* ix,int* jx,const char* _major_x,int nrx,int ncx,double* x,const char* _major_sub,int nrs,int ncs,double* y);
void print_matrix(int m,int n,double* x);

int luinv0(int n,int* ipiv,double* work,double* x);
int luinv(int n,int* ibuf,double* dbuf,double* x,double* xinv);
int ginv(double tol,int M, int N,double* x,double* xinv,double* singular_values);
int ldlinv(int n,double* x,double* xinv);

void C_dgesvd(int* jobu,int* jobv,int* nrx,int* ncx,double* x,double* s,double* u,double* vt,int* info);
void C_dgesdd(int* jobu,int* nrx,int* ncx,double* x,double* s,double* u,double* vt,int* info);
void C_singval_dgesvd(int* nrx,int* ncx,double* _x,double* s,int* info);

void dxd_(int* n, double* d1,double* _x,double* d2,double* _y);
void dxd(int* n, double* d1,double* _x,double* d2,double* _y);

void lower_trap(int nrx,int ncx,double* x,double* diag,int k,double* L);
void upper_trap(int nrx,int ncx,double* x,double* diag,int k,double* U);

void rowperm_ipiv(int* _n,int* ipiv,int* perm);
void invperm(int* n,int* p,int* ip);
void row_PL(int* _n, int* rowperm,int* iP,double* dP);

// void rcbind(double* x,int nrx,int ncx,int times,int each,double* y);
// void rrbind(double* x,int nrx,int ncx,int times,int each,double* y);
void rrbind(double* _x,int nrx,int ncx,int times,int each,int* vec_each,double* _y);
void rcbind(double* _x,int nrx,int ncx,int times,int each,int* vec_each,double* _y);

