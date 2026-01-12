#ifndef FCONE
# define FCONE
#endif


#include <R.h>
#include <R_ext/Lapack.h>
void C_dgetrf(int* M, int* N, double* A, int* LDA, int* IPIV, int*INFO){
    F77_CALL(dgetrf)(M, N, A, LDA, IPIV, INFO );
}
void C_dgetrs(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B,int* LDB,int*INFO){
    F77_CALL(dgetrs)("N", N, NRHS, A, LDA, IPIV, B, LDB, INFO FCONE);
}
void C_dgesv(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO){
    F77_CALL(dgesv)( N, NRHS, A, LDA, IPIV, B, LDB, INFO);
}
