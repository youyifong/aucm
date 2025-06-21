#include "minQuad_QP.h"

int iter_combn(int* info,int n, int m,int* x);
double* get_hypercube(int dim,double lo,double up);
// void next_hypercube_vertex(int i,int k,double lo,double up,double* x);
// int get_sub_solution(MNQD_QP* qp,int np,int* px,int nq,int* qx,double* _x,double* x);
int optimize_qp_exhaustive(MNQD_QP* qp, MNQD_QP* sub, double* solution,int flag,int trace);
int optimize_qp_exhaustive_fast(MNQD_QP* qp,MNQD_QP* sub, double* solution,int flag,int trace);

