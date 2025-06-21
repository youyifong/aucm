#ifndef MINQUAD_MNQD_QP_PROGRAM
#define MINQUAD_MNQD_QP_PROGRAM
#include <stdlib.h>
/* 
  solve the quadratic programming problem
  
  argmin 0.5x'Hx + b'x
  
  subject to eq_vec * x = eq_const
  
  xl <= x <= xu
 
  The linear constraint vector eq_vec can only have -1/+1 as entries as required by hideo 
*/


// Note: Each sub-problem in the decomposition algorithm in minQuad 
// needs to satisfy the equality constraint (if given)


typedef struct minQuad_quadratic_program {
  int  N;            /* allocated size, N >= n */
  int   n;           /* number of variables */
  double *H;          /* hessian of objective */
  double* iH;         /* inverse of Hessian */
  double *b;          /* linear term in objective */
  double *x0;         /* initial value for variables */
  double *xl,*xu;     /* box constraints */
  int n_constr;    
  double* mat_constr;  // (n_constr x n) matrix 'A' in column major order: lhs <= Ax <= rhs
  double* lhs_constr; /* length n_constr*/
  double* rhs_constr; /* length n_constr  */
  double* singular_values; /* determines 'rank' if 'rank' is not given */
  double* dbuf;      // buffer of length 'N'
  int* ibuf;         // buffer of length 'N' used in luinv()
  int rank; /* positive definite qp-optimizers can only handle this many variables */
  int max_order;
  double machine_double_eps;
}MNQD_QP;
#endif

double objective_QP(MNQD_QP* qp,double* x);
void free_QP(MNQD_QP* qp);
void print_QP(MNQD_QP* qp);
void get_sub_Q(int n1,int n2,int n1n2,int type_H,void* H,int nB,int* B,double* sub_Q);
int get_QP(int n1,int n2,int type_H,void* H,double* Ha,double* a,double* b,int nB,int* B,double* C,
    int n_constr,double* mat_constr,double* lhs_constr,double* rhs_constr,
    int do_inverse,int check_min_existence,int rank,MNQD_QP* qp);
int get_sub_QP(int nix,int* ix,int do_inverse,int check_min_existence,MNQD_QP* Qp,MNQD_QP* sub);
int init_QP(int n,int n_constr,double machine_double_eps,MNQD_QP* qp);
