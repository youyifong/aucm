#ifndef LOQO_h
#define LOQO_h
typedef enum {CONVERGED = 0,PRIMAL_DUAL_INFEASIBLE = 1,PRIMAL_INFEASIBLE = 2,
DUAL_INFEASIBLE = 3, PRIMAL_INFEASIBLE_SOLUTION = 4,
PRIMAL_UNBOUNDED = 5, DUAL_UNBOUNDED = 6,
SLOW_CONVERGENCE = 7,MEMORY_ALLOCATION_FAILURE = 8, LAPACK_ERROR = 9,LOGIC_ERROR = 10}LOQO_STATUS;
extern const char* LOQO_STATUS_NAMES[11];
typedef enum {LDL = 0,LU = 1,QR = 2}LAPACK_SOLVER;
#endif
			 
void loqo(int* _n,double* linear, double* hessian, double* l, double* u,
int* _m, double* constr_mat, double* constr_vec1,  double* constr_vec2, 
int* lin_solver,
double* sigfig_max, int* maxiter, double* margin, double* bound, double* inf,
int* verb,int* _dbuffer_size,double* _dbuffer,
double* primal,double* dual,double* _primal_obj,double* _dual_obj,
int* counter,int* convergence,int* error_code);

void get_buffer_size_loqo(int solver,int update_solver,int m,int n,int* dsize,int* isize);
