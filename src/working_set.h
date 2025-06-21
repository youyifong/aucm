


// functions to select working set. the naming convention is: 
// 'r' - random, 'v' - from violators, 'v2' - from two different sets of violators if available at current iteration
// 'w' - weighted selection, 'g' - weights proportional to gradient
int getWorkingSet_rv2wg(int nB,int* B,int nV0,int nVC,int n1n2,double* df,int* v,int* which_v,int* perm,double* prob);
int getWorkingSet_rvwg(int nB,int* B,int nV0,int nVC,int n1n2,double* df,int* v,int* which_v,int* perm,double* prob);
int getWorkingSet_rv(int nB,int* B,int nV0,int nVC,int n1n2,int* v,int* which_v,int* buf);
int getWorkingSet_v(int nB,int* B,int nV0,int nVC,int n1n2,double* a,double * df,int* v,int* which_v,int* perm,double* dbuf);
int getWorkingSet_greedy_0C(int nB,int* B,int n1n2,double* a,double * df,double* C,int* v,int* which_v,int* perm,double* dbuf);
int getWorkingSet_v2(int nB,int* B,int nV0,int nVC,int n1n2,double* a,double * df,int* v,int* which_v,int* perm,double* dbuf);
int getWorkingSet_rv2(int nB,int* B,int nV0,int nVC,int n1n2,int* v,int* which_v,int* perm);
//int getWorkingSet_rvwr(int nB,int* B,int nV0,int nVC,int n1n2,double* df,int* v,int* which_v,int* ibuf,int* ibuf2,double* dbuf);
  
// identify violators at current iteration, getKKTViolators_0C is used in code
int getKKTViolators_0C(int n, double* a, double* df, double* C,int* nV0,int* nVC,int* v);
//int getKKTViolators_0C(int n, double* a, double* df, double* C,int* nV0,int* nVC,int* v,int* which_v,int* which_v0,int* which_vc);
int getKKTViolators(int n, double* a, double* grad, double* C,int* v);

// KKT terminating criterion: getEpsilonKKTV() is used in the code 
double getEpsilonKKTV(int n1n2,double* a,double* grad,double* C);
double getEpsilonKKT(int n1n2,double* a,double* grad, double* C);  
