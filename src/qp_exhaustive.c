#include "qp_exhaustive.h"
#include "minQuad.h"
#include "matrix.h"
#include <R_ext/Lapack.h>
//#include <R_ext/Applic.h>

#include <R.h>
#define PRINTF Rprintf

#define	MAX(A,B)	((A) > (B) ? (A) : (B))
#define	MIN(A,B)	((A) < (B) ? (A) : (B))
#define PROJ(L,X,U) MAX((L),MIN((X),(U))) // median(L,X,U)

#include <math.h>
#include <float.h>
#define _Inf DBL_MAX


// or define a 'bit' type via bit fields for true and false and replace code pertaining to 'IS_FEASIBLE_X' in optimize_QP()
// struct {unsigned int false:1; unsigned int true:1;} boolean;
// boolean.false = 1; // 1-bit value of 0 or 1
#include <stdbool.h> //C99
#ifndef bool
typedef enum { false, true } bool;
#endif




int is_subset(int nx,int* x,int ny,int* y){
    if(nx == 0) return 1;
    if(!x) return 1;
    if(ny == 0) return 0;
    if(!y) return 0;
    
    for(int i = 0;i < nx;i++){
        int xi_matched = 0;
        for(int j = 0;j < ny;j++){
            if(x[i] == y[j]){xi_matched = 1;break;}
        }
        if(!xi_matched) return 0;
    }
    return 1;    
}



// R's combn() changed to an 'iterator'
// Originally written in R hence 1-based indexing adapted to C's 0-based indexing
// works on set {1,..,n} as opposed R's {a1,..,an} 
// output: 'x' with side effects in 'info' needed to at next call
// get next (info,x) by passing previous (info,x) 
void R_iter_combn(int* _n, int* _m,int* info,int* x){
//int comb_iter(int n, int m,int* info,int* x){

    int n = *_n;
    int m = *_m;
    int i = info[0];
    int e = info[1];
    int h = info[2];
    // at iteration 1, initialize as in R and return (1,2,..,m)    
    if(i <= 1){
        for(int j = 0;j < m;j++)x[j] = j + 1;
        info[0] = 2;
        info[1] = 0;
        info[2] = m;
        return;
    }

    if(x[0] == n - m + 1) return;
    if (e < n - h) {
        h = 1;
        e = x[m - 1];
        int j = 1;
        x[m - h + j - 1] = e + j;        
    }
    else {
        e = x[(m - 1) - h];
        h++;
        for(int j = 1;j < h + 1;j++) x[m - h + j - 1] = e + j;
    }
    i++;
    info[0] = i;info[1] = e;info[2] = h;
    return;
} 

// generate C(n,m) combinations
// return in 'x' the 'next' m-combination of {1,..,n}
// pass info = {0,0,0 initially} to get first 'x'
// at second call pass the 'pair' ('info','x') from the previous call
// as the necessary state to generate the next 'pair' (info,x)
int iter_combn(int* info,int _n, int _m,int* cx){

    int _i = info[0],_e = info[1],_h = info[2];
    if(_i <= 1){
        for(int _j = 0;_j < _m;_j++)cx[_j] = _j + 1;
        info[0] = 2;
        info[1] = 0;
        info[2] = _m;
        return 1;
    }else{
        if(cx[0] == _n - _m + 1) return 0;
        if (_e < _n - _h) {
            _h = 1;
            _e = cx[_m - 1];
            int _j = 1;
            cx[_m - _h + _j - 1] = _e + _j;        
        }
        else {
            _e = cx[(_m - 1) - _h];
            _h++;
            for(int _j = 1;_j < _h + 1;_j++) cx[_m - _h + _j - 1] = _e + _j;
        }
        _i++;
        info[0] = _i;info[1] = _e;info[2] = _h;
    }
    return 1;
} 

// iterator to traverse vertices of hypercube
void R_next_hypercube_vertex(int* _i,int* _k,double* _lo,double* _up,double* x){
    int i = *_i;int k = *_k;double lo = *_lo;double up = *_up;
    if(k <= 0) return;
    int r = i - 1;    
    for(int j = 0; j < k; j++) {
        int r0 = r / 2; 
        if(r - 2*r0 == 0) x[j] = lo; else x[j] = up;
        r = r0;
    }
}

// pass i = 1 .. 2^k to return 'next' vertex in 'x' with labels (lo,up)
void next_hypercube_vertex(int i,int k,double lo,double up,double* x){
    if(k <= 0) return;
    int r = i - 1;    
    for(int j = 0; j < k; j++) {
        int r0 = r / 2; 
        if(r - 2*r0 == 0) x[j] = lo; else x[j] = up;
        r = r0;
    }
}

// m x n matrix
double* get_hypercube(int n,double lo,double up){
    if(n <= 0) return NULL;
    int m = 1 << n;
    double* x = (double*)calloc(m * n , sizeof(double));
    if(!x) return NULL;
    for(int i = 0;i < m;i++){
        int r = i;
        for(int j = 0; j < n; j++) {
            int r0 = r / 2; 
            if(r - 2*r0 == 0) x[i*n + j] = lo; else x[i*n + j] = up;
            r = r0;
        }
    }
    return x;
}
//matrix(.C('R_get_hypercube',as.integer(3),0.0,1.0,double(3*2^3))[[4]],3,8)
void R_get_hypercube(int *_n,double *_lo,double *_up,double* x){
    int n = *_n;double lo = *_lo,up = *_up;
    if(n <= 0) return;
    int m = 1 << n;
    for(int i = 0;i < m;i++){
        int r = i;
        for(int j = 0; j < n; j++) {
            int r0 = r / 2; 
            if(r - 2*r0 == 0) x[i*n + j] = lo; else x[i*n + j] = up;
            r = r0;
        }
    }
}


//  A = Q[B,B]
// bB = b[B]+Q[B,-B]a[-B]
// _A = inv(A[px,px])
//  b = {bB[px] + A[px,-px]x[-px])}
//  y = - _A _b
int get_sub_solution(MNQD_QP* qp,MNQD_QP* sub,int np,int* px,int nq,int* qx,double* _x,double* x,int trace){
// linear term in restricted objective function
    double (*H)[qp->n] = (double (*)[qp->n])qp->H;    
    for(int k = 0; k < qp->n; k++)x[k] = _x[k];
    for(int i = 0; i < np; i++) {
        int ki = px[i];
        sub->b[i] = qp->b[ki]; 
        for(int j = 0; j < nq; j++){
            int kj = j;if(qx)kj = qx[j];
            sub->b[i] += H[ki][kj] * x[kj];
        }        
    }          
//  x[px] = -inv(Q[B[px],B[px]])bB[px]    
    double x0[np],none = -1.0,zero = 0.0;int inc = 1;
    F77_CALL(dsymv)("U", &np, &none, sub->iH, &np, sub->b, &inc, &zero, &x0[0], &inc FCONE);
    for(int i = 0; i < np; i++)x[px[i]] = x0[i];
    // double (*ih)[np] = (double (*)[np])(sub->iH);
    // for(int i = 0; i < np; i++) {
        // x[px[i]] = 0.0;
        // for(int j = 0; j < np; j++) x[px[i]] -= ih[i][j] * sub->b[j];
    // }    
    if(trace){
        PRINTF("sub-problem\n");
        PRINTF("px ");for(int j = 0; j < np; j++)PRINTF("%i ",px[j]);PRINTF("\n");        
        if(qx){PRINTF("qx ");for(int j = 0; j < nq; j++)PRINTF("%i ",qx[j]);PRINTF("\n");}        
//        PRINTF("subinverse\n");print_matrix(np,np,sub->iH);
        PRINTF("b ");for(int j = 0; j < np; j++)PRINTF("%+.3e ",sub->b[j]);PRINTF("\n");        
        PRINTF("x ");for(int j = 0; j < qp->n; j++)PRINTF("%+.3e ",x[j]);PRINTF("\n");
        for(int i = 0; i < np; i++) x0[i] = _x[i];
        PRINTF("sub_f(x)=%.8f\n",objective_QP(sub,x0));
        PRINTF("sub-problem end\n");        
    }    
    for(int k = 0; k < qp->n; k++){ 
        if((x[k] < qp->xl[k]) || (qp->xu[k] < x[k])){
            return 0;
        }
    }  
    return 1;
}



int optimize_qp_exhaustive_fast(MNQD_QP* qp,MNQD_QP* sub, double* solution,int flag,int trace)
{
    int n = qp->n;
    double* x0 = NULL;
    double* x = solution;
    double best_fx = _Inf;
    double x1[n];
    unsigned long nHV = 1<<n; // process 2^n vertices of hypercube  
    double sqrt_eps = sqrt(qp->machine_double_eps);
    int C_n_m_FEASIBLE[n];memset(&C_n_m_FEASIBLE[0],0,n * sizeof(int));

        
    if(trace) PRINTF("\n");    
    if(trace > 1) print_QP(qp);
    
//  unrestricted estimate    
//    double zero = 0.0,none = -1.0;int inc = 1;
//    F77_CALL(dsymv)("U", &n, &none, qp->iH, &n, qp->b, &inc, &zero, x, &inc);
    // for(int i = 0; i < qp->n; i++) {
        // x[i] = 0.0;
        // for(int j = 0; j < qp->n; j++) x[i] -= qp->iH[i * n + j] * qp->b[j];
    // }     
    // if(trace > 1){
        // PRINTF("Unrestricted estimate: f(x)=%+.3e x=",objective_QP(qp,x));
        // for(int j = 0; j < n; j++) PRINTF("%+.5e ",x[j]);PRINTF("\n");
    // }
    // int nv = 0;
    // int vx[n];int _which_vx[n];int* which_vx = &_which_vx[0];
    // for(int k = 0; k < n; k++){ 
        // if((x[k] < qp->xl[k]) || (qp->xu[k] < x[k])){
            // x[k] = PROJ(qp->xl[k],x[k],qp->xu[k]);
            // vx[k] = 1;
            // which_vx[nv] = k;
            // nv++;
        // }
    // }
// //    if(!nv) return 1;    
// //    if(!nv) which_vx = NULL;
    
    
    // double fx = 0.0;
    // double (*H)[qp->n] = (double (*)[qp->n])qp->H;
    // for(int i = 0; i < qp->n; i++) {
        // if(x[i] <= qp->machine_double_eps) continue;
        // for(int j = 0; j < qp->n; j++) {
            // if(x[j] <= qp->machine_double_eps) continue;
            // fx   +=  x[i] * H[i][j] * x[j];
        // }
    // }        
    // fx *= 0.5;for(int i = 0; i < qp->n; i++) fx += qp->b[i] * x[i];
    // best_fx = fx; 
    // for(int k = 0; k < n; k++) x1[k] = x[k]; // store 'x' as it may be different from vertices of hypercube
        
// buffers    
    double _x[n];x = &_x[0];        
    double _x0[n];x0 = &_x0[0]; 
  
// iter_combn expects a combination in {1..n} to generate the next, except at first call
// (px,qx) are 0-based combinations
// {px} union {qx} = {zx} = {0..n-1}
// {px} intersect {qx} = empty    
// qx may be empty, px is never empty

    int zx[n];for(int k = 0;k < n;k++)zx[k] = k; 
    for(int m = qp->max_order; 1 <= m; m--)
    {
        int cx[m]; // stores m-combination that will be passed to hypercube vertex iterator
        int px[m]; // C-based index of m-combination, px = cx - 1
        int _qx[MAX(1,n - m)]; // compliment of px, NULL if px is the full length 'n' combination
        int* qx = n > m ? &_qx[0] : NULL; 

        if(trace > 1)PRINTF("Combinations Cn(%i,%i)\n",n,m);      
 
/////////////////////////// get a combination of C(n,m) /////////////////////////////////////////////// 
        int _m = m;int _n = n;
        int info[3] = {0,0,0}; // int info[3] = {2,0,_m};
        for(int _j = 0;_j < _m;_j++)cx[_j] = _j + 1;
        int end_combinations = 0;
        while(!end_combinations){                
            int _i = info[0],_e = info[1],_h = info[2];
            if(_i <= 1){
                info[0] = 2;info[1] = 0;info[2] = _m;
            }else{
                if(cx[0] == _n - _m + 1){ 
                    end_combinations =  1;
                    break;
                }
                if (_e < _n - _h) {
                    _h = 1;
                    _e = cx[_m - 1];
                    int _j = 1;
                    cx[_m - _h + _j - 1] = _e + _j;        
                }
                else {
                    _e = cx[(_m - 1) - _h];
                    _h++;
                    for(int _j = 1;_j < _h + 1;_j++) cx[_m - _h + _j - 1] = _e + _j;
                }
                _i++;
                info[0] = _i;info[1] = _e;info[2] = _h;
            }
//////////////////////////////////////////////////////////////////////////////////////////////

            // subtract 1 for C-based indexing
            for(int k = 0;k < m;k++)px[k] = cx[k]-1;             
            if(qx){
                for(int k = 0;k < m;k++)zx[px[k]] = -1;
                int l = 0;for(int k = 0;k < n;k++)if(zx[k]!=-1) qx[l++] = zx[k];
                for(int k = 0;k < m;k++)zx[px[k]] = px[k];
            }            
            if(trace > 1){
                PRINTF("info(i=%i e=%i h=%i)\n",info[0],info[1],info[2]);
                PRINTF("px ");for(int q = 0; q < m; q++)PRINTF("%i ",px[q]);PRINTF("\n"); 
                PRINTF("qx ");for(int q = 0; q < n-m; q++)PRINTF("%i ",qx[q]);PRINTF("\n");
            }         

//            if(!get_sub_QP(m,&px[0],1,0,qp,sub)){continue;}
////////////////////////// get qp subproblem //////////////////////////////////////////////////////////    
            sub->n = m;    
            for(int i = 0;i < sub->n;i++)sub->singular_values[i] = -1.0;
            double (*H)[qp->n] = (double (*)[qp->n])qp->H;
            double (*h)[sub->n] = (double (*)[sub->n])sub->H;
            for(int i = 0;i < sub->n;i++){
                for(int j = 0;j < sub->n;j++){
                    h[i][j] = H[px[i]][px[j]];
                }
            }
            int success = 0;
            if(sub->n == 1){
                if(sub->H[0] * sub->H[0] <= sub->machine_double_eps) continue;
                success = 1;
                sub->iH[0] = 1.0/sub->H[0];
            }else if(sub->n == 2){
                double (*ih)[sub->n] = (double (*)[sub->n])sub->iH;
                double det =  h[0][0] * h[1][1] - h[0][1] * h[1][0];
                if(det <= sub->machine_double_eps){
                    break;
                }
                ih[0][1] = - h[1][0] / det;
                ih[1][0] = - h[0][1] / det;
                ih[0][0] = + h[1][1] / det;
                ih[1][1] = + h[0][0] / det;
                success = 1;     
            }
            if((!success) || (sub->n > 2)){ 
                if(luinv(sub->n,sub->ibuf,sub->dbuf,sub->H,sub->iH))
                    ginv(sqrt_eps,sub->n,sub->n,sub->H,sub->iH,sub->singular_values);          
            }
////////////////////////////////////////////////////////////////////////////////////
           

            int nmHV = 1<<(n - m);
            for(int k = 0;k <= nmHV;k++){      
            
               if(k == nmHV){
                   for(int j = 0;j < n;j++)x0[j] = x1[j]; // process the unrestricted estimate which **may partially** be in interior of n-cube
               }else{ // get hypercube vertex, the binary representation of the integer 'k'
                   for(int j = 0; j < n; j++)x[j] = 0.0;
                   int kk = k;int j = 0;
                   do{x[j++] = (double)(kk & 1);kk = kk >> 1;}while(kk);                               
                   for(int j = 0;j < n-m;j++)x0[qx[j]] = x[j];               
               }
               
               if(trace > 2){PRINTF("\nsub-projections of hypercube vertex: ");for(int j = 0; j < n; j++)PRINTF("%+.3e ",x0[j]);PRINTF("\n");}                 
                
                          
///////////////////////////// get sub-solution //////////////////////////////                
                int success = 1;
//              if(!get_sub_solution(qp,sub,m,&px[0],n-m,qx,x0,x,MAX(0,trace - 1)))
            // linear term in restricted objective function
                int np = m;int nq = n-m;
                double (*H)[qp->n] = (double (*)[qp->n])qp->H;
                for(int j = 0; j < qp->n; j++)x[j] = x0[j];
                for(int i = 0; i < np; i++) {
                    sub->b[i] = qp->b[px[i]]; 
                    for(int j = 0; j < nq; j++)
                        sub->b[i] += H[px[i]][qx[j]] * x[qx[j]];
                }          
            //  x[px] = -inv(Q[B[px],B[px]])bB[px]    
                // double x_m[np],none = -1.0,zero = 0.0;int inc = 1;
                // F77_CALL(dsymv)("U", &np, &none, sub->iH, &np, sub->b, &inc, &zero, &x_m[0], &inc);
                // for(int i = 0; i < np; i++)x[px[i]] = x_m[i];
                double (*ih)[np] = (double (*)[np])(sub->iH);
                for(int i = 0; i < np; i++) {
                    x[px[i]] = 0.0;
                    for(int j = 0; j < np; j++) x[px[i]] -= ih[i][j] * sub->b[j];
                }    
                if(trace){
                    PRINTF("sub-problem\n");
                    PRINTF("px ");for(int j = 0; j < np; j++)PRINTF("%i ",px[j]);PRINTF("\n");        
                    if(qx){PRINTF("qx ");for(int j = 0; j < nq; j++)PRINTF("%i ",qx[j]);PRINTF("\n");}        
            //        PRINTF("subinverse\n");print_matrix(np,np,sub->iH);
                    PRINTF("b ");for(int j = 0; j < np; j++)PRINTF("%+.3e ",sub->b[j]);PRINTF("\n");        
                    PRINTF("x ");for(int j = 0; j < qp->n; j++)PRINTF("%+.3e ",x[j]);PRINTF("\n");
                    PRINTF("sub-problem end\n");        
                }    
                for(int j = 0; j < qp->n; j++){ 
                    if((x[j] < qp->xl[j]) || (qp->xu[j] < x[j])){
                        success = 0;
                        break;
                    }
                }                
///////////////////////////////////////////////////////////////////////                
                               
                if(!success){
                    if(trace > 2){
                        PRINTF("projection out of bounds: f(x)=%+.3e x=",objective_QP(qp,x));
                        for(int j = 0; j < n; j++)PRINTF("%+.3e ",x[j]);PRINTF("\n");
                    } 
                    continue;  
                }else{
                    //double fx = objective_QP(qp,x);
                    double fx = 0.0;
                    double (*H)[qp->n] = (double (*)[qp->n])qp->H;
                    for(int i = 0; i < qp->n; i++) {
                        if(x[i] <= qp->machine_double_eps) continue;
                        for(int j = 0; j < qp->n; j++) {
                            if(x[j] <= qp->machine_double_eps) continue;
                            fx   +=  x[i] * H[i][j] * x[j];
                        }
                    }        
                    fx *= 0.5;for(int i = 0; i < qp->n; i++) fx += qp->b[i] * x[i];
                    if(fx < best_fx){
                        best_fx = fx;
                        for(int j = 0; j < n; j++)solution[j] = x[j];
                    }
                    C_n_m_FEASIBLE[m-1]++;
                    if(trace > 2){
                        PRINTF("projection in bounds f(x)=%+.3e x=",fx);
                        for(int j = 0; j < n; j++)PRINTF("%+.3e ",x[j]);PRINTF("\n");
                    } 
//                    break;
                }                
//                if(C_n_m_FEASIBLE[m-1]) break;

            }//end hypercube loop   
        }// end all length-m combinations of 'n'
//        for(int j = 0; j < n; j++)if(j < m) best_px[j] = px[j];
//        if(C_n_m_FEASIBLE[m-1]) break;
    }// end C(n,1)..C(n,n) combinations       

    if(trace > 1)PRINTF("Hypercube vertices: \n ");
    for(int i = 0;i < nHV;i++){
        for(int j = 0; j < n; j++)x[j] = 0.0;
        unsigned int ii = (unsigned int)i;
        int j = 0;
        do{x[j++] = (double)(ii & 1);ii = ii >> 1;}while(ii);

        double fx = 0.0;// = objective_QP(qp,x);
        double (*H)[qp->n] = (double (*)[qp->n])qp->H;
        for(int i = 0; i < qp->n; i++) {
            if(x[i] <= qp->machine_double_eps) continue;
            for(int j = 0; j < qp->n; j++) {
                if(x[j] <= qp->machine_double_eps) continue;
                fx   +=  x[i] * H[i][j] * x[j];
            }
        }        
        fx *= 0.5;for(int i = 0; i < qp->n; i++) fx += qp->b[i] * x[i];
        
        if(trace > 1){PRINTF("f(x)=%+.3e x=",fx);for(int q = 0; q < n; q++)PRINTF("%+.3e ",x[q]);PRINTF("\n");} 
        if(fx < best_fx){
            best_fx = fx;
            for(int j = 0; j < n; j++)solution[j] = x[j];
        }
    }    
    if(trace){
//        PRINTF("Number of Feasible points ");for(int j = 0;j < n;j++) PRINTF("C(%i,%i)=%i ",n,j+1,C_n_m_FEASIBLE[j]);PRINTF("\n");
        PRINTF("projected unrestricted estimate: f(x)=%+.3f x=",x1[n]);for(int j = 0; j < n; j++)PRINTF("%.3f ",x1[j]);PRINTF("\n");
        PRINTF("best estimate: f(x)=%+.3e x=",best_fx);for(int j = 0; j < n; j++)PRINTF("%+.3e ",solution[j]);PRINTF("\n");
    } 
    return 1;
}       




int optimize_qp_exhaustive(MNQD_QP* qp, MNQD_QP* sub,double* solution,int flag,int trace)
{
    int n = qp->n;
    double* x0 = NULL;
    double* x = solution;
    double best_fx = _Inf;
    double x1[n];
    unsigned long nHV = 1<<n; // process 2^n vertices of hypercube  
    int C_n_m_FEASIBLE[n];memset(&C_n_m_FEASIBLE[0],0,n * sizeof(int));
    
    if(trace) PRINTF("\n");    
    if(trace > 1) print_QP(qp);

     double zero = 0.0,one = -1.0;int inc = 1;
// //    int best_px[n];int n_best_px;
// //    int best_qx[n];int n_best_qx;
    
// // unrestricted estimate    
    // F77_CALL(dsymv)("U", &n, &one, qp->iH, &n, qp->b, &inc, &zero, x, &inc);
    // // // for(int i = 0; i < n; i++) {
        // // // x[i] = 0.0;
        // // // for(int j = 0; j < n; j++) x[i] -= qp->iH[i * n + j] * qp->b[j];
    // // // }     
    // if(trace > 1){
        // PRINTF("Unrestricted estimate: f(x)=%+.3e x=",objective_QP(qp,x));
        // for(int j = 0; j < n; j++) PRINTF("%+.5e ",x[j]);PRINTF("\n");
    // }
    // int nv = 0;
    // int vx[n];int _which_vx[n];int* which_vx = &_which_vx[0];
    // for(int k = 0; k < n; k++){ 
        // if((x[k] < qp->xl[k]) || (qp->xu[k] < x[k])){
            // x[k] = PROJ(qp->xl[k],x[k],qp->xu[k]);
            // vx[k] = 1;
            // which_vx[nv] = k;
            // nv++;
        // }
    // }    
    // if(!nv) which_vx = NULL;
 
// //    best_fx = objective_QP(qp,x); // current best estimate in (x,best_fx)  
    // one = 1.0;zero = 0.0;inc = 1;
	// F77_CALL(dsymv)("U", &qp->n, &one, qp->H, &qp->n, x, &inc, &zero, qp->dbuf, &inc);
	// best_fx = 0.5*F77_CALL(ddot)(&qp->n, x, &inc, qp->dbuf, &inc) + F77_CALL(ddot)(&qp->n, x, &inc, qp->b, &inc);

    // for(int k = 0; k < n; k++) x1[k] = x[k]; // store 'x' as it may be different from vertices of hypercube
        
// buffers    
    double _x[n];x = &_x[0];        
    double _x0[n];x0 = &_x0[0]; 
    
// iter_combn expects a combination in {1..n} to generate the next, except at first call
// (px,qx) are 0-based combinations
// {px} union {qx} = {zx} = {0..n-1}
// {px} intersect {qx} = empty    
// qx may be empty, px is never empty

    int zx[n];for(int k = 0;k < n;k++)zx[k] = k;
    for(int m = qp->max_order; 1 <= m; m--)
    {
        int cx[m]; // stores m-combination that will be passed to hypercube vertex iterator
        int px[m]; // C-based index of m-combination, px = cx - 1
        int _qx[MAX(1,n - m)]; // compliment of px, NULL if px is the full length 'n' combination
        int* qx = n > m ? &_qx[0] : NULL; 

        if(trace > 1)PRINTF("Combinations Cn(%i,%i)\n",n,m);  
        // "backward" exclusion - only if outer m-loop is in order m = 0,1,..,n 
        // if cannot project jointly for any(!) combination in dim = p-1 then cannot jointly project in dim p > p - 1.
        // if(!(exclude % 2)){
            // if((m > 2) && (C_n_m_FEASIBLE[(m-1)-1] == 0)){
                // break;
            // }
        // }           
        int info[3] = {0,0,0};
        while(iter_combn(&info[0],n,m,&cx[0])){
            // subtract 1 for C-based indexing
            for(int k = 0;k < m;k++)px[k] = cx[k]-1;             
            if(qx){
                for(int k = 0;k < m;k++)zx[px[k]] = -1;
                int l = 0;for(int k = 0;k < n;k++)if(zx[k]!=-1) qx[l++] = zx[k];
                for(int k = 0;k < m;k++)zx[px[k]] = px[k];
            }            
            if(trace > 1){
                PRINTF("info(i=%i e=%i h=%i)\n",info[0],info[1],info[2]);
                PRINTF("px ");for(int q = 0; q < m; q++)PRINTF("%i ",px[q]);PRINTF("\n"); 
                PRINTF("qx ");for(int q = 0; q < n-m; q++)PRINTF("%i ",qx[q]);PRINTF("\n");
            }         

//           if(!is_subset(nv,which_vx,n - m,&qx[0])) continue;
//           if(!is_subset(n - m,&qx[0],nv,which_vx)) continue;
//             if(m < m0)if(!is_subset(m,&px[0],n_best_px,&best_px[0])) continue;
//             if(m < m0)if(!is_subset(n-m,qx,n_best_qx,&best_qx[0])) continue;
//           if(!is_subset(nv,which_vx,m,&px[0])) continue;

            if(!get_sub_QP(m,&px[0],1,0,qp,sub)){continue;}
       
// Process the unrestricted-projected vertex            

            // if(trace > 2){PRINTF("\nsub-projections of PROJ(unrestricted) = ");for(int q = 0; q < n; q++)PRINTF("%+.3e ",x1[q]);PRINTF("\n");} 
            // if(!get_sub_solution(qp,sub,m,&px[0],n-m,qx,&x1[0],x,MAX(0,trace - 1))){
                // if(trace > 2){PRINTF("projection out of bounds: f(x)=%+.3e x=",objective_QP(qp,x));for(int q = 0; q < n; q++)PRINTF("%+.3e ",x[q]);PRINTF("\n");}                 
              
            // }else{
// //                double fx = objective_QP(qp,x);
                // one = 1.0;zero = 0.0;inc = 1;
                // F77_CALL(dsymv)("U", &qp->n, &one, qp->H, &qp->n, x, &inc, &zero, &buf[0], &inc);
                // double fx = 0.5*F77_CALL(ddot)(&qp->n, x, &inc, &buf[0], &inc) + F77_CALL(ddot)(&qp->n, x, &inc, qp->b, &inc);
                
                // if(fx < best_fx){
                    // best_fx = fx;
                    // for(int j = 0; j < n; j++)solution[j] = x[j];                
// //                    for(int j = 0; j < m; j++)best_px[j] = px[j];n_best_px = m;                                    
// //                    if(qx){for(int j = 0; j < n-m; j++)best_qx[j] = qx[j];n_best_qx = n-m;}else n_best_qx = 0;                                    
                // }
                // C_n_m_FEASIBLE[m-1]++;
                // if(trace > 2){PRINTF("projection in bounds: f(x)=%+.3e x=",fx);for(int q = 0; q < n; q++)PRINTF("%+.3e ",x[q]);PRINTF("\n");}
            // }
            


            int nmHV = 1<<(n - m);
            for(int i = 0;i < nmHV;i++){
                next_hypercube_vertex(i+1,n - m,qp->xl[0],qp->xu[0],x);                
                for(int j = 0;j < n-m;j++)x0[qx[j]] = x[j];                           
                if(trace > 2){PRINTF("\nsub-projections of hypercube vertex: ");for(int q = 0; q < n; q++)PRINTF("%+.3e ",x0[q]);PRINTF("\n");}                                 
                if(!get_sub_solution(qp,sub,m,&px[0],n-m,qx,x0,x,MAX(0,trace - 1))){
                    if(trace > 2){
                        PRINTF("projection out of bounds: f(x)=%+.3e x=",objective_QP(qp,x));
                        for(int q = 0; q < n; q++)PRINTF("%+.3e ",x[q]);PRINTF("\n");
                    } 
                    continue;  
                }
                
//                double fx = objective_QP(qp,x);
                one = 1.0;zero = 0.0;inc = 1;
                F77_CALL(dsymv)("U", &qp->n, &one, qp->H, &qp->n, x, &inc, &zero, qp->dbuf, &inc FCONE);
                double fx = 0.5*F77_CALL(ddot)(&qp->n, x, &inc, qp->dbuf, &inc) + F77_CALL(ddot)(&qp->n, x, &inc, qp->b, &inc);
                if(fx < best_fx){
                    best_fx = fx;
                    for(int j = 0; j < n; j++)solution[j] = x[j];
//                    for(int j = 0; j < n; j++)best_x0[j] = x0[j];                    
                }
                C_n_m_FEASIBLE[m-1]++;
                if(trace > 2){
                    PRINTF("projection in bounds f(x)=%+.3e x=",fx);
                    for(int q = 0; q < n; q++)PRINTF("%+.3e ",x[q]);PRINTF("\n");
                } 
//                    break;
                
//                if(C_n_m_FEASIBLE[m-1]) break;
            }//end hypercube loop  
//            if(m > qp->rank)
//            if(C_n_m_FEASIBLE[m-1]) break;
            
        }// end all 'n0 choose m' combinations
//        for(int j = 0; j < n; j++)if(j < m) best_px[j] = px[j];
//        if(C_n_m_FEASIBLE[m-1]){free_QP(sub);return 1;}
    }// end C(n,1)..C(n,n) combinations       
    if(trace > 1)PRINTF("Hypercube vertices: \n ");
    for(int i = 0;i < nHV;i++){
        next_hypercube_vertex(i+1,n,qp->xl[0],qp->xu[0],x);        
//        double fx = objective_QP(qp,x);
        double one = 1.0,zero = 0.0;int inc = 1;
        F77_CALL(dsymv)("U", &qp->n, &one, qp->H, &qp->n, x, &inc, &zero, qp->dbuf, &inc FCONE);
        double fx = 0.5*F77_CALL(ddot)(&qp->n, x, &inc, qp->dbuf, &inc) + F77_CALL(ddot)(&qp->n, x, &inc, qp->b, &inc);
        if(trace > 1){PRINTF("f(x)=%+.3e x=",fx);for(int q = 0; q < n; q++)PRINTF("%+.3e ",x[q]);PRINTF("\n");} 
        if(fx < best_fx){
            best_fx = fx;
            for(int j = 0; j < n; j++)solution[j] = x[j];
        }
    }    
    if(trace){
        PRINTF("Number of Feasible points ");for(int i = 0;i < n;i++) PRINTF("C(%i,%i)=%i ",n,i+1,C_n_m_FEASIBLE[i]);PRINTF("\n");
        PRINTF("projected unrestricted estimate: f(x)=%+.3f x=",x1[n]);for(int q = 0; q < n; q++)PRINTF("%.3f ",x1[q]);PRINTF("\n");
        PRINTF("best estimate: f(x)=%+.3e x=",best_fx);for(int q = 0; q < n; q++)PRINTF("%+.3e ",solution[q]);PRINTF("\n");
    } 
    return 1;
}     

// int optimize_qp_exhaustive_xout(MNQD_QP* qp, double* solution,double* hypercube,int flag,int trace)
// {
    // int n = qp->n;
    // double* x0 = NULL;
    // double* x = solution;
    // double best_fx = _Inf;
    // double x1[n];
    // unsigned long nHV = 1<<n; // process 2^n vertices of hypercube  
    // int C_n_m_FEASIBLE[n];memset(&C_n_m_FEASIBLE[0],0,n * sizeof(int));
    // MNQD_QP sub;
    
    // if(trace) PRINTF("\n");    
    // if(trace > 1) print_QP(qp);

    // memset(&sub,0,sizeof(MNQD_QP));
    // if(!get_sub_QP(n,NULL,0,0,qp,&sub))return 0; //linear term in restricted obj. fn depends on 'x', set in get_sub_solution    

    // double zero = 0.0,one = -1.0;int inc = 1;
    // double buf[qp->n];
    // int best_px[n];int n_best_px;
    // int best_qx[n];int n_best_qx;       
// // buffers    
    // double _x[n];x = &_x[0];        
    // double _x0[n];x0 = &_x0[0]; 
    

    // int zx[n];for(int k = 0;k < n;k++)zx[k] = k;
    // int m0 = qp->rank;
    // for(int m = m0; 1 <= m; m--)
    // {
        // int cx[m]; // stores m-combination that will be passed to hypercube vertex iterator
        // int px[m]; // C-based index of m-combination, px = cx - 1
        // int _qx[MAX(1,n - m)]; // compliment of px, NULL if px is the full length 'n' combination
        // int* qx = n > m ? &_qx[0] : NULL; 

        // if(trace > 1)PRINTF("Combinations Cn(%i,%i)\n",n,m);  
        // int info[3] = {0,0,0};
        // while(iter_combn(&info[0],n,m,&cx[0])){
            // // subtract 1 for C-based indexing
            // for(int k = 0;k < m;k++)px[k] = cx[k]-1;             
            // if(qx){
                // for(int k = 0;k < m;k++)zx[px[k]] = -1;
                // int l = 0;for(int k = 0;k < n;k++)if(zx[k]!=-1) qx[l++] = zx[k];
                // for(int k = 0;k < m;k++)zx[px[k]] = px[k];
            // }            
            // if(trace > 1){
                // PRINTF("info(i=%i e=%i h=%i)\n",info[0],info[1],info[2]);
                // PRINTF("px ");for(int q = 0; q < m; q++)PRINTF("%i ",px[q]);PRINTF("\n"); 
                // PRINTF("qx ");for(int q = 0; q < n-m; q++)PRINTF("%i ",qx[q]);PRINTF("\n");
            // }         
            // if(!get_sub_QP(m,&px[0],1,0,qp,&sub)){continue;}

            // int nmHV = nHV;
// //            if(n > m) 
            // nmHV = 1<<(n - m);
            // for(int i = 0;i < nmHV;i++){
                // next_hypercube_vertex(i+1,n - m,qp->xl[0],qp->xu[0],x);                
                // for(int j = 0;j < n-m;j++)x0[qx[j]] = x[j];                           
            // // for(int i = 0;i < nHV;i++){
                // // next_hypercube_vertex(i+1,n,qp->xl[0],qp->xu[0],x0);


                // if(trace > 2){PRINTF("\nsub-projections of hypercube vertex: ");for(int q = 0; q < n; q++)PRINTF("%+.3e ",x0[q]);PRINTF("\n");}                                 
                // if(!get_sub_solution(qp,&sub,m,&px[0],n-m,qx,x0,x,MAX(0,trace - 1))){
                    // if(trace > 2){
                        // PRINTF("projection out of bounds: f(x)=%+.3e x=",objective_QP(qp,x));
                        // for(int q = 0; q < n; q++)PRINTF("%+.3e ",x[q]);PRINTF("\n");
                    // } 
                    // continue;  
                // }
                
// //                double fx = objective_QP(qp,x);
                // one = 1.0;zero = 0.0;inc = 1;
                // F77_CALL(dsymv)("U", &qp->n, &one, qp->H, &qp->n, x, &inc, &zero, &buf[0], &inc);
                // double fx = 0.5*F77_CALL(ddot)(&qp->n, x, &inc, &buf[0], &inc) + F77_CALL(ddot)(&qp->n, x, &inc, qp->b, &inc);
                // if(fx < best_fx){
                    // best_fx = fx;
                    // for(int j = 0; j < n; j++)solution[j] = x[j];
// //                    for(int j = 0; j < n; j++)best_x0[j] = x0[j];                    
                // }
            // }//end hypercube loop              
        // }// end all 'n choose m' combinations
    // }// end C(n,1)..C(n,n) combinations       
    // if(trace > 1)PRINTF("Hypercube vertices: \n ");
    // for(int i = 0;i < nHV;i++){
        // next_hypercube_vertex(i+1,n,qp->xl[0],qp->xu[0],x);        
// //        double fx = objective_QP(qp,x);
        // double buf[qp->n];double one = 1.0,zero = 0.0;int inc = 1;
        // F77_CALL(dsymv)("U", &qp->n, &one, qp->H, &qp->n, x, &inc, &zero, &buf[0], &inc);
        // double fx = 0.5*F77_CALL(ddot)(&qp->n, x, &inc, &buf[0], &inc) + F77_CALL(ddot)(&qp->n, x, &inc, qp->b, &inc);
        // if(trace > 1){PRINTF("f(x)=%+.3e x=",fx);for(int q = 0; q < n; q++)PRINTF("%+.3e ",x[q]);PRINTF("\n");} 
        // if(fx < best_fx){
            // best_fx = fx;
            // for(int j = 0; j < n; j++)solution[j] = x[j];
        // }
    // }    

    // free_QP(&sub);
    // return 1;
// }     

// void ui(int* _q,int* x){
    // unsigned int ii = (unsigned int)*_q;
    // int j = 0;
    // do{x[j++] = (int)(ii & 1);ii = ii >> 1;}while(ii);
// }
// void ull(int* _q,int* x){
    // unsigned long long ii = (unsigned long long)*_q;
    // int j = 0;
    // do{x[j++] = (int)(ii & 1);ii = ii >> 1;}while(ii);
// }
// void u(int* _q,int* x){
    // unsigned ii = (unsigned)*_q;
    // int j = 0;
    // do{x[j++] = (int)(ii & 1);ii = ii >> 1;}while(ii);
// }
// void i(int* _q,int* x){
    // int ii = *_q;
    // int j = 0;
    // do{x[j++] = (ii & 1);ii = ii >> 1;}while(ii);
// }
// library(aucm)
 // q <- as.integer(3);qq <- 0;if(q)qq <- integer(1 + log(q)/log(2))
// .C('i',q,qq)
// .C('ui',q,qq)
// .C('ull',q,qq)
// .C('u',q,qq)
