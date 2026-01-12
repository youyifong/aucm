#include <string.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <R.h>
#define _Inf DBL_MAX
#define intCEIL(X,Y) (1 + (((X) - 1) / (Y))) //ceil(X/Y)  
#define PRINTF Rprintf

///////////////////////////////////////////////////////////////////////////////
/////////////////// Can define or comment out the following ///////////////////
// #define DO_VLA 
// a very bad idea to define this for any n = n1 + n2 s.t. the size of 'x' in 
// double x[n1*n2] is greater than a few kilobytes, unless your system has
// large amount of memory dedicated to the STACK

// // // Kahan-summing for updating 'Ka' 
//#define DO_KAHAN

//#define MEASURE_ELAPSED_TIME
#ifdef MEASURE_ELAPSED_TIME
#include <time.h>
#endif
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
 
#define DIM_B 2

extern void get_machine_double_eps(double* eps);

typedef struct ptr_mnqd{
double **H; //to conveniently access either 'K' or for 'Q' depending on whether mem_efficient = 0/1
double *Ka,*df;
#ifdef DO_KAHAN
double* Ka_bits;
#endif
}PTR_MNQD2;
PTR_MNQD2 state2;

void free_minQuad_2(){
    if(state2.H){free(state2.H);state2.H = NULL;}
    if(state2.Ka){free(state2.Ka);state2.Ka = NULL;}
    if(state2.df){free(state2.df);state2.df = NULL;}
#ifdef DO_KAHAN
    if(state2.Ka_bits){free(state2.Ka_bits);state2.Ka_bits=NULL;}
#endif    
}


//output :B = (p : max{-df[p] | a[p] < C} , q : min{-df[q] | a[q] > 0})
//output :'epsilon'
void workingSet_2(int n1n2,double* a,double * df,double C,int* B,double* epsilon) {
    int _i = 0;
    int _j = 0;
    double M_C = -_Inf;
    double m_0 = +_Inf;
    for(int p = 0; p < n1n2; p++) {
        if(a[p] <   C) {
            if(-df[p] > M_C) {
                _i  = p;
                M_C = -df[p];
            }
        }
        if(a[p] > 0.0) {
            if(-df[p] < m_0) {
                _j  = p;
                m_0 = -df[p];
            }
        }
    }
    *B = _i;
    *(B+1) = _j;
    *epsilon = M_C - m_0;
}



void minQuad_2(
    //input
    int* _n1n2,        // n1n2 = n_diseased * n_non_diseased from R
    double* b,         // length n1n2
    double* _C,        // length 1
    double* _H,
    int* mem_efficient,//  1/0 <-> _H = K_nxn or _H = Q_n1n2xn1n2 
    int* _n1,
    int* _n2,
    int* maxit,        // max no. of iterations, >= 1
    double* tol,       // tolerance criterion
	double* machine_double_eps,	
    double* a0,        // initial vector of alpha's, length n1n2
    int* _trace,       // 0,1,2: print at each iteration (2), before returning to caller (1) no print(0)
    int* return_working_set, //save working-set at each iteration  up to 'maxit'?
    //output
    double* a,         // return final alpha of length 'n1n2' to caller
    double* v,         // return value of objective function to caller
    int* iter,         // return the final number of iteration to caller
    int* convergence,  // did alg. converge (Y/N) <-> (0/1)
    double * epsilon,  // return value in to termination condition to caller
    int* working_set   // save ((B_1),(B_2),..,(B_maxit)) B_k = (p,q), length = 2 * maxit
)
{
    int trace = *_trace;
    int p,q,q1,q2,done,ip,jp,iq,jq;
    int n1n2 = *_n1n2;
    int n1 = *_n1;
    int n2 = *_n2;
    int n = n1+n2;
    double C = *_C;

// working set
    int sum_vB = -1;
    int B[DIM_B]= {0,0};// current working set (_i,_j)
    double aB[DIM_B] = {0.0,0.0};
    double aB1[DIM_B] = {0.0,0.0};
    double daB[DIM_B] = {0.0,0.0};
    double aB0[DIM_B] = {0.0,0.0};
    double bB[DIM_B] = {0.0,0.0};
    // in general hypercube of dimension 2^DIM_B x DIM_B
    double B_CUBE[4][DIM_B] = {{0.0,0.0},{0.0,1.0},{1.0,1.0},{1.0,0.0}};    
    double det_Q_BB = 1.0;
    double Q_BB[DIM_B][DIM_B];
    double iQ_BB[DIM_B][DIM_B]; // inverse of 'Q[B,B]'
   
    
    #ifdef MEASURE_ELAPSED_TIME
    clock_t t0_init = clock();
    clock_t t1_init,t0_loop,t1_loop,t0_end,t1_end;
    #endif    

	double MACHINE_DOUBLE_EPS = *machine_double_eps;
	if(MACHINE_DOUBLE_EPS + 1.0 != 1.0)get_machine_double_eps(&MACHINE_DOUBLE_EPS);
    
    memset(&state2,0,sizeof(state2));// all pointers set to NULL
    double* df = (double*)calloc(n1n2 , sizeof(double));if(!df){free_minQuad_2();return;}
    state2.df = df;
    double* Ka = (double*)calloc(n1n2 , sizeof(double));if(!Ka){free_minQuad_2();return;}
    state2.Ka = Ka;
    #ifdef DO_KAHAN
    double* Ka_bits = (double*)calloc(n1n2 , sizeof(double));if(!Ka_bits){free_minQuad_2();return;}
    state2.Ka_bits = Ka_bits;
    #endif    
    
#ifdef DO_VLA
    double (*K)[n]    = (double (*)[n])_H;	
    double (*Q)[n1n2] = (double (*)[n1n2])_H;	

#else
    double** Q = NULL;double** K = NULL;
    if(!(*mem_efficient)) {
        if(n1n2 <= 0){PRINTF("Cannot allocate size n1n2 = 0 memory.\n");done = 1;}
        Q = (double**) calloc(n1n2,sizeof(double*));
        if(!Q){
				//PRINTF("minQuad: Unable to allocate %lld(bytes)\n",n1n2*sizeof(double*));
				//free_minQuad_2();return;}
        for(p = 0; p < n1n2; p++)Q[p] = _H + p*n1n2;
        state2.H = Q;
    }else{
        K = (double**) calloc(n,sizeof(double*));if(!K){free_minQuad_2();return;}
        for(p = 0; p < n; p++)K[p] = _H + p*n;
        state2.H = K;
    }

#endif
    
// Initialization
    *convergence = 1;
    *iter = 0;
    *epsilon = 2.0 * (*tol);
	for(p = 0; p < n1n2; p++)a[p] = a0[p];//memcpy(a,a0,n1n2 * sizeof(double));
//  Rprintf("init: value(%.1f) eps(%.1f) reltol(%.1f) convergence(%p)\n1n2",*v,*epsilon,*reltol,*convergence);
    if(C != 1.0){
        for(int i = 0; i < 4; i++) 
            for(int j = 0; j < 2; j++) 
                B_CUBE[i][j] *= C;
    } 
// Caching Ka
    if(!(*mem_efficient)) {
        for(p = 0; p < n1n2; p++) {
            Ka[p] = 0.0;
            for(q = 0; q < n1n2; q++){
            #ifdef DO_KAHAN
                double volatile y = Q[p][q] * a[q] - Ka_bits[p];
                double volatile t = Ka[p] + y;
                Ka_bits[p] = t - Ka[p] - y;
                Ka[p] = t;                                    
            #else
                  Ka[p] += Q[p][q] * a[q];
            #endif      
            }
        }
    }else {
        for(p = 0; p < n1n2; p++) {
            Ka[p] = 0.0;
#ifdef intCEIL
            ip=intCEIL(p+1,n2);
#else
            ip=(int)ceil(((double)p+1.0)/(double)n2)
#endif
            jp=(p+1)-(ip-1)*n2;  
            for(q = 0; q < n1n2; q++) {
//                  ip=ceil((p+1.0)/n2); jp=(p+1)-(ip-1)*n2; iq=ceil((q+1.0)/n2); jq=(q+1)-(iq-1)*n2;                     
#ifdef intCEIL
                iq=intCEIL(q+1,n2); 
#else
                iq=(int)ceil(((double)q+1.0)/(double)n2);
#endif
                jq=(q+1)-(iq-1)*n2; 

            #ifdef DO_KAHAN
               double volatile y = a[q] * (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]) - Ka_bits[p];
               double volatile t = Ka[p] + y;                
               Ka_bits[p] = t - Ka[p] - y;
               Ka[p] = t;    
            #else   
                Ka[p] += a[q] * (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]);
            #endif
            }
        }
    }



    #ifdef MEASURE_ELAPSED_TIME
    t1_init = clock();
    t0_loop = t1_init;
    #endif    
    done = 0;
    while(!done) {
        (*iter)++;
        for(p = 0; p < n1n2; p++) df[p] = 2.0 * Ka[p] - b[p];
        
        
        workingSet_2(n1n2,a,df,C,&(B[0]),epsilon); //B,epsilon set
        // int _i = 0;
        // int _j = 0;
        // double M_C = -_Inf;
        // double m_0 = +_Inf;
        // for(p = 0; p < n1n2; p++) {
            // if(a[p] <   C) {
                // if(-df[p] > M_C) {
                    // _i  = p;
                    // M_C = -df[p];
                // }
            // }
            // if(a[p] > 0.0) {
                // if(-df[p] < m_0) {
                    // _j  = p;
                    // m_0 = -df[p];
                // }
            // }
        // }
        // B[0] = _i;
        // B[1] = _j;
        // *epsilon = M_C - m_0;
        
        if(*return_working_set) {
            for(q1 = 0; q1 < DIM_B; q1++) working_set[(*iter - 1)*DIM_B + q1] = B[q1] + 1; //'+1' for 1-based indexing
        }
        if(*epsilon <= *tol) {
            done = 1;
            *convergence = 0;
            if(trace > 1)trace = 1;
            break;
        }

////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////
        // get Q[B][B]
        if(!(*mem_efficient)) {
            for(q1 = 0; q1 < DIM_B; q1++)
                for(q2 = 0; q2 < DIM_B; q2++)
                    Q_BB[q1][q2] = Q[B[q1]][B[q2]];
        }else{        
            for(q1 = 0; q1 < DIM_B; q1++) {
                p = B[q1];
    #ifdef intCEIL
                ip=intCEIL(p+1,n2);
    #else
                ip=(int)ceil(((double)p+1.0)/(double)n2)
    #endif
                jp=(p+1)-(ip-1)*n2;  
                for(q2 = 0; q2 < DIM_B; q2++) {
                    q=B[q2];
    //              ip=ceil((p+1.0)/n2); jp=(p+1)-(ip-1)*n2; iq=ceil((q+1.0)/n2); jq=(q+1)-(iq-1)*n2;                     
    #ifdef intCEIL
                    iq=intCEIL(q+1,n2); 
    #else
                    iq=(int)ceil(((double)q+1.0)/(double)n2);
    #endif
                    jq=(q+1)-(iq-1)*n2; 
                    Q_BB[q1][q2] = (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]);
                }
            }
        }
        
        // b_B = b[B]- 2 K_BN %*% alpha_N
        for(q1 = 0; q1 < DIM_B; q1++) {
            bB[q1] = 0.0;
            for(q2 = 0; q2 < DIM_B; q2++) {
                if(a[B[q2]] != 0.0) {
                    bB[q1] +=  Q_BB[q1][q2] * a[B[q2]];
                }
            }
            bB[q1] -= Ka[B[q1]];
            bB[q1] *= 2.0;
            bB[q1] += b[B[q1]];
        }

////////////////////////////////////////////////////////////////////////////////////////////
        // new candidate 'a_B', alpha for working set B
        // a_B is solution beta_hat in Q[B][B] beta = bB
        // bB = b_B* - 2 * Q[B,-B] %*% alpha[-B]
        // X = 2 * Q[B,B]
        // aB 2x2 solution of Xbeta = y, aB = 0.5Q_BB^(-1)bB
        // iQ_BB = Q_BB^(-1):
        det_Q_BB =  Q_BB[0][0] * Q_BB[1][1] - Q_BB[0][1] * Q_BB[1][0];
        if(det_Q_BB == 0.0) {
            Rprintf("Error: det(Q_BB) = 0\n");
            break;
        }
        iQ_BB[0][1] = - Q_BB[1][0];
        iQ_BB[1][0] = - Q_BB[0][1];
        iQ_BB[0][0] = + Q_BB[1][1];
        iQ_BB[1][1] = + Q_BB[0][0];
        for(q1 = 0; q1 < DIM_B; q1++) {
            aB[q1] = 0.0;
            for(q2 = 0; q2 < DIM_B; q2++) aB[q1] += iQ_BB[q1][q2]*bB[q2];
            aB[q1] *= 0.5 / det_Q_BB;
        }
////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////
        
        int COMPUTE_RESTRICTED_OBJECTIVE_FUNCTION = 1;
// 8 'squares' in the plane violating square 5 specified by {(0,0)',(0,C)',(C,0)',(C,C)'}
// 1 2 3
// 4 5 6
// 7 8 9
//for (DIM_B = 2,sum_vB = 2), the stmt "if((0.0 <= d_nq) && (d_nq <= C))" is true for at most one of them, hence the 'break' stmt.
        int vB[DIM_B]; //logical indicates if a coordinate of aB 'violates', p.e. outside of [0,C]
        sum_vB = 0; //number of violating coordinates in alpha_B,{0,1,2}
        double m_q = NAN;
        double d_nq = NAN;
        for(q1 = 0; q1 < DIM_B; q1++) {
            vB[q1] = (int) ((aB[q1] < 0.0) | (C < aB[q1]));
            sum_vB += vB[q1]; //{0,1,2}
        }


        if(sum_vB == 0) { //square 5, no coordinate of alpha_B violates
            for(q1 = 0; q1 < DIM_B; q1++) aB1[q1] = aB[q1];
            COMPUTE_RESTRICTED_OBJECTIVE_FUNCTION = 0;
        } else { // at least 1 coordinate of alpha_B violates
            for(q1 = 0; q1 < DIM_B; q1++) {
                if(vB[q1]) { //violating coordinate
                    int nq = !q1;
                    m_q = aB[q1];
                    if(m_q > C) m_q = C;
                    else if(m_q < 0.0) m_q = 0.0;
                    d_nq = (0.5*bB[nq] - m_q * Q_BB[nq][q1]) / Q_BB[nq][nq];
                    if((0.0 <= d_nq) && (d_nq <= C)) {
                        aB1[nq] = d_nq;
                        aB1[q1]  = m_q;
                        COMPUTE_RESTRICTED_OBJECTIVE_FUNCTION = 0;
                        break; //Shuxin: cannot get here at DIM_B(=2) times, only 0 or 1 time
                    }
                }
            }
        }

        if(COMPUTE_RESTRICTED_OBJECTIVE_FUNCTION){
        // reduced objective function f_B(x) evaluated on a square
        // f_B(x) = x' Q_BB x - b_B' x
        // find alpha_min = x in { (0,0)',(0,C)',(C,0)',(C,C)' } s.t. f_B(x) is minimum

            double fx_min = _Inf;            
            for(q = 0;q < 4;q++){
                double fx = 0.0;
                double* x = (double*)&(B_CUBE[q][0]);
                for(int i = 0; i < DIM_B; i++) {
                    if(fabs(x[i]) <= MACHINE_DOUBLE_EPS) continue;
                    for(int j = 0; j < DIM_B; j++) {
                        if(fabs(x[j]) <= MACHINE_DOUBLE_EPS) continue;
                        fx   +=  x[i] * Q_BB[i][j] * x[j];
                    }
                }        
                for(int i = 0; i < DIM_B; i++) fx -= bB[i] * x[i];
                if(fx < fx_min){
                    fx_min = fx;
                    for(int j = 0; j < DIM_B; j++)aB1[j] = x[j];
                }
            }
        }
        
// update alpha and Ka = Q %*% alpha at B
        for(q1 = 0; q1 < DIM_B; q1++) {
            aB0[q1]  = a[B[q1]];
            a[B[q1]] = aB1[q1];
            daB[q1]  = aB1[q1] - aB0[q1];
        }


        
        if(!(*mem_efficient)) {
            for(q1 = 0; q1 < DIM_B; q1++){
			    if(daB[q1] == 0.0) continue;
                q = B[q1];
                for(p=0; p < n1n2; p++){
                #ifdef DO_KAHAN
                    double volatile y = Q[p][q] * daB[q1] - Ka_bits[p];
                    double volatile t = Ka[p] + y;
                    Ka_bits[p] = t - Ka[p] - y;
                    Ka[p] = t;                                      
                #else
                    Ka[p] += Q[p][q] * daB[q1];
                #endif                      
                }
            }
        }else{
        
            for(q1 = 0; q1 < DIM_B; q1++) {
					if(daB[q1] == 0.0) continue;
                    q = B[q1];
#ifdef intCEIL
                    iq=intCEIL(q+1,n2); 
#else
                    iq=(int)ceil(((double)q+1.0)/(double)n2);
#endif
                    jq=(q+1)-(iq-1)*n2;

                    for(p=0; p < n1n2; p++) {
//                            ip=ceil((p+1.0)/n2); jp=(p+1)-(ip-1)*n2; iq=ceil((q+1.0)/n2); jq=(q+1)-(iq-1)*n2;                     
#ifdef intCEIL
                        ip=intCEIL(p+1,n2);
#else
                        ip=(int)ceil(((double)p+1.0)/(double)n2)
#endif
                        jp=(p+1)-(ip-1)*n2;  
                        
                    #ifdef DO_KAHAN
                        double volatile y = -Ka_bits[p] + daB[q1] * (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]);
                        double volatile t = Ka[p] + y;
                        Ka_bits[p] = t - Ka[p] - y;
                        Ka[p] = t;
                    #else
                        Ka[p] += daB[q1] * (K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1]);
                    #endif    
                    }
            }
        }
        
        if(*iter == *maxit) {
            done = 1;
            *convergence = 1;
            if(trace > 1)trace = 1;
            break;
        }
        if(trace > 1) PRINTF("i(%i) #V(%i) B(%i,%i) e(%+.2e) M(-df_C)(%+.3e) m(-df_0)(%+.3e) det(%.10f) da[B](%+.3e,%+.3e) a[B](%+.3e,%+.e) aB1(%+.3e,%+.3e) aB(%+.3e,%+.3e) bB(%+.3e,%+.3e)\n",
                                  *iter,sum_vB,B[0]+1,B[1]+1,*epsilon,-df[B[0]],-df[B[1]],det_Q_BB,daB[0],daB[1],a[B[0]],a[B[1]],aB1[0],aB1[1],aB[0],aB[1],bB[0],bB[1]);

    } // end while loop
    #ifdef MEASURE_ELAPSED_TIME
    t1_loop = clock();t0_end = t1_loop;
    #endif    

// objective function evaluated at parameter
    *v = 0.0;for(q = 0; q < n1n2; q++)*v += a[q] * (Ka[q] - b[q]);

    if(trace == 1) PRINTF("i(%i) #V(%i) B(%i,%i) e(%+.2e) M(-df_C)(%+.3e) m(-df_0)(%+.3e) det(%.10f) da[B](%+.3e,%+.3e) a[B](%+.3e,%+.e) aB1(%+.3e,%+.3e) aB(%+.3e,%+.3e) bB(%+.3e,%+.3e)\n",
                               *iter,sum_vB,B[0]+1,B[1]+1,*epsilon,-df[B[0]],-df[B[1]],det_Q_BB,daB[0],daB[1],a[B[0]],a[B[1]],aB1[0],aB1[1],aB[0],aB[1],bB[0],bB[1]);


    free_minQuad_2();

    #ifdef MEASURE_ELAPSED_TIME
    t1_end = clock();    
    PRINTF("minQuad elapsed initialization : %f seconds\n", (double)(t1_init - t0_init) / CLOCKS_PER_SEC);
    PRINTF("minQuad elapsed loop           : %f seconds\n", (double)(t1_loop - t0_loop) / CLOCKS_PER_SEC);
    PRINTF("minQuad elapsed exit           : %f seconds\n", (double)(t1_end - t0_end) / CLOCKS_PER_SEC); 
    #endif    
}


