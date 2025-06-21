#include "working_set.h"
#include "sort.h"
#include "random.h"

#include <R.h>
#define PRINTF Rprintf

#include <math.h>
#include <float.h>
#define _Inf DBL_MAX


#define	MAX(A,B)	((A) > (B) ? (A) : (B))
#define	MIN(A,B)	((A) < (B) ? (A) : (B))


// output nV - number of violators
// output v  - violators
// v[p] = 1 if KKT conditions are not met, 0 otherwise
int getKKTViolators(int n, double* a, double* df, double* C,int* v) 
{
    int nV = 0;
    for(int p = 0; p < n; p++) {      
        v[p] = 0;    
        if(((a[p] < C[p]) && (df[p] < 0.0)) || ((a[p] > 0.0) && (df[p] > 0.0))){ 
            v[p] = 1; 
            nV++;
        }
    }  
    return(nV);
}

// v[i] = 0 if pth entry (a[p],df[p]) does not violate KKT conditions
// v[i] = 1 if pth entry is s.t. (0 < a[p],0 < df[p]) 
// v[i] = 2 if pth entry is s.t. (C > a[p],0 > df[p]) 
// output nV - number of violators
// output nV0 - number of violators for 0 < a[p]
// output nVC - number of violators for C > a[p]
// output v   - violators
int getKKTViolators_0C(int n, double* a, double* df, double* C,int* nV0,int* nVC,int* v) 
//int getKKTViolators_0C(int n, double* a, double* df, double C,int* nV0,int* nVC,int* v,int* which_v,int* which_v0,int* which_vc) 
{

    *nV0 = 0;*nVC = 0;
    for(int p = 0; p < n; p++) {      
       v[p] = 0;     
       if((a[p] < C[p]) && (df[p] < 0.0)){
        v[p] = 2; 
//        if(which_v)which_v[*nV0 + *nVC] = p;
//        if(which_vc)which_vc[*nVC] = p;
        (*nVC)++;
       }else if((a[p] > 0.0) && (df[p] > 0.0)){
        v[p] = 1; 
//        if(which_v)which_v[*nV0 + *nVC] = p;
//        if(which_v0)which_v0[*nV0] = p;        
        (*nV0)++;
       }
    }  
    return(*nVC + *nV0);
}

    
//return epsilon according to termination criterion: df_max - df_min
//where (df_max,df_min) = (max(p){-df[p] | a[p] < C} , min(q){-df[q] | a[q] > 0})
double getEpsilonKKT(int n1n2,double* a,double* df, double* C){    
    double M_C = -_Inf;
    double m_0 = +_Inf;
    for(int p = 0; p < n1n2; p++) {
        if(a[p] <   C[p]){
            if(-df[p] > M_C) M_C = -df[p];
        }
        if(a[p] > 0.0){ 
            if(-df[p] < m_0) m_0 = -df[p];
        }
    }
    return(M_C - m_0);
}


//if one of V0 or V1 is empty, we return largest abs(df) of the other
double getEpsilonKKTV(int n1n2,double* a,double* df,double* C){    
    double M_C = -_Inf;
    double m_0 = +_Inf;
    
    for(int p = 0; p < n1n2; p++) 
        if((a[p] < C[p]) && (df[p] < 0.0) && (M_C < -df[p])) M_C = -df[p];
        
    for(int p = 0; p < n1n2; p++) 
        if((a[p] > 0.0) && (df[p] > 0.0) && (m_0 > -df[p])) m_0 = -df[p];
        
    if(M_C == -_Inf) M_C = 0.0;
    if(m_0 ==  _Inf) m_0 = 0.0;
 
    return(M_C - m_0);
}






// randomly select WithOut Replacement(WOR) 2 from violators (v1,v2) from V = V1 union V2
// v[p] = 1 if {p: 0 < a[p],0 < df[p]} U {p: C > a[p],0 > df[p]}
// 'rv' - r(random) v(1 set of violators)
int getWorkingSet_rv(int nB,int* B,int nV0,int nVC,int n1n2,int* v,int* which_v,int* buf) {
    int p;    
    int nV = nV0 + nVC;
//  PRINTF("Random working set selection n(%i) nV(%i) nVC(%i) nV0(%i): ",n1n2,nV,nVC,nV0);
    if(nV == 0) return 0;
    if(nV == 1){for(int p = 0; p < n1n2; p++)if(v[p]){for(int q = 0; q < nB; q++)B[q] = p;return nV;}};
    if(nV <= nB){
//        PRINTF("Warning: Number of violators(%i) is less than dimension of working set(%i).\n",nV,nB);
        int q = 0;for(int p = 0; p < n1n2; p++)if(v[p])B[q++] = p;
        return nV;        
    }          
    if(nV == 2){
        int i = 0;
        for(int p = 0; p < n1n2; p++) if(v[p]) B[i++] = p;   
        return 2;   
    }    
    int nC = 0,n0 = 0;
    if(nVC == 0){
        n0 = nB;
    }else if(nV0 == 0){
        nC = nB;
    }else{
        int half = nB / 2;
        nC = MIN(half,nVC);
        n0 = MIN(half,nV0);
        if(nC < half) n0 = nB - nC;else if(n0 < half) nC = nB - n0;
   }       
//        PRINTF("n0 %i nC %i nV0 %i nVC %i nB %i\n",n0,nC,nV0,nVC,nB);

    int* V = which_v;//  'true' index of violators, R's 'which()'       
    int q = 0;for(p = 0; p < n1n2; p++) if(v[p])V[q++] = p;
    SampleNoReplace(nB, nV, B, buf);
    for(p = 0; p < nB; p++)B[p] = V[B[p] - 1];   // -1 for C-indexing   
    return nB;
}

int getWorkingSet_rv2(int nB,int* B,int nV0,int nVC,int n1n2,int* v,int* which_v,int* perm) {
//int getWorkingSet_rv2(int n1n2,int* v,int nV0,int nVC,int nB,int* B) {

    int p;    
    int nV = nV0 + nVC;
//  PRINTF("Random working set selection n(%i) nV(%i) nVC(%i) nV0(%i): ",n1n2,nV,nVC,nV0);
    if(nV == 0) return 0;
    if(nV == 1){for(int p = 0; p < n1n2; p++)if(v[p]){for(int q = 0; q < nB; q++)B[q] = p;return nV;}};
    if(nV <= nB){
//        PRINTF("Warning: Number of violators(%i) is less than dimension of working set(%i).\n",nV,nB);
        int q = 0;for(int p = 0; p < n1n2; p++)if(v[p])B[q++] = p;
        return nV;        
    }          
    if(nV == 2){
        int i = 0;
        for(int p = 0; p < n1n2; p++) if(v[p]) B[i++] = p;   
        return 2;   
    }    
    int nC = 0,n0 = 0;
    if(nVC == 0){
        n0 = nB;
    }else if(nV0 == 0){
        nC = nB;
    }else{
        int half = nB / 2;
        nC = MIN(half,nVC);
        n0 = MIN(half,nV0);
        if(nC < half) n0 = nB - nC;else if(n0 < half) nC = nB - n0;        
    }       
//        PRINTF("n0 %i nC %i nV0 %i nVC %i nB %i\n",n0,nC,nV0,nVC,nB);

//    int perm[MAX(nV0,nVC)];
//    int V[MAX(nV0,nVC)];
    int* V = which_v;
      
    if(nC){  
        int q = 0;for(p = 0; p < n1n2; p++)if(v[p] == 2)V[q++] = p;
        SampleNoReplace(nC, nVC, B, perm);
        for(p = 0;p < nC;p++)B[p] = V[B[p] - 1];    
    }
    if(n0){
        int q = 0;for(p = 0; p < n1n2; p++)if(v[p] == 1)V[q++] = p;
        SampleNoReplace(n0, nV0, B+nC, perm);
        for(p = 0;p < n0;p++)B[p + nC] = V[B[p + nC] - 1];  
    }
    return n0 + nC;
}



// randomly WOR select 2 from violators (v1,v2) from V = V1 U V2 with weights proportional to |df| from
// v[p] > 0 if {p: 0 < a[p],0 < df[p]} U {p: C > a[p],0 > df[p]}
// 'rvwg' - r(random) w(weighted) v(1 set of violators)
int getWorkingSet_rvwg(int nB,int* B,int nV0,int nVC,int n1n2,double* df,int* v,int* which_v,int* perm,double* prob) {

    int p;    
    int nV = nV0 + nVC;
//  PRINTF("Random working set selection n(%i) nV(%i) nVC(%i) nV0(%i): ",n1n2,nV,nVC,nV0);
    if(nV == 0) return 0;
    if(nV == 1){for(int p = 0; p < n1n2; p++)if(v[p]){for(int q = 0; q < nB; q++)B[q] = p;return nV;}};
    if(nV <= nB){
//        PRINTF("Warning: Number of violators(%i) is less than dimension of working set(%i).\n",nV,nB);
        int q = 0;for(int p = 0; p < n1n2; p++)if(v[p])B[q++] = p;
        return nV;        
    }          
    if(nV == 2){
        int i = 0;
        for(int p = 0; p < n1n2; p++) if(v[p]) B[i++] = p;   
        return 2;   
    }    
   
//  'true' index of violators, R's 'which()'       
//    int V[nV];int* _V = &(V[0]); 
    int* V = which_v;
    int q = 0;for(p = 0; p < n1n2; p++) if(v[p]){V[q++] = p;}
    double s = 0.0; 
    for(p = 0; p < nV; p++)prob[p] = fabs(df[V[p]]);
  
    for(p = 0; p < nV; p++)s += prob[p];
    if(s == 0.0){
        for(p = 0; p < nV; p++)prob[p] = 1.0 / (double)nV;
    }else{
        for(p = 0; p < nV; p++)prob[p] /= s;
    }       
    ProbSampleNoReplace(nV,prob,perm,nB,B);
    for(p = 0; p < nB; p++)B[p] = V[B[p] - 1];    
    return nB;
}

// randomly (WOR) select violators (v1,v2,...) from (V1,V2) with weights based on |df|
// 'rv2wg' - r(random) w(weighted) v2(2 sets of violators)
// v[p] > 0 if V0[p] or VC[p]
// VC: v[p] == 2 if {-df[p] > 0,a[p] < C} 
// V0: v[p] == 1 if {-df[p] < 0,a[p] > 0} 
int getWorkingSet_rv2wg(int nB,int* B,int nV0,int nVC,int n1n2,double* df,int* v,int* which_v,int* perm,double* prob) {

    int p;    
    double s;
    int nV = nV0 + nVC;
//  PRINTF("Random working set selection n(%i) nV(%i) nVC(%i) nV0(%i): ",n1n2,nV,nVC,nV0);
    if(nV == 0) return 0;
    if(nV == 1){for(int p = 0; p < n1n2; p++)if(v[p]){for(int q = 0; q < nB; q++)B[q] = p;return nV;}};
    if(nV <= nB){
//        PRINTF("Warning: Number of violators(%i) is less than dimension of working set(%i).\n",nV,nB);
        int q = 0;for(int p = 0; p < n1n2; p++)if(v[p])B[q++] = p;
        return nV;        
    }          
    if(nV == 2){
        int i = 0;
        for(int p = 0; p < n1n2; p++) if(v[p]) B[i++] = p;   
        return 2;   
    }    
    int nC = 0,n0 = 0;
    if(nVC == 0){
        n0 = nB;
    }else if(nV0 == 0){
        nC = nB;
    }else{
        int half = nB / 2;
        nC = MIN(half,nVC);
        n0 = MIN(half,nV0);
        if(nC < half) n0 = nB - nC;else if(n0 < half) nC = nB - n0;        
    }   
    if(n0 + nC > nB){
		PRINTF("Logical error in getWorkingSet_rv2wg: n0(%i) + nC(%i) > nB(%i)\n",n0,nC,nB);
		return 0;
	}    
//   PRINTF("n0 %i nC %i nV0 %i nVC %i nB %i\n",n0,nC,nV0,nVC,nB);

    int* V = which_v;

    if(nC){  
//        int perm[nVC];double pC[nVC];int vC[nVC]; //  'true' index of violators, R's 'which()'  
        int q = 0;for(p = 0; p < n1n2; p++)if(v[p] == 2)V[q++] = p;

        s = 0.0;    
        for(p = 0; p < nVC; p++){
            prob[p] = MAX(0.0,-df[V[p]]); //negative should never occur, but just int case 
            s += prob[p];
        }
        if(s == 0.0)
            for(p = 0; p < nVC; p++)prob[p] = 1.0 / (double)nVC;
        else
            for(p = 0; p < nVC; p++)prob[p] /= s;
        ProbSampleNoReplace(nVC,prob,perm,nC,B);
		
        // for(p = 0;p < nC;p++){
			// if((p > nB-1) | (B[p] < 0) | (n1n2-1 < B[p])){
				// PRINTF("nV0(%i) nVC(%i) n0(%i) nC(%i) p=%i B[p]=%i V[B[p]-1] = V[%i] = %i\n",nV0,nVC,n0,nC,p,B[p], B[p]-1,V[B[p] - 1]);
				// return 0;
			// }
		// }
		
        for(p = 0;p < nC;p++)B[p] = V[B[p] - 1];    
    }
    
    if(n0){
//        int perm0[nV0];double p0[nV0];int v0[nV0];//  'true' index of violators
        int q = 0;for(p = 0; p < n1n2; p++)if(v[p] == 1)V[q++] = p;   
        s = 0.0;    
        for(p = 0; p < nV0; p++){
            prob[p] = MAX(df[V[p]],0.0); 
            s += prob[p];
        }
        if(s == 0.0) 
            for(p = 0; p < nV0; p++)prob[p] = 1.0 / (double)nV0;
        else 
            for(p = 0; p < nV0; p++)prob[p] /= s;  
			
        ProbSampleNoReplace(nV0,prob,perm,n0,B+nC); //B[0] was set from 'vC'

        // for(p = 0;p < n0;p++){
			// if((p+nC > nB-1) | (B[p+nC] < 0) | (n1n2-1 < B[p+nC])){
				// PRINTF("nV0(%i) nVC(%i) n0(%i) nC(%i) p+nC=%i B[p+nC]=%i V[B[p+nC]-1] = V[%i] = %i\n",nV0,nVC,n0,nC,p+nC,B[p+nC], B[p+nC] - 1,V[B[p+nC] - 1]);
				// return 0;
			// }
		// }		
//        for(p = 0;p < n0;p++)PRINTF("n0=%i p+nC=%i B[p+nC]=%i V[B[p+nC]-1] = V[%i] = %i\n",n0,p+nC,B[p+nC], B[p+nC] - 1,V[B[p+nC] - 1]);    
        for(p = 0;p < n0;p++)B[p+nC] = V[B[p+nC] - 1];    
    }
    return n0 + nC;
}


//output :B = (p,q) s.t. |df[p]| > |df[q] > ... among all violators
// if MIN(nB,nV) = 2 : B = (p,q) s.t. |df[p]| > |df[q] > all else) 
// if MIN(nB,nV) > 2 : B = (p1,p2,...) s.t. |df[p1]| > |df[p2] > among violators in V1 and V2
int getWorkingSet_v(int nB,int* B,int nV0,int nVC,int n1n2,double* a,double * df,int* v,int* which_v,int* perm,double* dbuf) {


    int nV = nV0 + nVC;
//  PRINTF("Random working set selection n(%i) nV(%i) nVC(%i) nV0(%i): ",n1n2,nV,nVC,nV0);
    if(nV == 0) return 0;
    if(nV == 1){for(int p = 0; p < n1n2; p++)if(v[p]){for(int q = 0; q < nB; q++)B[q] = p;return nV;}};
    if(nV <= nB){
//        PRINTF("Warning: Number of violators(%i) is less than dimension of working set(%i).\n",nV,nB);
        int q = 0;for(int p = 0; p < n1n2; p++)if(v[p])B[q++] = p;
        return nV;        
    }          
    if(nV == 2){
        int i = 0;
        for(int p = 0; p < n1n2; p++) if(v[p]) B[i++] = p;   
        return 2;   
    }    
    if(nB == 2){
        double M[2]= {-_Inf,-_Inf};
        for(int p = 0; p < n1n2; p++) {  
            if(!v[p]) continue;   
            double d = fabs(df[p]);
            if(d > M[0]){
                M[0] = d;
                B[0] = p;
            }else if((M[1] < d) && (d < M[0])){
                M[1] = d;
                B[1] = p;
            }   
        }   
        return 2;
   }

    int* V = which_v;
    int q = 0;for(int p = 0; p < n1n2; p++)if(v[p])V[q++] = p;
    for(int p = 0; p < nV; p++)dbuf[p]=fabs(df[V[p]]);  
    for(int p = 0; p < nV; p++)perm[p]=p+1;  
    revsort(&dbuf[0], &perm[0], nV);
    for(int p = 0; p < nB; p++)B[p] = V[perm[p] - 1];
    return nB;
       
}


  
//output :
// if MIN(nB,nV) = 2 : B = (p,q) = (p : max{-df[p] | a[p] < C} , q : min{-df[q] | a[q] > 0}) 
// if MIN(nB,nV) > 2 : B = 
int getWorkingSet_greedy_0C(int nB,int* B,int n1n2,double* a,double * df,double* C,int* v,int* which_v,int* perm,double* x) {

    if(nB < 1) return 0;

    if(nB <= 2){
        double M_C = -_Inf;
        double m_0 = +_Inf;  
        B[0] = -1;
        B[1] = -1;
        for(int p = 0; p < n1n2; p++) {     
            if(a[p] < C[p]) {
                if(-df[p] > M_C) {
                    B[0]  = p;
                    M_C   = -df[p];
                }
            }
            if(a[p] > 0.0) {
                if(-df[p] < m_0) {
                    B[1]  = p;
                    m_0 = -df[p];
                }
            }
        }    
        if(nB == 1){if(B[0] == -1)B[0] = B[1];}
        if((B[0] == -1) && (B[1] == -1)) return 0;
        if((B[0] == -1) || (B[1] == -1)) return 1;
        return nB;    
    }
    
    int p;
//  double x[n1n2];int perm[n1n2];int V[n1n2];
    int* V = which_v;
    
    int nVC = 0,nV0 = 0;
    int n0 = nB / 2;
    int nC = nB - n0;
    
   
    int q = 0;
    for(p = 0; p < n1n2; p++){ 
        if(a[p] < C[p]){
            V[q++] = p;
            nVC++;
        }
    }            
    
    for(p = 0; p < nVC; p++)x[p] = -df[V[p]];
    for(p = 0; p < nVC; p++)perm[p] = p+1;  
    revsort(&x[0], &perm[0], nVC); 
    for(p = 0; p < nC; p++)B[p] = V[perm[p] - 1];


    q = 0;
    for(p = 0; p < n1n2; p++){ 
        if(a[p] > 0.0){
            V[q++] = p;
            nV0++;
        }
    }
    for(p = 0; p < nV0; p++)x[p] = df[V[p]];
    for(p = 0; p < nV0; p++)perm[p] = p+1;  
    revsort(&x[0], &perm[0], nV0);
    for(p = 0; p < n0; p++)B[p+nC] = V[perm[p]-1];
    return n0 + nC;
      
}


//output :
// if MIN(nB,nV) = 2 : B = (p,q) = (p : max{-df[p] | a[p] < C} , q : min{-df[q] | a[q] > 0}) 
// if MIN(nB,nV) > 2 : B = most extreme violators from V1 and V2
int getWorkingSet_v2(int nB,int* B,int nV0,int nVC,int n1n2,double* a,double * df,int* v,int* which_v,int* perm,double* x) {

           
    int nV = nV0 + nVC;
    if(nV == 0) return 0;
    if(nV == 1){for(int p = 0; p < n1n2; p++)if(v[p]){for(int q = 0; q < nB; q++)B[q] = p;return nV;}};
    if(nV <= nB){
//        PRINTF("Warning: Number of violators(%i) is less than dimension of working set(%i).\n",nV,nB);
        int q = 0;for(int p = 0; p < n1n2; p++)if(v[p])B[q++] = p;
        return nV;        
    }          
    if(nV == 2){
        int i = 0;
        for(int p = 0; p < n1n2; p++) if(v[p]) B[i++] = p;   
        return 2;   
    }    
    int nC = 0,n0 = 0;
    if(nVC == 0){
        n0 = nB;
    }else if(nV0 == 0){
        nC = nB;
    }else{
        int half = nB / 2;
        nC = MIN(half,nVC);
        n0 = MIN(half,nV0);
        if(nC < half) n0 = nB - nC;else if(n0 < half) nC = nB - n0;
   }       
//        PRINTF("n0 %i nC %i nV0 %i nVC %i nB %i\n",n0,nC,nV0,nVC,nB);

    if(n0){
        int* B_0 = B;
        double* M_0 = x;//
		//double M_0[n0];
		for(int p = 0; p < n0; p++)M_0[p] = _Inf;
        int size_0 = 0;
        int which_max_0 = 0;   
        for(int p = 0; p < n1n2; p++) {       
            if(v[p] == 1) {
                if(size_0 == n0){
                    if(-df[p] < M_0[which_max_0]){
                        M_0[which_max_0] = -df[p];
                        B_0[which_max_0] = p;
                        which_max_0 = 0;
                        for(int k = 1;k < n0;k++)
                            if(M_0[which_max_0] > M_0[k])
                                which_max_0 = k;
                    }
                }else{
                    B_0[size_0] = p;
                    M_0[size_0++] = -df[p];
                    if(size_0 == n0){
                        which_max_0 = 0;
                        for(int k = 1;k < n0;k++)
                            if(M_0[which_max_0] > M_0[k]) 
                                which_max_0 = k;
                    }
                }
            }
        }    
    }
     
    if(nC){
        int* B_C = B + n0;
        double* M_C = x;//
		//double M_C[nC];
		for(int p = 0; p < nC; p++)M_C[p] = -_Inf; 
        int size_C = 0;
        int which_min_C = 0;   
        for(int p = 0; p < n1n2; p++) {       
            if(v[p] == 2) {
                if(size_C == nC){
                    if(-df[p] > M_C[which_min_C]){
                        M_C[which_min_C] = -df[p];
                        B_C[which_min_C] = p;
                        which_min_C = 0;
                        for(int k = 1;k < nC;k++)
                            if(M_C[which_min_C] < M_C[k]) 
                                which_min_C = k;
                    }
                }else{
                    B_C[size_C] = p;                
                    M_C[size_C++] = -df[p];
                    if(size_C == nC){
                        which_min_C = 0;
                        for(int k = 1;k < nC;k++)
                            if(M_C[which_min_C] < M_C[k]) 
                                which_min_C = k;
                    }
                }
            }
        }    
    }
    return n0 + nC;
   
/*
// Doing a full sort is better if size of working set ~ n1n2 but not of 'fixed' sized working sets   
    int* V = which_v;
    int p;
//    double x[MAX(nVC,nV0)];int perm[MAX(nVC,nV0)];int V[MAX(nVC,nV0)];
    if(nC){
        int q = 0;
        for(p = 0; p < n1n2; p++) if(v[p] == 2)V[q++] = p;
        for(p = 0; p < nVC; p++)x[p] = -df[V[p]];
        for(p = 0; p < nVC; p++)perm[p] = p+1;  
        revsort(&x[0], &perm[0], nVC); 
        for(p = 0; p < nC; p++)B[p] = V[perm[p] - 1];
    }
    if(n0){
        int q = 0;
        for(p = 0; p < n1n2; p++) if(v[p] == 1)V[q++] = p;
        for(p = 0; p < nV0; p++)x[p] = df[V[p]];
        for(p = 0; p < nV0; p++)perm[p] = p+1;  
        revsort(&x[0], &perm[0], nV0);
        for(p = nC; p < n0+nC; p++)B[p] = V[perm[p-nC]-1];
    }
    return n0 + nC;
*/	
}


// // weights based on rank(fabs(df))    
// // v[p] > 0 if {p: 0 < a[p],0 < df[p]} U {p: C > a[p],0 > df[p]}
// int getWorkingSet_rvwr(int nB,int* B,int nV0,int nVC,int n1n2,double* df,int* v,int* which_v,int* ibuf,int* ibuf2,double* dbuf) {

    // int p;    
    // double s = 0.0; 
    // int nV = nV0 + nVC;
// //  PRINTF("Random working set selection n(%i) nV(%i) nVC(%i) nV0(%i): ",n1n2,nV,nVC,nV0);
    // if(nV == 0) return 0;
    // if(nV == 1){for(int p = 0; p < n1n2; p++)if(v[p]){for(int q = 0; q < nB; q++)B[q] = p;return nV;}};
    // if(nV <= nB){
// //        PRINTF("Warning: Number of violators(%i) is less than dimension of working set(%i).\n",nV,nB);
        // int q = 0;for(int p = 0; p < n1n2; p++)if(v[p])B[q++] = p;
        // return nV;        
    // }          
    // if(nV == 2){
        // int i = 0;
        // for(int p = 0; p < n1n2; p++) if(v[p]) B[i++] = p;   
        // return 2;   
    // }    
    // int nC = 0,n0 = 0;
    // if(nVC == 0){
        // n0 = nB;
    // }else if(nV0 == 0){
        // nC = nB;
    // }else{
        // int half = nB / 2;
        // nC = MIN(half,nVC);
        // n0 = MIN(half,nV0);
        // if(nC < half) n0 = nB - nC;else if(n0 < half) nC = nB - n0;
   // }       
// //        PRINTF("n0 %i nC %i nV0 %i nVC %i nB %i\n",n0,nC,nV0,nVC,nB);

// //    int perm[nV];double pV[nV];int r[nV];int V[nV];
    // int* perm = ibuf;int* r = ibuf2;int* V = which_v;double* pV = dbuf;
    
    // int q = 0;for(p = 0; p < n1n2; p++) if(v[p]) V[q++] = p;

    // for(p = 0; p < nV; p++)pV[p] = fabs(df[V[p]]);
    // rank((double*)&pV[0],(int*)&r[0], &nV);
    // s = 0.5 * (double)nV * ((double)nV + 1.0);
    // for(p = 0; p < nV; p++)pV[p] = (double)r[p] / s;  
    // ProbSampleNoReplace(nV,&pV[0],&(perm[0]),nB,B);   

    // for(p = 0; p < nB; p++)B[p] = V[B[p] - 1];    
    // return nB;
   
// }

/*
# R code to illustrate the indexing used in working set selection methods
# where we sample/select a subset of df[] of length 'n', sort it then find 
# corresponding elements in {1..n}
df <- c(10,3,8,1,2,5,4,7,6,9)
v <- c(0,1,1,0,1,0,1,0,0,1)
V <- numeric(sum(v));
j=1;
for(i in 1:length(v)){if(v[i]){V[j] = i;j=j+1}}
x <- df[as.logical(v)]
y <- .C("R_heap_sort",as.double(x),as.integer(1:length(x)), as.integer(length(x)))
df[V[y[[2]]]]
*/
