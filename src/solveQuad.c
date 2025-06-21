#include <string.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <R.h>

#define PRINTF printf
#ifdef Rprintf
#undef PRINTF
#define PRINTF Rprintf
#endif


//#define __DEBUG__SOLVEQUAD   
#ifdef __DEBUG__SOLVEQUAD    
  // for columns to align, all data arguments are floating type, column names are strings
  // MAX_STRING = MAX_LEFT and MAX_RIGHT is s.t. MAX_LEFT + MAX_RIGHT >= the length of the longest possible formatted number.
  // Thus if %f.2 then -1234.123 has length 8, if %e #digits(LHS) + #digits(RHS) of a %f number 
  #define MAX_STRING 10
  #define MAX_LEFT 10
  #define MAX_RIGHT 1
  #include <stdio.h>
  #include <R.h>
  double ushape(double a,double b,double c,double x){return a*x*x + b*x + c;}
#endif

// Q is (n x n) matrix from R, represented as n*n length array
// Q = as.double(t(Q)) from R is in row-order so that Qij = Q[i,j] = Q[i*n + j]
// a is an array of length n
// B = {\i}, dim(Q_{iB}) = (1 x (n-1)), dim(a_B) = ((n-1) x 1)  
// NOTE: 'malloced' arrays, "memcpy(a0,a,n*sizeof(a));"   
//      'VLA' arrays        "memcpy(&(a[0]),&(a[0]),n * sizeof(double));"

  
// working set size = 1 	
void solve_Quad(
    double* _C,
	int* _n1n2,
    double* b, //n1n2
    int* bigQ,
    double* _Q,//n1n2 x n1n2    
    int* _n1,
    int* _n2,
    double* _K,//(n1+n2) x (n1+n2)
	double* _a,//n1n2 input alpha 
    int* maxit,double* _tol,int* trace,  //input
	double* a,//n1n2 output 
    double* v, int* iter,int* convergence,double * epsilon  //output
	)
{        
    int p,q,done;
    double C = *_C;
    int n1n2 = *_n1n2;
    int n1 = *_n1;
    int n2 = *_n2;
    int n = n1+n2;
    double c[n1n2];  
    double a0;
    double da;
    double v0;
    #ifdef __DEBUG__SOLVEQUAD    
    double c0;  
    double raw_a;
    #endif
    double reltol = *_tol;
    double** Q = NULL;//this matrix can be rather big hence allocate it on the heap    
    if(*bigQ) {
     if(n1n2 <= 0){Rprintf("Allocate size n1n2 = 0 memory?\n");return;}
      Q = (double**) calloc(n1n2,sizeof(double*));
      if(!Q){Rprintf("solveQuad: Unable to allocate %i(bytes)\n",n1n2*sizeof(double*));return;}      
      for(p = 0; p < n1n2; p++)Q[p] = _Q + p*n1n2;
    }
    double (*K)[n]   = (double (*)[n])_K;
    double Q_pp = 0.0;
    double Q_pq = 0.0;
    
    // initialization
    *v = 0.0;*convergence = 1;*iter = 0;
    if(reltol < 0.0) reltol = 0.0;
    *epsilon = 2.0 * (reltol);
    memcpy(a,_a,n1n2 * sizeof(double));

    int ip,jp,iq,jq;
    for(p = 0;p < n1n2;p++){
        ip=(int)ceil((p+1.0)/n2); jp=(p+1)-(ip-1)*n2; 
        if(!Q){
            Q_pp = K[ip-1][ip-1]+K[jp+n1-1][jp+n1-1]-K[ip-1][jp+n1-1]-K[jp+n1-1][ip-1];
        }else Q_pp = Q[p][p];
        c[p] = -Q_pp * a[p];//c_i = < Q_iB, a_B >   
        for(q = 0; q < n1n2; q++) {
            if(!Q){
                iq=(int)ceil((q+1.0)/n2); jq=(q+1)-(iq-1)*n2;                     
                Q_pq = K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1];
            } else Q_pq = Q[p][q];
            c[p] += Q_pq * a[p];
             *v  += Q_pq * a[p] * a[q];//a'Qa               
        }
        *v -= b[p] * a[p];//a'Qa - b'a
    }    
    

    *iter = 0;   
    done = 0;
    while(!done){
        (*iter)++; 
        v0 = *v;   

        #ifdef __DEBUG__SOLVEQUAD        
        char gC[MAX_STRING];sprintf(gC,"g(%+.1f)",C);
        PRINTF("%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n",
               MAX_STRING,"v", MAX_STRING,"dv",MAX_STRING,"a0", MAX_STRING,"a",MAX_STRING,"da",
               MAX_STRING,"g(a0)", MAX_STRING,"g(a)",MAX_STRING,"g(0)",MAX_STRING,gC,
               MAX_STRING,"coef(x^2)",MAX_STRING,"coef(x)",MAX_STRING,"b",MAX_STRING,"c",MAX_STRING,"c0",
               MAX_STRING,"ra",MAX_STRING,"g(ra)");
        #endif  

        
        for(p = 0;p < n1n2;p++){      
            int ip=ceil((p+1.0)/n2); int jp=(p+1)-(ip-1)*n2; 
            if(!Q){
                Q_pp = K[ip-1][ip-1]+K[jp+n1-1][jp+n1-1]-K[ip-1][jp+n1-1]-K[jp+n1-1][ip-1];
            }else Q_pp = Q[p][p];


            ///////////////////////////////////////////////////////
            // new a[p] = median({0, C, (0.5b[p] - c[p]) / Qii} )
              a0 = a[p];
              a[p] = (0.5*b[p] - c[p]) / Q_pp;
              #ifdef __DEBUG__SOLVEQUAD    
                c0 = c[p];
                raw_a = a[p]; 
              #endif
              if(a[p] > C) a[p] = C;else if(a[p] < 0.0) a[p] = 0.0;       
              da = a[p] - a0;      
            ///////////////////////////////////////////////////////

            ///////////////////////////////////////////////////////
            // update c and value      
              if(da != 0.0){         
                for(q = 0;q < n1n2;q++){
                    if(!Q){
                        iq=ceil((q+1.0)/n2); jq=(q+1)-(iq-1)*n2;                     
                        Q_pq = K[ip-1][iq-1]+K[jp+n1-1][jq+n1-1]-K[ip-1][jq+n1-1]-K[jp+n1-1][iq-1];
                    } else Q_pq = Q[p][q];
                    c[q] +=  Q_pq * da;
                }
                c[p] -= Q_pp * da;
                *v += (Q_pp * (a[p] + a0) + 2.0*c[p] - b[p]) * da; 
              }
            ///////////////////////////////////////////////////////
              
            #ifdef __DEBUG__SOLVEQUAD        
              PRINTF("%+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e\n",
                     MAX_LEFT,MAX_RIGHT,*v, 
                     MAX_LEFT,MAX_RIGHT,dv, 
                     MAX_LEFT,MAX_RIGHT,a0,
                     MAX_LEFT,MAX_RIGHT,a[p],
                     MAX_LEFT,MAX_RIGHT,da,
                     MAX_LEFT,MAX_RIGHT,ushape(Q_pp,2.0*c[p] - b[p],0.0,a0),
                     MAX_LEFT,MAX_RIGHT,ushape(Q_pp,2.0*c[p] - b[p],0.0,a[p]),
                     MAX_LEFT,MAX_RIGHT,ushape(Q_pp,2.0*c[p] - b[p],0.0,0.0),             
                     MAX_LEFT,MAX_RIGHT,ushape(Q_pp,2.0*c[p] - b[p],0.0,C),
                     MAX_LEFT,MAX_RIGHT,Q_pp,
                     MAX_LEFT,MAX_RIGHT,2.0*c[p] - b[p],
                     MAX_LEFT,MAX_RIGHT,b[p],
                     MAX_LEFT,MAX_RIGHT,c[p],
                     MAX_LEFT,MAX_RIGHT,c0,                        
                     MAX_LEFT,MAX_RIGHT,raw_a,                        
                     MAX_LEFT,MAX_RIGHT,ushape(Q_pp,2.0*c[p] - b[p],0.0,raw_a)
                    );              
            #endif
        }//end for{p}

        //    if(v0 == 0) *epsilon = fabs(*v);else *epsilon = fabs((*v - v0)/v0);
          *epsilon = fabs(*v - v0);

        #ifdef __DEBUG__SOLVEQUAD    
        PRINTF("%-13s %-13s %-13s %-13s\n","iter","eps","v","v0");    
        PRINTF("%+10E %+10E %+10E %+10E\n\n",(double)*iter,*epsilon,*v,v0);
        #endif    
        if(*trace > 1) Rprintf("iter(%i) v(%+.20e) e(%+.20e) \n",*iter,*v,*epsilon);

        if(*epsilon <= *_tol){
          done = 1;*convergence = 0;
          break;		
        }else if(*iter == *maxit){ 
          done = 1;*convergence = 1;break;
        }
    } // end while loop  
       
    if(*trace) Rprintf("iter(%i) v(%+.3e) e(%+.3e) \n",*iter,*v,*epsilon);	  

    if(Q){free(Q);Q = NULL;}
}


// indexing Q[i,j] = Q[i*n + j]  
// working set size = 1 	
void solve_Quad1(
	int* _n,double* Q,double* b,double* _C,int* maxit,double* _reltol,double *alpha0,int* trace,  //input
	double* a, double* v, int* iter,int* convergence,double * epsilon  //output
	)
{        
  int i,j,done;
  int n = *_n;
  double c[n];  
  double a0;
  double delta_i;
  double v0;
  double C = *_C;
#ifdef __DEBUG__SOLVEQUAD    
  double c0;  
  double raw_a;
#endif
  double reltol = *_reltol;
  
  // initialization
  *v = 0.0;*convergence = 1;*iter = 0;
  if(reltol < 0.0) reltol = 0.0;
  *epsilon = 2.0 * (reltol);
  memcpy(a,alpha0,n * sizeof(double));

  for(i = 0;i < n;i++){
    double* Q_i = Q + i * n;
    c[i] = -Q_i[i] * a[i];
    for(j = 0;j < n;j++){
      c[i] += Q_i[j] * a[j]; //c = < Q_iB, a_B >   
      *v   += Q_i[j] * a[i] * a[j];//a'Qa
    }
    *v -= b[i] * a[i];//a'Qa - b'a
  }

  *iter = 0;   
  done = 0;
  while(!done){
    (*iter)++; 
    v0 = *v;   

#ifdef __DEBUG__SOLVEQUAD        
    char gC[MAX_STRING];sprintf(gC,"g(%+.1f)",C);
    PRINTF("%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n",
           MAX_STRING,"v", MAX_STRING,"dv",MAX_STRING,"a0", MAX_STRING,"a",MAX_STRING,"da",
           MAX_STRING,"g(a0)", MAX_STRING,"g(a)",MAX_STRING,"g(0)",MAX_STRING,gC,
           MAX_STRING,"coef(x^2)",MAX_STRING,"coef(x)",MAX_STRING,"b",MAX_STRING,"c",MAX_STRING,"c0",
           MAX_STRING,"ra",MAX_STRING,"g(ra)");
#endif  
    for(i = 0;i < n;i++){      
      double Qii = *(Q + i * n + i);  
///////////////////////////////////////////////////////
// new a[i] = median({0, C, (0.5b[i] - c[i]) / Qii} )
      a0 = a[i];
      a[i] = (0.5*b[i] - c[i]) / Qii;
      #ifdef __DEBUG__SOLVEQUAD    
        c0 = c[i];
        raw_a = a[i]; 
      #endif
      if(a[i] > C) a[i] = C;else if(a[i] < 0.0) a[i] = 0.0;       
      delta_i = a[i] - a0;      
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
// update c and value      
      double dv = 0.0; 
      if(delta_i != 0.0){         
        for(j = 0;j < n;j++){
          double Qji = *(Q + j*n + i);
          c[j] +=  Qji * delta_i;
        }
        c[i] -= Qii * delta_i;
        dv = (Qii * (a[i] + a0) + 2.0*c[i] - b[i]) * delta_i; 
        *v += dv;       
      }
///////////////////////////////////////////////////////
      
#ifdef __DEBUG__SOLVEQUAD        
      PRINTF("%+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e\n",
             MAX_LEFT,MAX_RIGHT,*v, 
             MAX_LEFT,MAX_RIGHT,dv, 
             MAX_LEFT,MAX_RIGHT,a0,
             MAX_LEFT,MAX_RIGHT,a[i],
             MAX_LEFT,MAX_RIGHT,delta_i,
             MAX_LEFT,MAX_RIGHT,ushape(Qii,2.0*c[i] - b[i],0.0,a0),
             MAX_LEFT,MAX_RIGHT,ushape(Qii,2.0*c[i] - b[i],0.0,a[i]),
             MAX_LEFT,MAX_RIGHT,ushape(Qii,2.0*c[i] - b[i],0.0,0.0),             
             MAX_LEFT,MAX_RIGHT,ushape(Qii,2.0*c[i] - b[i],0.0,C),
             MAX_LEFT,MAX_RIGHT,Qii,
             MAX_LEFT,MAX_RIGHT,2.0*c[i] - b[i],
             MAX_LEFT,MAX_RIGHT,b[i],
             MAX_LEFT,MAX_RIGHT,c[i],
             MAX_LEFT,MAX_RIGHT,c0,                        
             MAX_LEFT,MAX_RIGHT,raw_a,                        
             MAX_LEFT,MAX_RIGHT,ushape(Qii,2.0*c[i] - b[i],0.0,raw_a)
            );              
#endif
    }//end for{i}
    
//    if(v0 == 0) *epsilon = fabs(*v);else *epsilon = fabs((*v - v0)/v0);
      *epsilon = fabs(*v - v0);
	
#ifdef __DEBUG__SOLVEQUAD    
    PRINTF("%-13s %-13s %-13s %-13s\n","iter","eps","v","v0");    
    PRINTF("%+10E %+10E %+10E %+10E\n\n",(double)*iter,*epsilon,*v,v0);
#endif    
	if(*trace > 1) Rprintf("i(%i) obj(%+.20e) eps(%+.20e) \n",*iter,*v,*epsilon);

    if(*iter == *maxit){ 
      done = 1;*convergence = 1;
    }else if((1 < *iter) & (*epsilon <= reltol)){
      done = 1;*convergence = 0;
    }else{;}
  } // end while loop  
//	if(*trace == 1) Rprintf("i(%i) obj(%+.3e) eps(%+.3e) \n",*iter,*v,*epsilon);	  
}

// working set size = 1 
// ever-slightly faster due to pointer arithmetic p++ vs p[i]
void solve_Quad_2(
        int* _n,double* _Q,double* _b,double* _C,int* maxit,double* _reltol,double* alpha0,int* trace,int* update,//input
        double* _a, double* v, int* iter,int* convergence,double * epsilon  //output
        )
{        
  int i,j,done;
  int n = *_n;
  double _c[n];
  double *ci,*bi,*ai;
  double *Qii,*Qij;
  double a0;
  double delta_i;
  double v0;
  double C = *_C;
#ifdef __DEBUG__SOLVEQUAD    
  double c0;  
  double raw_a;
#endif
  double reltol = *_reltol;    
  // initialization
  *v = 0.0;*convergence = 1;*iter = 0;
  if(reltol < 0.0) reltol = 0.0;
  *epsilon = 2.0 * (reltol);
  memcpy(_a,alpha0,n * sizeof(double));
  Qii = Qij = _Q;
  ai = _a;
  bi = _b;
  ci = &(_c[0]);//&[0] is needed for VLA-type, as opposed to [0] 
  *v = 0.0;
  for(i = 0;i < n;i++){
    *ci = -(*(Qii)) * (*ai); //Q = Q_i
    double* aj = _a;
    for(j = 0;j < n;j++){
      *ci += (*Qij) * (*aj); //c = < Q_iB, a_B >   
      *v  += (*Qij) * (*ai) * (*aj);//a'Qa
      Qij++;aj++;
    };
    *v -= (*bi) * (*ai);//a'Qa - b'a
    Qii += n+1;ai++;bi++;ci++;
   };
  
  *iter = 0;   
  done = 0;
  while(!done){
    (*iter)++; 
    v0 = *v;   

#ifdef __DEBUG__SOLVEQUAD        
    char gC[MAX_STRING];sprintf(gC,"g(%+.1f)",C);
    PRINTF("%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n",
           MAX_STRING,"v", MAX_STRING,"dv",MAX_STRING,"a0", MAX_STRING,"a",MAX_STRING,"da",
           MAX_STRING,"g(a0)", MAX_STRING,"g(a)",MAX_STRING,"g(0)",MAX_STRING,gC,
           MAX_STRING,"coef(x^2)",MAX_STRING,"coef(x)",MAX_STRING,"b",MAX_STRING,"c",MAX_STRING,"c0",
           MAX_STRING,"ra",MAX_STRING,"g(ra)");
#endif

    Qii  = _Q;
    ai = _a;
    bi = _b;
    ci = &_c[0];

    for(i = 0;i < n;i++){              
///////////////////////////////////////////////////////
// new a[i] = median({0, C, (0.5b[i] - c[i]) / Qii} )
      a0 = *ai;
      *ai = (0.5*(*bi) - (*ci)) / *Qii;
      #ifdef __DEBUG__SOLVEQUAD    
        c0 = *ci;
        raw_a = *ai; 
      #endif
      if(*ai > C) *ai = C;else if(*ai < 0.0) *ai = 0.0;       
      delta_i = *ai - a0;      
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
// update c and value      
      double dv = 0.0; 
      if(delta_i != 0.0){         
        double* Qji = _Q + i;
        double* cj = &(_c[0]); 
        for(j = 0;j < n;j++){
          *cj +=  (*Qji) * delta_i;
          Qji +=n;cj++;
        }
        *ci -= *(Qii) * delta_i;
        dv = (*(Qii) * (*ai + a0) + 2.0*(*ci) - (*bi)) * delta_i; 
        *v += dv;       
      }
///////////////////////////////////////////////////////
     
#ifdef __DEBUG__SOLVEQUAD        
      PRINTF("%+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e %+*.*e\n",
             MAX_LEFT,MAX_RIGHT,*v, 
             MAX_LEFT,MAX_RIGHT,dv, 
             MAX_LEFT,MAX_RIGHT,a0,
             MAX_LEFT,MAX_RIGHT,*ai,
             MAX_LEFT,MAX_RIGHT,delta_i,
             MAX_LEFT,MAX_RIGHT,ushape(*(Qii),2.0*(*ci) - (*bi),0.0,a0),
             MAX_LEFT,MAX_RIGHT,ushape(*(Qii),2.0*(*ci) - (*bi),0.0,*ai),
             MAX_LEFT,MAX_RIGHT,ushape(*(Qii),2.0*(*ci) - (*bi),0.0,0.0),             
             MAX_LEFT,MAX_RIGHT,ushape(*(Qii),2.0*(*ci) - (*bi),0.0,C),
             MAX_LEFT,MAX_RIGHT,*(Qii),
             MAX_LEFT,MAX_RIGHT,2.0*(*ci) - (*bi),
             MAX_LEFT,MAX_RIGHT,(*bi),
             MAX_LEFT,MAX_RIGHT,*ci,
             MAX_LEFT,MAX_RIGHT,c0,                        
             MAX_LEFT,MAX_RIGHT,raw_a,                        
             MAX_LEFT,MAX_RIGHT,ushape(*(Qii),2.0*(*ci) - *bi,0.0,raw_a)
            );              
#endif
    Qii += n + 1;ai++;bi++;ci++;
    }//end for{i}
    
    if(v0 == 0.0) *epsilon = fabs(*v);else *epsilon = fabs((*v - v0)/v0);
#ifdef __DEBUG__SOLVEQUAD    
    PRINTF("%-13s %-13s %-13s %-13s\n","iter","eps","v","v0");    
    PRINTF("%+10E %+10E %+10E %+10E\n\n",(double)*iter,*epsilon,*v,v0);
#endif    

    if(*iter == *maxit){ 
      done = 1;*convergence = 1;
    }else if((1 < *iter) & (*epsilon <= reltol)){
      done = 1;*convergence = 0;  
    }else{;}
  } // end while loop  
}
