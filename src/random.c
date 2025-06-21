#include "sort.h"
#include <R.h>
#ifndef SCYTHE_COMPILE_DIRECT
#include <Rmath.h>
#define RUNIF runif
#else
#include "mersenne.h" 
#include "rng.h" 
scythe::mersenne _rng; 
double RUNIF(double,double) { return _rng.runif(); }
#endif

// taken from R's src/main/random.c
// 'perm' \in {1,..,n} is pre-allocated of length n
void ProbSampleNoReplace(int n, double *p, int *perm,int nans, int *ans)
{
    double rT, mass, totalmass;
    int i, j, k, n1;
    
#ifndef SCYTHE_COMPILE_DIRECT    
    GetRNGstate();    
#endif
    /* Record element identities */
    for (i = 0; i < n; i++)
	perm[i] = i + 1;

    /* Sort probabilities into descending order */
    /* Order element identities in parallel */
    revsort(p, perm, n);

    /* Compute the sample */
    totalmass = 1.0;
    for (i = 0, n1 = n-1; i < nans; i++, n1--) {
        rT = totalmass * RUNIF(0.0,1.0);
        mass = 0.0;
        for (j = 0; j < n1; j++) {
            mass += p[j];
            if (rT <= mass)
            break;
        }
        ans[i] = perm[j];
        totalmass -= p[j];
        for(k = j; k < n1; k++) {
            p[k] = p[k + 1];
            perm[k] = perm[k + 1];
        }
    }
#ifndef SCYTHE_COMPILE_DIRECT    
    PutRNGstate();    
#endif	
}

// taken from R's src/main/random.c
// 'x' is preallocated of length n
// output: 'y' of length 'k' preallocated
void SampleNoReplace(int k, int n, int *y, int *x)
{
    int i, j;
#ifndef SCYTHE_COMPILE_DIRECT    
    GetRNGstate();    
#endif
    for (i = 0; i < n; i++) x[i] = i;
    for (i = 0; i < k; i++) {
    	j = (int)((double)n * RUNIF(0.0,1.0));
 	    y[i] = x[j] + 1;
 	    x[j] = x[--n];
    }
#ifndef SCYTHE_COMPILE_DIRECT    
    PutRNGstate();    
#endif
}

// taken from R's src/main/random.c
/* Unequal probability sampling; with-balancedment case */
void ProbSampleReplace(int n, double *p, int *perm, int nans, int *ans)
{
    double rU;
    int i, j;
    int nm1 = n - 1;

#ifndef SCYTHE_COMPILE_DIRECT    
    GetRNGstate();    
#endif

    /* record element identities */
    for (i = 0; i < n; i++)
	perm[i] = i + 1;

    /* sort the probabilities into descending order */
    revsort(p, perm, n);

    /* compute cumulative probabilities */
    for (i = 1 ; i < n; i++)
	p[i] += p[i - 1];

    /* compute the sample */
    for (i = 0; i < nans; i++) {
	rU = RUNIF(0.0,1.0);
	for (j = 0; j < nm1; j++) {
	    if (rU <= p[j])
		break;
	}
	ans[i] = perm[j];
    }
#ifndef SCYTHE_COMPILE_DIRECT    
    PutRNGstate();    
#endif    
}
// taken from R's src/main/random.c
/* Equal probability sampling; with-balancedment case */
void SampleReplace(int k, int n, int *y)
{
#ifndef SCYTHE_COMPILE_DIRECT    
    GetRNGstate();    
#endif
    int i;
    for (i = 0; i < k; i++) y[i] = (int)((double)n * RUNIF(0.0,1.0)) + 1;
#ifndef SCYTHE_COMPILE_DIRECT    
    PutRNGstate();    
#endif   
}

// ties are not handled
void rank(double *a, int *rank, int* _n){
    int n = *_n;
    double _a[n];
    for(int i = 0;i < n;i++) rank[i] = i+1;
    for(int i = 0;i < n;i++)_a[i] = -a[i];
    revsort((double*)&(_a[0]),rank,n);////since revsort() sorts in descending order only
    for(int i = 0;i < n;i++)_a[i] = -(double)rank[i];
    revsort((double*)&(_a[0]),rank,n); 
}  

