// functions lifted from R's source code, except for rank()
void rank(double *a, int *rank, int* _n);
void ProbSampleNoReplace(int n, double *p, int *perm,int nans, int *ans);
void ProbSampleReplace(int n, double *p, int *perm, int nans, int *ans);
void SampleNoReplace(int k, int n, int *y, int *x);
void SampleReplace(int k, int n, int *y);

