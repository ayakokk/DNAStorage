void ClearVect(unsigned long *V, int len);
void ClearArray(int **V, int lr, int lc);
void PrintVect(const unsigned char *V, int len);
void PrintVect(const unsigned long *V, int len, const char *pre, const char *post);
void PrintVect(const int           *V, int len, const char *pre, const char *post);
void PrintVect(const unsigned int  *V, int len, const char *pre, const char *post);
void PrintVect(const double        *V, int len, const char *pre, const char *post);
void PrintVect4(const unsigned char *V, int len);
void PrintVectRatio(const unsigned long *V, int len);
void PrintArray(const double **V, int l0, int l1);
void PrintArray(const int    **V, int l0, int l1);
void GenRndVect(int *V, int mod, int len);
void ConvIntVect(int u, unsigned char *V, int len); // unsigneg int -> 4-ary vector
int  ConvVectInt(const unsigned char *V, int len);  // 4-ary vector -> unsigned int
void ConvToBinVect(unsigned char *V, unsigned int x, int len); // V = binary expression of x
unsigned int ConvFromBinVect(const unsigned char *V, int len); //    (inverse)
int  SymbolCnt(const unsigned char *V, unsigned char val, int len);
int    GCweight(    const unsigned char *V, int len); 
int    GCbalance(   const unsigned char *V, int len);          // GC-balance wG + wC - wA- wT 
double GCbalanceWin(const unsigned char *V, int len, int ws);  // GC-balance wG + wC - wA- wT (windowed,abs,ave)
void   GCbalanceCnt(const unsigned char *V, int len, int ws, unsigned long *BLC);
int    MaxRunLength(const unsigned char *V, int len); // maximum run-length
double AveRunLength(const unsigned char *V, int len); // average run-length
void   CntRunLength(const unsigned char *V, int len, unsigned long *RLC, int RLmax); // set counter
void RandomSel(class bmatrix *CW, int num, class bmatrix *F);  // select num elements from F=0 
void Swap(int *a, int *b);
int  min(int a, int b);
int  max(int a, int b);
int  min(const int *V, int len);
int  max(const int *V, int len);
int  ArgMin(const double *V, int len);
int  ArgMax(const double *V, int len);
int    Sum(const int    *V, int len);
double Sum(const double *V, int len);
int    HammingDist(const int *V0, const int *V1, int len);
int    HammingDist(const unsigned char *V0, const unsigned char *V1, int len);
double Normalize(double *V, int len);  // returns sum
double Entropy(const int *X, int len);
double ErrEntropy(const int *X, const double **P, int len, int Q);
