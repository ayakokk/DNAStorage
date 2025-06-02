void ClearVect(unsigned long *V, int len);
void PrintVect(const unsigned char *V, int len);
void PrintVect4(const unsigned char *V, int len);
void ConvIntVect(int u, unsigned char *V, int len); // unsigneg int -> 4-ary vector
int  ConvVectInt(const unsigned char *V, int len);  // 4-ary vector -> unsigned int
void ConvToBinVect(unsigned char *V, unsigned int x, int len); // V = binary expression of x
unsigned int ConvFromBinVect(const unsigned char *V, int len); //    (inverse)
int  SymbolCnt(const unsigned char *V, unsigned char val, int len);
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
int  HammingDist(const int *V0, const int *V1, int len);
int  HammingDist(const unsigned char *V0, const unsigned char *V1, int len);
void Normalize(double *V, int len);
