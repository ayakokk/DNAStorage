// vector/array
void PrintVect( const int           *V, int len, const char *pre, const char *post);
void PrintVect( const unsigned int  *V, int len, const char *pre, const char *post);
void PrintVect2(const unsigned char *V, int len, const char *pre, const char *post);
void PrintVect2(const bool          *V, int len, const char *pre, const char *post);
void PrintSP(int len);
void SetVal(bool *V, bool val, int len);
void ClearArray(int **V, int lr, int lc);
//---
int MaxRunLength( const unsigned char *V, int len); 
int HammingWeight(const unsigned char *V, int len);
int HammingWeight(const bool          *V, int len); // count true
int HammingDist(const unsigned char *V0, const unsigned char *V1, int len);
int HammingDist(const bool          *V0, const bool          *V1, int len);

// conversion
void ConvIntV4(unsigned char *V, unsigned int u, int len); // unsigned int -> 4-ary vector
unsigned int ConvV4Int(const unsigned char *V, int len);   // 4-ary vector -> unsigned int
void ConvIntV2(unsigned char *V, unsigned int x, int len); // unsigned int -> binary vector
unsigned int ConvV2Int(const unsigned char *V, int len);   // binary vector-> unsigned int

// max,min,mod
int max(int v0, int v1);
int min(int v0, int v1);
int argmax(const int *V, int len);
int argmin(const int *V, int len);
int max(const int *V, int len);
int min(const int *V, int len);
int mod(int a, int m);   // a mod m
