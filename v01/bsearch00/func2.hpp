void PrintVect(const unsigned char *V, int len, const char *pre, const char *post);
void PrintVect(const int           *V, int len, const char *pre, const char *post);
int  max(const int *V, int len);
int  binom(int n, int k);
void ConvUintVect(unsigned char *V, unsigned int u, unsigned int len);  // V <- u
unsigned int ConvVectUint(const unsigned char *V,   unsigned int len);  // u <- V
