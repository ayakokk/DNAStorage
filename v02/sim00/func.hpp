int max(int a, int b);
int min(int a, int b);
int argmax(const double *V, int len);
int argmin(const double *V, int len);
double max(const double *V, int len);
double min(const double *V, int len);

void PrintVect(const double *V, int len, const char *pre, const char *post);
void PrintVectX(const int   *V, int len, const char *pre, const char *post);
void PrintVectB(const unsigned char *V, int len, const char *pre, const char *post);
void RandVect(int *V, int len, int Vmin, int Vmax);

long VectToLong(const unsigned char *V, int len);
void LongToVect(unsigned char *V, long val, int len);

void ReadConstraints(const char *fn, int *Rho, int *ell, int *Delta);
void HardDecision(int *V, const double **P, int N, int Q);
int  HammingDist(const int *V0, const int *V1, int len);
