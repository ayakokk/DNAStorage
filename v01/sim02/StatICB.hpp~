class StatICB {
private:
  int    nu,beta,beta4p,Bsum;
  int    HdThNum;
  double HdThStep;
  double *HdThList;
  int    *B;                 // info size (bits)
  unsigned long ****HdCnt;   // counter [HdThNum][nu][X:beta4p][Y:beta4p+1]
  //-----
  void   InitHdCnt();
  double CalcIxy(const unsigned long **C, int numX, int numY);
  void   SetY(int *Y, const double **P, double Pth, int Nseg);
  double Entropy(const double *P, int len);
  int    ArgMax(const double *V, int len);
  void   ClearVect(unsigned long *V, int len);
public:
  StatICB(int _nu, int _beta, const int *_B);
  ~StatICB();
  void count(const int *X, const double **P, int Nseg);
  void PrintIxy();
  void dump();
};
