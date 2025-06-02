class StatICB {
private:
  int    nu,beta,beta4p,Bsum;
  int    lmd;                 // number of quantization levels
  int    *B;                  // info size (bits)
  unsigned long ***HdCnt;    // counter [nu][X:beta4p][Y:beta4p*lmd]
  //-----
  void   InitHdCnt();
  double CalcIxy(const unsigned long **C, int numX, int numY);
  void   SetY( int *Y, const double **P, int Nseg);
  void   SetHp(int *Hp,const double **P, int Nseg);
  double Entropy(const double *P, int len);
  int    ArgMax(const double *V, int len);
  void   ClearVect(unsigned long *V, int len);
public:
  StatICB(int _nu, int _beta, const int *_B, int _lmd);
  ~StatICB();
  void count(const int *X, const double **P, int Nseg);
  void PrintIxy();
  void dump();
};
