class FBalg4 {
private:
  int    Nseg,Dmax;
  int    beta,beta4p,N;
  int    DmaxW;     // Dmax*2+1
  double Pid,Ps,Pth;
  class Gtable *GTB;
  double **Pxin, **Pxout;  // Px:    [Nseg][beta4p]
  double **PF, **PB;       // Drift: [Nseg+1][DmaxW]
  unsigned char *Y;        // received word [N+Dmax]
  unsigned char *dbgX;     // (dbg)X: [N]
  //---
  void ArrayCopy(double **Qt, const double **Qs, int len0, int len1);
  void init(const double **Pin, const unsigned char *Yin, int N2);
  void forward( int idx, int N2);   // idx -> idx+1
  void backward(int idx, int N2);   // idx+1 -> idx
  void downward(int idx, int N2);   // st Pxout
  double CalcLmdY( int idx, int N2);               // calc lambda: I(Y)
  double CalcLmdXY(int idx, int N2, const int *X); // calc lambda: I(X,Y)
public:
  FBalg4(int _Nseg, int _Dmax, const char *fnGtable);
  ~FBalg4();
  int    getBeta();
  double getPid();
  double getPs();
  double getPth();
  void   calc(double **Pout, const double **Pin, const unsigned char *Yin, int N2);
  double calcIxy(const double **Pin, const unsigned char *Yin, int N2, const unsigned char *Xin, int nu, const int *B);
  void   dump(int idx);
  void   dump();
  void   dumpFB();
  void   dbgSetX(const unsigned char *X);
};
