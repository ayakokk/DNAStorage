class CIDSchannel {
private:
  int    N;           // block length
  int    Kc, Nseq;    // composite resolution & Num. seq.
  int    Dmin,Dmax;   // drift range
  double Pi,Pd,Ps;    // error prob
  unsigned char **XV; // [Nseq][N]: channel input
  int    **DV;        // [Nseq][N+1]: drift vector
  //-----
  int max(int a, int b);
  int min(int a, int b);
  int sum(const int *V, int len);
  unsigned char inv(unsigned char a);
  unsigned char CompSel(const int *dist);  // select according to dist[4]
public:
  CIDSchannel(int _N, double _Pi, double _Pd, double _Ps, int _Dmin, int _Dmax, int _Kc, int _Nseq);
  ~CIDSchannel();
  void   transmit(unsigned char **YV, int *N2, const int **XC);           // YV[Nseq][N+Dmax], N2[Nseq], XC[N][4]  
  int    transmitOne(unsigned char *Y, int *DD, const unsigned char *X);  // returns length of Y
  void   funcFd(unsigned char **YV2, const unsigned char **YV);           // YV2 <- YV:  binarize [Nseq][N+Dmax]
  int    GetN();
  int    GetDmin();
  int    GetDmax();
  double GetPi();
  double GetPd();
  double GetPs();
  void   GetDV(int           **dbgDV);  //(dbg) [Nseq][N+1]
  void   GetXV(unsigned char **dbgXV);  //(dbg) [Nseq][N]
};
