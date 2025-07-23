class IDSchannel {
private:
  int    N;          // block length
  int    Dmin,Dmax;  // drift range
  double Pi,Pd,Ps;   // error prob
  int    *DR;        // [N+1]: drift vector
  //-----
  int max(int a, int b);
  int min(int a, int b);
  unsigned char inv(unsigned char a);
public:
  IDSchannel(int _N, double _Pi, double _Pd, double _Ps);
  ~IDSchannel();
  int transmit(unsigned char *Y, const unsigned char *X);  // returns length of Y
  int    GetN();
  int    GetDmin();
  int    GetDmax();
  double GetPi();
  double GetPd();
  double GetPs();
  void   GetDR(int *dbgDR);  //(dbg)
};
