class SLFBAdec {
private:
  int    Ns;            // block length (symbol)
  int    Nu;            // symbol size (bit)
  int    Nb;            // = Ns*Nu (bit)
  int    Q;             // number of inner codewords |C|
  int    Dmax,Dmin;     // channel
  int    Drng;          // = Dmax-Dmin+1
  int    Nu2max,Nu2min; // recv symbol length
  long   Nu2p;          // = pow(2,Nu)
  double Pi,Pd,Ps;      // channel
  double **GD;          // [Drng][Drng]: P(d_Nu | d_0)
  double ***GX;         // [Nu2][y][x]:  P(y|x,d0,d1)
  class InnerCodebook *ICB;
  class ChannelMatrix *ECM;   // (encoding CM)
  class IDSchannel    *CH;
  //-----
  int    max(int a, int b);
  int    min(int a, int b);
  double Psub(unsigned char a, unsigned char b);  // substitution prob
  void ClearMat(double **X, int M0, int N0);
  void CopyMat(double **Y, const double **X, int M0, int N0);  // Y=X
  void PrintMat(const double **X, int M0, int N0, const char *pre, const char *post);
  void MultMat(double **Y, const double **X0, const double **X1, int M0, int N0, int M1, int N1); // Y = X0*X1
  void PrintVect(const unsigned char *V, int len, const char *pre, const char *post);
  void PrintVect(const int           *V, int len, const char *pre, const char *post);
  long VectToLong(const unsigned char *V, int len);
  void LongToVect(unsigned char *V, long val, int len);
  //-----
  void SetGD();
  void DelGD();
  void SetGX();
  void DelGX();
  double CalcPyx(long y, long x, int ly, int lx);  // P(y|x)
public:
  SLFBAdec(class InnerCodebook *_ICB, class ChannelMatrix *_ECM, class IDSchannel *_CH);
  ~SLFBAdec();
};
