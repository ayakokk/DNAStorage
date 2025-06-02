class SLFBAdec {
private:
  int    Ns;           // block length (symbol)
  int    Nu,NuA,NuB;   // symbol size (bit), 2-seg lengths
  int    Nb;           // = Ns*Nu (bit)
  int    Q;            // number of inner codewords |C|
  int    Dmax,Dmin;    // channel
  int    Drng;         // = Dmax-Dmin+1
  double Pi,Pd,Ps;     // channel
  double **GDA,**GDB;  // [Drng][Drng]: P(d_Nu{A,B} | d_0)
  class InnerCodebook *ICB;
  class ChannelMatrix *ECM;   // (encoding CM)
  class IDSchannel    *CH;
  //-----
  void ClearMat(double **X, int M0, int N0);
  void CopyMat(double **Y, const double **X, int M0, int N0);  // Y=X
  void PrintMat(const double **X, int M0, int N0, const char *pre, const char *post);
  void MultMat(double **Y, const double **X0, const double **X1, int M0, int N0, int M1, int N1); // Y = X0*X1
  long VectToLong(const unsigned char *V, int len);
  void LongToVect(unsigned char *V, long val, int len);
  //-----
  void SetGD();
  void DelGD();
public:
  SLFBAdec(class InnerCodebook *_ICB, class ChannelMatrix *_ECM, class IDSchannel *_CH);
  ~SLFBAdec();
};
