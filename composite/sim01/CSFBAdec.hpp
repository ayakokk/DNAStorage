class CSFBAdec {
private:
  int    Ns;            // block length (symbol)
  int    Nu;            // symbol size (bit)
  int    Nb;            // = Ns*Nu (bit)
  int    Nseq;          // Number of sequences (=r)
  int    Q;             // number of inner codewords |C|
  int    Dmax,Dmin;     // channel
  int    Drng;          // = Dmax-Dmin+1
  int    Nu2max,Nu2min; // recv symbol length
  int    NumItr;        // number of decoding iterations
  long   Nu2p;          // = pow(2,Nu)
  double Pi,Pd,Ps,Pt;   // channel
  double **GD;          // [Drng][Drng]: P(d_Nu | d_0)
  double ***GX;         // [Nu2][y][x]:  P(y|x,d0,d1)
  double ***PF, ***PB;  // [Nseq][Ns+1][Drng]: FG forward/backward
  double ***PU, ***PD;  // [Nseq][Ns][Q]:      FG up/down
  double **PI, **PO;    // [Ns][Q]:            FG input/output
  double **PU0;         // [Ns][Q]:            FG PI->CM->PU0
  double **PM;          // [Ns][Q]:            FG memory
  unsigned char **Yin;  // [Nseq][Nb+Dmax]: received word
  class InnerCodebook *ICB;
  class ChannelMatrix *ECM;   // (encoding CM)
  class IDSchannel    *CH;
  //-----
  int    max(int a, int b);
  int    min(int a, int b);
  int    max(const int *a, int len);
  int    min(const int *a, int len);
  double Psub(unsigned char a, unsigned char b);  // substitution prob
  void ClearMat(double **X,  int M0, int N0);
  void ClearMat(double ***X, int M0, int N0, int K0);
  void SetMat(double **X, int M0, int N0, double val);
  void CopyMat(double **Y, const double **X, int M0, int N0);  // Y=X
  void PrintMat(const double **X, int M0, int N0, const char *pre, const char *post);
  void MultMat(double **Y, const double **X0, const double **X1, int M0, int N0, int M1, int N1); // Y = X0*X1
  void PrintVect(const unsigned char *V, int len, const char *pre, const char *post);
  void PrintVect(const int           *V, int len, const char *pre, const char *post);
  void PrintVect(const double        *V, int len, const char *pre, const char *post);
  long VectToLong(const unsigned char *V, int len);
  void LongToVect(unsigned char *V, long val, int len);
  void normalize(double *V, int len);
  //-----
  void SetGD();
  void DelGD();
  void SetGX();
  void DelGX();
  double CalcPyx(long y, long x, int ly, int lx);  // P(y|x)
  double GetGX(int Nu2, long y, long xi);
  void SetFG();
  void DelFG();
  void InitFG(const unsigned char **RW, const double **Pin, int *Nb2);
  void CalcPU0(int idx);                    // PU0 (input->up)
  void CalcPU( int idx);                    // PU
  void CalcPM( int idx);                    // PM
  void CalcPO( int idx);                    // PO
  void CalcPF( int jdx, int idx, int Nb2);  // PF
  void CalcPB( int jdx, int idx, int Nb2);  // PB
  void CalcPD( int jdx, int idx, int Nb2);  // PD
  void CalcPDf(int jdx, int idx, int Nb2);  // PD (forward)
  void CalcPDb(int jdx, int idx, int Nb2);  // PD (backward)
public:
  CSFBAdec(class InnerCodebook *_ICB, class ChannelMatrix *_ECM, class IDSchannel *_CH, int _Nseq, int _NumItr);
  ~CSFBAdec();
  void Decode(double **Pout, const unsigned char **RW, int *Nb2, const int *dbgIW, const double **Pin);
  void Decode(double **Pout, const unsigned char **RW, int *Nb2, const int *dbgIW);  // (uniform prior)
  void PrintNode(int jdx, int idx, int iw);
  void PrintNode(const int *dbgIW);
};
