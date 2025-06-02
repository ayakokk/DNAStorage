class CPScodec {
private:
  int    Ns;            // block length (symbol)
  int    Nu;            // symbol size (bit)
  int    Nb;            // = Ns*Nu (bit)
  int    Nseq;          // Number of sequences (=r)
  int    Kc;            // resolution (=k)
  int    Dmax,Dmin;     // channel
  int    Drng;          // = Dmax-Dmin+1
  int    QB, QB2;       // number of composite symbols, QB2=QB/2
  int    **CSL;         // [QB][4]: composite symbol list
  double Pi,Pd,Ps,Pt;   // channel
  double **GS1;         // [4][QB]:    p(y|s)*Pt     (trans)
  double ***GS2;        // [4][4][QB]: p(y0,y1|s)*Pi (ins)
  double ***PF, ***PB;  // [Nseq][Nb+1][Drng]: FG forward/backward
  double ***PU, ***PD;  // [Nseq][Nb][QB]:     FG up/down
  double **PM;          // [Nb][QB]:           FG memory
  double **PXU,**PXD;   // [Nb][QB]:           FG green
  double **PIB;         // [Nb][2]:            FG b_i input
  double **PIC, **POC;  // [Nb][QB2]:          FG C_i input/output
  unsigned char **Yin;  // [Nseq][Nb+Dmax]: received word
  //-----
  int    max(int a, int b);
  int    min(int a, int b);
  int    max(const int *a, int len);
  int    min(const int *a, int len);
  double Psub(unsigned char a, unsigned char b);  // substitution prob
  void ClearVect(double  *X, int M0);
  void ClearMat(double  **X, int M0, int N0);
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
  void   SetCSL(const char *fn);
  void   DelCSL();
  void   SetGS();
  void   DelGS();
  void   SetFG();
  void   DelFG();
  double Psub(int x, int y);
  double Psbc(int s, int b, int c);
  void   CalcPXU(int idx);
  void   CalcPOC(int idx);
  void   CalcPDf(int jdx, int idx, int Nb2);
  void   CalcPDb(int jdx, int idx, int Nb2);
  void   CalcPD( int jdx, int idx, int Nb2);
  void   CalcPF( int jdx, int idx, int Nb2);
  void   CalcPB( int jdx, int idx, int Nb2);
  void   CalcPU( int idx);
  void   CalcPXD(int idx);
public:
  CPScodec(int _Nb, double _Pi, double _Pd, double _Ps, int _Dmin, int _Dmax,
	   int _Nu, int _Nseq, const char *fn);
  ~CPScodec();
  void   Encode(int **CCW, const int **IWB, const unsigned char *CWA); // CCW[Nb][4]: composite CW
  void   Decode(double **Pout, const unsigned char **RW, int *Nb2, const int **dbgIW, const double **Pin);
  void   Decode(double **Pout, const unsigned char **RW, int *Nb2, const int **dbgIW);  // (uniform prior)
  void   InitFG();
  void   InitFG(const unsigned char **RW, const double **Pin, int *Nb2);
  void   SetFGD(const double ***PrD);  // PF,PB <- Pin[Nseq][Ns+1][Drng]
  void   SetFGB(const double **PrB2);  // PIB <- PrB2[Nb][2]
  void   PrintCSL();
  void   PrintGS();
  void   PrintFG(int idx, int dbgIWi);
  void   PrintFG(int idx);
  void   PrintFG(const int **dbgIW);
  void   PrintFG();
  int    Get_Ns();
  int    Get_Nu();
  int    Get_Nb();
  int    Get_Nseq();
  int    Get_Kc();
  int    Get_Dmin();
  int    Get_Dmax();
  int    Get_Drng();
  int    Get_QB();
  double Get_Pi();
  double Get_Pd();
  double Get_Ps();
  double Get_Pt();
};
