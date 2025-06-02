class InnerCodec {
private:
  int    RunK, Delta;   // RL, balance
  double Pi, Pd, Ps;    // IDS channel
  int    Dmin,Dmax;     // drift max & min
  int    N;             // block lenth
  int    Nmax;          // = N+Dmax: max block length
  int    NumSTR;        // = 4*RunK+1
  int    NumSTB;        // = 2*Delta+2
  int    NumD;          // = Dmax-Dmin+1;
  int    **STR;         // [4*RunK+1][4]:  state transition (run length)
  int    **STB;         // [2*Delta+2][4]: state transition (balance)
  double ****GT1;       // [Y:4]       [X:4][P:NumSTR][Q:NumSTB]: G-table(L=1)
  double *****GT2;      // [Y0:4][Y1:4][X:4][P:NumSTR][Q:NumSTB]: G-table(L=2)
  double **PU, **PD;    // [N][4]: factor graph Xi up/down
  double ****PF;        // [N+1][NumSTR][NumSTB][NumD]: factor graph forward
  double ****PB;        // [N+1][NumSTR][NumSTB][NumD]: factor graph forward
  int    *Y;            // [Nmax]: received wird
  //-----
  int    BitGet(int  val, int pos);
  void   BitSet(int *val, int pos, int sval);
  void   Normalize(double *V, int len);
  void   Normalize(double ***V, int len0, int len1, int len2);
  void   ClearVect(double *V, int len);
  void   PrintVect(const int    *V, int len, const char *pre, const char *post);
  void   PrintVect(const double *V, int len, const char *pre, const char *post);
  double SumVect(const double *V, int len);
  //-----
  int    CalcL(int p, int q);
  double FuncPd(int d0, int d1);  // P(d1|d0)
  double FuncPpqx(int p0, int q0, int x, int p1, int q1);
  void   SetSTR();
  void   SetSTB();
  void   SetGT();
  void   AllocFG();
  void   DelFG();
  void   InitFG(const int *RW, const double **Px, int N2);
  void   ClearNode(int idx);
  void   ClearNode();
  void   Forward(int idx, int N2);
public:
  InnerCodec(int _RunK, int _Delta, double _Pi, double _Pd, double _Ps, int _N, int _Dmin, int _Dmax);
  ~InnerCodec();
  void Encode(int *CW, const int *IW);
  void Decode(double **Pout, const int *RW, const double **Px, int N2, const int *DbgCW);
  void PrintSTR();
  void PrintSTB();
  void PrintFenc();
  void PrintGT1();
  void PrintGT2();
  void PrintNode(int idx, int DbgCW);
  void PrintNode(const int *DbgCW);
};
