class InnerCodec{
private:
  int  ell, rlk, eps;  // ell:window size / rlk; max run-length / eps: local-balance
  long numQ;           // number of states: pow(4,ell)
  int  **STT;          // [numQ][4]: state transition table
  bool *Vflg;          // [numQ]: valid state
  //-----
  void ConvLongStat(int *Qx, long val, int len);   // S <= state(val)
  long ConvStatLong(const int *Qx, int len);       // returns LongVal(S)
  void ConvCheck(int len);
  int  MaxRunLen(const int *Qx, int len);          // max run-length in Qx
  int  LocalBalance(const int *Qx, int len);       // num23-num01
  long NextState(long Q, int a);                  // t(Q,a)
  void SetTable(); // STT, Vflg
  //-----
  void PrintVect(const int *V, int len, const char *pre, const char *post);
  void SetVect(bool *V, long len, bool val);
  void SetVect(int  *V, long len, int  val);
public:
  InnerCodec(int _ell, int _rlk, int _eps);
  ~InnerCodec();
};
