// A=0 T=1 G=2 C=3
class InnerCodebook{
private:
  int beta,nu;
  int beta4p;
  int *B, *B2p;        // [nu]
  class bmatrix *CWM;  // codeword map [beta4p][nu]
  int **EncMap;        // encode map [B2p][nu]
  int **DecMap;        // decode map [beta4p][nu]
  //-----
  void FuncTest();
  void PrintParam();
public:
  InnerCodebook(int _beta, int *_B, int _nu);    // empty
  InnerCodebook(const char *fn);                 // read from file
  InnerCodebook(class InnerCodebook *ICB); // copy
  ~InnerCodebook();
  int  Get_beta();
  int  Get_nu();
  int  Get_beta4p();
  void Get_B(int *_B);
  void Get_B2p(int *_B2p);
  int  GetNumCW(int cb);
  void GetNumCW(int *NumCW);
  void Encode(unsigned char *ICout, const int *ICin, int Nseg);
  void Decode(int *IDout, const unsigned char *IDin, int Nseg);
  void GenRndInfo(int *ICin, int Nseg);
  void GetPrior(double **Px, int Nseg);
  int  CWMget(int idx, int cb);
  void CWMset(int idx, int cb, int val);
  void CWMset(class bmatrix *CW, int cb);
  void CWMclear();
  void CWMset();  // all
  void CWMwrite(const char *fn);
  void GenMap();  // generate EncMap, DecMap
  bool check();
  void PrintCWM();
  void dump();
};

