// A=0 T=1 G=2 C=3
class InnerCodebook{
private:
  int beta,B,nu;
  int beta4p,B2p;
  class bmatrix *CWM;  // codeword map [beta4p][nu]
  int **EncMap;        // encode map [B2p][nu]
  int **DecMap;        // decode map [beta4p][nu]
  //-----
  int  GetNumCW(int cb);
  void FuncTest();
public:
  InnerCodebook(int _beta, int _B, int _nu);  // empty
  InnerCodebook(const char *fn);      // read from file
  ~InnerCodebook();
  int  Get_beta();
  int  Get_B();
  int  Get_nu();
  int  Get_beta4p();
  int  Get_B2p();
  void Encode(unsigned char *ICout, const int *ICin, int Nseg);
  void Decode(int *IDout, const unsigned char *IDin, int Nseg);
  void GetPrior(double **Px, int Nseg);
  int  CWMget(int idx, int cb);
  void CWMset(int idx, int cb, int val);
  void CWMset(class bmatrix *CW, int cb);
  void CWMwrite(const char *fn);
  void GenMap();  // generate EncMap, DecMap
  bool check();
  void dump();
};

