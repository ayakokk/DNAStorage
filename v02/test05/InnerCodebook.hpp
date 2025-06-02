class InnerCodebook{
private:
  int  Nu;               // block length
  long Nu2p;             // 2^Nu
  long numCW;            // number of codewords
  bool FlgUnique;        // CW: uniqueness
  bool FlgInvertible;    // CW: invertible
  bool FlgBalanced;      // CW: balanced
  int  Rho;              // run-length
  int  ell,Delta;        // local-balance
  int  L0;               // ceil( (ell-1)/Nu )
  unsigned char **CW;    // [numCW][Nu]: codewords
  int  *RLmax;           // [numCW]: max run-length (intra-codeword)
  int  *RLleft;          // [numCW]: left-end RL
  int  *RLright;         // [numCW]: right-end RL
  int  **Wleft;          // [numCW][Nu]: weight left
  int  **Wright;         // [numCW][Ni]: weight right
  int  *CWL;             // [Nu2p]: CW list
  int  Rst01, Rst10;     // Reset symbol: CW[Rst01]=0101... CW[Rst10]=1010...
  //-----
  int  argmax(const int *V, int len);
  int  argmin(const int *V, int len);
  int  max(const int *V, int len);
  int  min(const int *V, int len);
  void VectInv(unsigned char *VI, const unsigned char *V, int len);  // binary vector invert 
  long VectToLong(const unsigned char *V, int len);
  void LongToVect(unsigned char *V, long val, int len);
  bool CheckVectConv();                                    // VectToLong -> LongToVect
  void PrintVect(const unsigned char *V, int len, const char *pre, const char *post);
  void PrintVect(const int           *V, int len, const char *pre, const char *post);
  bool IsEqual(   const unsigned char *V0, const unsigned char *V1, int len);  // (V0==V1)?
  bool IsInvEqual(const unsigned char *V0, const unsigned char *V1, int len);  // (V0==~V1)?
  int  Balance01(const unsigned char *V, int len);         // -Nu(all-0) to +Nu(all-1)
  int  BalanceSW(const unsigned char *V, int len, int ws); // sliding window (ws:window size)
  int  MaxRunLen(const unsigned char *V, int len);
  int  RunLenF(const unsigned char *V, int pos, int len);  // V[pos] => V[len-1]
  int  RunLenB(const unsigned char *V, int pos, int len);  // V[0] <= V[pos]
  void ReadFile(const char *fn);
  void SetFlg();         // FlgUnique FlgInvertible
  void SetCWL();         // CW -> CWL
  void SetWLR(const unsigned char *V, int *WL, int *WR, int len);  // left/right weight
  void SetRstSymb();                                               // reset symbol (or -1)
  bool EncodeRunLen( const unsigned char *V, int idx, int len);    // idx,len: in symbol
  bool EncodeBalance(const unsigned char *V, int idx, int len);    // sliding window: idx->idx+Nu-1
public:
  InnerCodebook(const char *fn, int _Rho, int _ell, int _Delta);
  ~InnerCodebook();
  void Encode(unsigned char *CV, const int *IV, int Ns);   // Ns: information vector length (symbol)
  void PrintCodebook();
  int  Get_Nu();
  int  Get_Nu2p();
  int  Get_numCW();
  bool Get_FlgUnique();
  bool Get_FlgInvertible();
};
