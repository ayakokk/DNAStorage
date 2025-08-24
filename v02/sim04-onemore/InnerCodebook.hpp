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
  int  *CWL;             // [Nu2p]: CW list (count)
  int  *CWI;             // [Nu2p]: CW->index(or -1) (inverse map)
  int  Rst01, Rst10;     // Reset symbol: CW[Rst01]=0101... CW[Rst10]=1010...
  //-----
  int  argmax(const int *V, int len);
  int  argmin(const int *V, int len);
  int  max(const int *V, int len);
  int  min(const int *V, int len);
  void VectInv(unsigned char *VI, const unsigned char *V, int len);  // binary vector invert
  void VectInv4(unsigned char *VI, const unsigned char *V, int len);  // 4元上位ビット反転
  void Extract4toUpper(unsigned char *upper, const unsigned char *V4, int len);  // 4元→上位ビット抽出 
  long VectToLong(const unsigned char *V, int len);
  void LongToVect(unsigned char *V, long val, int len);
  long VectToLong4(const unsigned char *V4, int len);  // 4元用: 上位ビット抽出してVectToLong
  bool CheckVectConv();                                    // VectToLong -> LongToVect
  void PrintVect(const unsigned char *V, int len, const char *pre, const char *post);
  void PrintVect(const int           *V, int len, const char *pre, const char *post);
  bool IsEqual(   const unsigned char *V0, const unsigned char *V1, int len);  // (V0==V1)?
  bool IsInvEqual(const unsigned char *V0, const unsigned char *V1, int len);  // (V0==~V1)?
  bool IsEqual4(   const unsigned char *V0, const unsigned char *V1, int len);  // 4元用: 上位ビット比較
  bool IsInvEqual4(const unsigned char *V0, const unsigned char *V1, int len);  // 4元用: 上位ビット反転比較
  int  Balance01(const unsigned char *V, int len);         // -Nu(all-0) to +Nu(all-1)
  int  Balance014(const unsigned char *V4, int len);       // 4元用: 上位ビットバランス計算
  int  BalanceSW(const unsigned char *V, int len, int ws); // sliding window (ws:window size)
  int  BalanceSW4(const unsigned char *V4, int len, int ws); // 4元用: 上位ビットスライディングウィンドウバランス
  int  MaxRunLen(const unsigned char *V, int len);
  int  MaxRunLen4(const unsigned char *V4, int len);      // 4元用: 上位ビットランレングス計算
  int  RunLenF(const unsigned char *V, int pos, int len);  // V[pos] => V[len-1]
  int  RunLenB(const unsigned char *V, int pos, int len);  // V[0] <= V[pos]
  int  RunLenF4(const unsigned char *V4, int pos, int len);  // 4元用: 上位ビット前方ランレングス
  int  RunLenB4(const unsigned char *V4, int pos, int len);  // 4元用: 上位ビット後方ランレングス
  void ReadFile(const char *fn);
  void SetFlg();         // FlgUnique FlgInvertible
  void SetCWL();         // CW -> CWL
  void SetCWI();         // CW -> CWI
  void SetWLR(const unsigned char *V, int *WL, int *WR, int len);  // left/right weight
  void SetWLR4(const unsigned char *V4, int *WL, int *WR, int len);  // 4元用: 上位ビット重み計算
  void SetRstSymb();                                               // reset symbol (or -1)
  bool EncodeRunLen( const unsigned char *V, int idx, int len);    // idx,len: in symbol
  bool EncodeBalance(const unsigned char *V, int idx, int len);    // sliding window: idx->idx+Nu-1
  bool EncodeRunLen4( const unsigned char *V4, int idx, int len);   // 4元用: 上位ビット制約チェック
  bool EncodeBalance4(const unsigned char *V4, int idx, int len);   // 4元用: 上位ビット制約チェック
public:
  InnerCodebook(const char *fn, int _Rho, int _ell, int _Delta);
  ~InnerCodebook();
  void Encode(unsigned char *CV, const int *IV, int Ns);   // Ns: information vector length (symbol)
  int  CWindex(const unsigned char *V);   // returns CWI[V]
  void PrintCodebook();
  void Get_CW(unsigned char *V, int idx);
  int  Get_Nu();
  int  Get_Nu2p();
  int  Get_numCW();
  bool Get_FlgUnique();
  bool Get_FlgInvertible();
};
