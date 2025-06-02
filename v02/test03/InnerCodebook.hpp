class InnerCodebook{
private:
  int N;                 // block length
  int numCW;             // number of codewords
  bool FlgUnique;        // CW uniqueness
  bool FlgInvertible;    // CW invertible 
  unsigned char **CW;    // [numCW][N]: codewords
  //-----
  void PrintVect(const unsigned char *V, int len, const char *pre, const char *post);
  bool IsEqual(   const unsigned char *V0, const unsigned char *V1, int len);  // (V0==V1)?
  bool IsInvEqual(const unsigned char *V0, const unsigned char *V1, int len);  // (V0==~V1)?  
  void ReadFile(const char *fn);
  void SetFlg();         // FlgUnique FlgInvertible
public:
  InnerCodebook(const char *fn);
  ~InnerCodebook();
  void PrintCodebook();
  bool GetFlgUnique();
  bool GetFlgInvertible();
};
