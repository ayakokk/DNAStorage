class ChannelMatrix {
private:
  int M, N;             // M-input N-output
  unsigned long **cnt;  // [M][N]: count
  unsigned long *cntX;  // [M]: X count
  unsigned long *cntY;  // [N]: Y count
  unsigned long cntAll; // all count
  double **Pxy;         // [M][N]: cnt/cntAll
  //-----
  void   PrintVect(const unsigned long *V, int len, const char *pre, const char *post);
  unsigned long sum(const unsigned long *V, int len);
  double entropy(const unsigned long *V, int len);
public:
  ChannelMatrix(int _M, int _N);
  ChannelMatrix(const char *fnPxy);
  ~ChannelMatrix();
  void   clear();
  void   countup(int x, int y);
  void   PrintCnt();
  void   PrintPxy();
  void   WritePxy(const char *fn);
  double Hx();     // H(X)
  double Hxy();    // H(X|Y)
  double Ixy();    // I(X;Y)
  void   GetPy(double *Py);
  int    GetM();
  int    GetN();
  double GetPxy(int x, int y);
};

bool ChannelMatrix_IsEqual(class ChannelMatrix *CM0, class ChannelMatrix *CM1); // compares Pxy
