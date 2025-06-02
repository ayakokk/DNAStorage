class MutualInfo {
private:
  int M, N;             // M-input N-output
  unsigned long **cnt;  // [M][N]: count
  unsigned long *cntX;  // [M]: X count
  unsigned long *cntY;  // [N]: Y count
  //-----
  void   PrintVect(const unsigned long *V, int len, const char *pre, const char *post);
  unsigned long sum(const unsigned long *V, int len);
  double entropy(const unsigned long *V, int len);
public:
  MutualInfo(int _M, int _N);
  ~MutualInfo();
  void   clear();
  void   countup(int x, int y);
  void   PrintCnt();
  double Hx();     // H(X)
  double Hxy();    // H(X|Y)
  double Ixy();    // I(X;Y)
  void   GetPy(double *Py);
};
