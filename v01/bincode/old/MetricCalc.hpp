class MetricCalc {
private:
  int beta, beta2p;
  int **L;   // calculation lattice [beta+1][beta+1]
  int **MT;  // metric [beta2p][beta2p]
  //-----
  int  min(int a, int b);
  void PrintVect(unsigned char *V, int len, const char *pre, const char *post);
  bool CheckConv();
  void CalcLval(int s, int t, const unsigned char *V0, const unsigned char *V1);
  void SetMT();
  void DumpL();
  void DumpMT();
public:
  MetricCalc(int _beta);
  ~MetricCalc();
  int calc(const unsigned char *V0, const unsigned char *V1);
  void ConvUintVect(unsigned char *V, unsigned int u, unsigned int len);  // V <- u
  unsigned int ConvVectUint(const unsigned char *V,   unsigned int len);  // u <- V
};
