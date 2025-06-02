class MetricCalc {
private:
  int    beta, beta2p;
  double Pi, Pd, Ps, Pt;
  double *Qi;   // insertion  [0,1,2]
  double *Qd;   // deleteion  [0] 
  double *Qs;   // substition [0,1]
  double **L;   // calculation lattice [beta+1][beta+1]
  double **MT;  // metric [beta2p][beta2p]
  //-----
  int  min(int a, int b);
  int  HammingDist(unsigned char v0, unsigned char v1);  // one symbol
  int  HammingDist(const unsigned char *V0, const unsigned char *V1, int len);
  void PrintVect(unsigned char *V, int len, const char *pre, const char *post);
  bool CheckConv();
  void CalcLval(int s, int t, const unsigned char *V0, const unsigned char *V1);
  void SetD();
  void SetMT();
  void DumpL();
  void DumpMT();
public:
  MetricCalc(int _beta, double _Pi, double _Pd, double _Ps);
  ~MetricCalc();
  double calc( const unsigned char *V0, const unsigned char *V1);
  double calc1(const unsigned char *V0, const unsigned char *V1);
  void   ConvUintVect(unsigned char *V, unsigned int u, unsigned int len);  // V <- u
  unsigned int ConvVectUint(const unsigned char *V,   unsigned int len);  // u <- V
};
