class NBIDSchannel{
private:
  int    N,Q,Dmax;
  double Pid,Ps;
  int    *DV;            // drift vector [N+1]
  //-----
  void ClearVect(unsigned char *V, int len);
  int  Rnd01(double p1);  // p(1)=p1, p(0)=1-p1
  int  NextD(int d);      // next drift value
  unsigned char substitution(unsigned char z); // substitution error 
public:
  NBIDSchannel(int _N, int _Q, int _Dmax, double _Pid, double _Ps);
  ~NBIDSchannel();
  int transmit(unsigned char *Y, const unsigned char *X); // returns length of Y
};
