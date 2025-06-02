#include <list>
#include <algorithm>

struct Gtable_xp {
  int   x;
  float p;
  //-----
  bool operator<( const Gtable_xp& right ) const {
    return p==right.p ? x < right.x : p > right.p;
    //return p > right.p;
  }
};

class Gtable {
private:
  int    beta, beta4p, beta4pp;  // beta4p=4^beta beta4pp=4^2beta
  double Pid,Ps,Ptrn;            // Ptrn = 1-2*Pid
  double Pth;                    // cutoff threshold
  //-----
  double Psub(unsigned char x, unsigned char y);
  double CalcProb(const unsigned char *X, int alphaX, int betaX,
		  const unsigned char *Y, int alphaY, int betaY);
public:
  std::list<Gtable_xp> **LIST;   // list of (x,prb): [beta':0...2beta][y:0...4^beta']
  //-----
  Gtable(int _beta, double _Pid, double _Ps, double _Pth); // empty
  Gtable(const char *fn); // read from file
  ~Gtable();
  void   SetLIST();
  void   write(const char *fn);
  int    getBeta();
  double getPid();
  double getPs();
  double getPth();
  void   dump();
};
