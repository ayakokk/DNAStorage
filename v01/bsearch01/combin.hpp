#include <boost/math/special_functions/binomial.hpp>

class combin{
private:
  int N,K;
  unsigned long Num;
  unsigned char *V;   // weight-k pattern
  int HammingWeight(const unsigned char *X, int len);
public:
  combin(unsigned int _N, unsigned int _K);
  ~combin();
  unsigned long getNum();
  void initV();
  bool nextV();
  void getV(unsigned char *X);
  void printV();
};
