class EnumVect {
private:
  int N;   // vector length
  int *NE; // [N] number of elements
  int *V;  // [N] value
  //-----
  bool nextV(int pos);
  void PrintVect(const int *V, int len, const char *pre, const char *post);
public:
  EnumVect(int _N, const int *_NE);
  ~EnumVect();
  void initV();
  void getV(int *X);
  bool nextV();
};
