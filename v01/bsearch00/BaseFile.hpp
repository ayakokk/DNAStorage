class BaseFile {
private:
  int beta, nu;
  int *B;   // [nu] bit/block
  int *B2p; // [nu] pow(2,B[i]-beta)
  int **CW;  // [nu][B2p[i]]
  //---
  void PrintVect(const int *V, int len, const char *pre, const char *post);
public:
  BaseFile(const char *fn);
  BaseFile(int _beta, int _nu, const int *_B);  // CW empty 
  ~BaseFile();
  int  get_beta();
  int  get_nu();
  void get_B(  int *_B);
  void get_B2p(int *_B2p);
  void get_CW(int idx, int *_CW);
  void set_CW(int idx, const int *_CW);
  void write(const char *fn);
};
