class CodeParam{
private:
  int beta,nu;
  int *B;
public:
  CodeParam(const char *fn);
  ~CodeParam();
  int get_beta();
  int get_nu();
  void get_B(int *_B);
};
