class bmatrix{
private:
  int M, N;        // bits
  int Nb,Nbs;      // Nb=ceil(N/8), Nbs=sizeof(unsigned char)*Nb
  int Nfrac;       // N-(Nb-1)*8: length of last byte
  unsigned char **V;  // [Mb][Nb]: elements
  //-----
  void GetSize(int *Mx, int *Nx, const char *fn, bool bin);
  void bitset(unsigned char *u, int pos, int s);
  int  bitget(unsigned char u, int pos);
  bool bittest();
  int  byteweight(unsigned char u, int len);
public:
  bmatrix(int _M, int _N);
  bmatrix(const char *fn, bool bin); // bin[true=binary,false=alist]
  ~bmatrix();
  int  getM();
  int  getN();
  int  getNb();
  void clear();
  void randmat();
  void randvec(int w);  // random row vector of weight w
  void print();
  int  getV(int r, int c);
  void setV(int r, int c, int s);
  void setV(class bmatrix *H);
  void setV(class bmatrix *H, int rt, int ct, int rs, int cs, int lr, int lc); // V[rt][ct] = H->V[rs][cs] (lr*lc mat)
  void invV(int r, int c); // invert
  void getRV(unsigned char *RV, int r);        // RV = V[r]
  void getRV(class bmatrix *RV, int r);        // RV = V[r]
  void setRV(int r, const unsigned char *RV);  // V[r] = RV
  void setRV(int r, class bmatrix *RV);        // V[r] = RV
  void setIdentity(int r, int c, int l);       // V[r][c] = I_l (l*l identity matrix)
  void RowExchange(int r0, int r1);
  void RowAdd(int r0, int r1);                 // V[r0]=V[r0]+V[r1]
  void RowAdd(int r, const unsigned char *RV); // V[r]=V[r]+RV
  void ColExchange(int c0, int c1);
  void GenCyclicMat(class bmatrix *RV);        // generate cyclic matrix: RV=1st row
  void ShiftRight(int r);                      // right cyclic shift of r-th row
  void ShiftLeft(int r);                       // left cyclic shift of r-th row
  bool isIdentity();
  bool isZero();
  bool isRcyclic();
  bool isLcyclic();
  int  HammingWeight();
  int  HammingWeightOld();
  void setAnd(class bmatrix *H0, class bmatrix *H1);  // V = H0 & H1
  void write(const char *fn);
  void read( const char *fn);
  void write_alist(const char *fn);
  void read_alist( const char *fn);
};

int  bmatrix_compare(class bmatrix *H0, class bmatrix *H1);   // 0:H0=H1 1:H0!=H1
bool bmatrix_inverse(class bmatrix *HI, class bmatrix *H);    // HI = H^{-1}  (returns false if not full rank)
void bmatrix_transpose(class bmatrix *HT, class bmatrix *H);  // HT = H^T
void bmatrix_mul(class bmatrix *H, class bmatrix *H0, class bmatrix *H1); // H = H0*H1
void bmatrix_add(class bmatrix *H, class bmatrix *H0, class bmatrix *H1); // H = H0+H1
void bmatrix_and(class bmatrix *H, class bmatrix *H0, class bmatrix *H1); // H = H0&H1
int  bmatrix_HammingDist(class bmatrix *H0, class bmatrix *H1); // Hamming distance
