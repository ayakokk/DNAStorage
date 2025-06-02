class ConcatAB {
private:
  int  Nb, Nu, Ns;          // code length Nb = Nu * Ns
  int  Nseq;                // Num sequences
  int  Qa, Qb;              // Number symbols
  int  Dmin, Dmax, Drng;    // drift
  double ***PrD;            // [Nseq][Ns+1][Drng]: Pr(drift)
  double **PrB2;            // [Nb][2]:  Pr(b_i) binary
  double **PrBQ;            // [Ns][Qa]: Pr(b)   Q-ary
  class InnerCodebook *ICB;
  class CSFBAdec      *DECa;
  class CPScodec      *DECb;
  //----
  void ClearVect(double *V, int len);
  void Normalize(double *V, int len);
  void ConvPrbQB();
public:
  ConcatAB(class InnerCodebook *_ICB, class CSFBAdec *_DECa, class CPScodec *_DECb);
  ~ConcatAB();
  void AtoB();  // DECa -> DECb
  void Get_PrB2(double **PrB2out);
};

