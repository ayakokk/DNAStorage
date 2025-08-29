class SLFBAdec {
private:
  // エラー状態定数
  static const int ERROR_MATCH = 0;
  static const int ERROR_INSERTION = 1;
  static const int ERROR_DELETION = 2;
  static const int ERROR_SUBSTITUTION = 3;
  static const int NUM_ERROR_STATES = 4;
  
  int    Ns;            // block length (symbol)
  int    Nu;            // symbol size (bit)
  int    Nb;            // = Ns*Nu (bit)
  int    Q;             // number of inner codewords |C|
  int    Dmax,Dmin;     // channel
  int    Drng;          // = Dmax-Dmin+1
  int    Nu2max,Nu2min; // recv symbol length
  long   Nu2p;          // = pow(2,Nu)
  double Pi,Pd,Ps,Pt;   // channel
  double **GD;          // [Drng][Drng]: P(d_Nu | d_0)
  double ***GX;         // [Nu2][y][x]:  P(y|x,d0,d1)
  double **PF, **PB;    // [Ns+1][Drng]: FG forward/backward
  double **PU, **PD;    // [Ns][Q]:      FG up/down
  double **PI, **PO;    // [Ns][Q]:      FG input/output
  unsigned char *Yin;   // [Nb+Dmax]: received word
  
  // エラー状態関連データ構造
  double ***PFE, ***PBE; // [Ns+1][Drng][NUM_ERROR_STATES]: エラー状態を含む前進/後進確率
  double ***PE;          // [Ns][NUM_ERROR_STATES][NUM_ERROR_STATES]: エラー状態遷移確率
  int k_mer_length;      // k-merの長さ（将来の拡張用）
  class InnerCodebook *ICB;
  class ChannelMatrix *ECM;   // (encoding CM)
  class IDSchannel    *CH;
  //-----
  int    max(int a, int b);
  int    min(int a, int b);
  double Psub(unsigned char a, unsigned char b);  // substitution prob
  void ClearMat(double **X, int M0, int N0);
  void CopyMat(double **Y, const double **X, int M0, int N0);  // Y=X
  void PrintMat(const double **X, int M0, int N0, const char *pre, const char *post);
  void MultMat(double **Y, const double **X0, const double **X1, int M0, int N0, int M1, int N1); // Y = X0*X1
  void PrintVect(const unsigned char *V, int len, const char *pre, const char *post);
  void PrintVect(const int           *V, int len, const char *pre, const char *post);
  void PrintVect(const double        *V, int len, const char *pre, const char *post);
  long VectToLong(const unsigned char *V, int len);
  void LongToVect(unsigned char *V, long val, int len);
  void normalize(double *V, int len);
  //-----
  void SetGD();
  void DelGD();
  void SetGX();
  void DelGX();
  double CalcPyx(long y, long x, int ly, int lx);  // P(y|x)
  double GetGX(int Nu2, long y, long xi);
  void SetFG();
  void DelFG();
  void InitFG(const unsigned char *RW, const double **Pin, int Nb2);
  void CalcPU(int idx);           // PU
  void CalcPO(int idx);           // PO  
  void CalcPF(int idx, int Nb2);  // PF
  void CalcPB(int idx, int Nb2);  // PB
  void CalcPD(int idx, int Nb2);  // PD
  
  // エラー状態関連関数
  void SetFGE();                  // エラー状態拡張FGデータ構造の設定
  void DelFGE();                  // エラー状態拡張FGデータ構造の削除
  void InitFGE(const unsigned char *RW, const double **Pin, int Nb2); // エラー状態拡張FGの初期化
  void CalcPE(int idx);           // エラー状態遷移確率の計算
  void CalcPFE(int idx, int Nb2); // エラー状態を含む前進確率
  void CalcPBE(int idx, int Nb2); // エラー状態を含む後進確率
  void CalcPDE(int idx, int Nb2); // エラー状態を含む事後確率計算
public:
  SLFBAdec(class InnerCodebook *_ICB, class ChannelMatrix *_ECM, class IDSchannel *_CH);
  ~SLFBAdec();
  void Decode(double **Pout, const unsigned char *RW, int Nb2, const int *dbgIW, const double **Pin);
  void Decode(double **Pout, const unsigned char *RW, int Nb2, const int *dbgIW);  // (uniform prior)
  void PrintNode(int idx, int iw);
  void PrintNode(const int *dbgIW);
  // テーブル出力機能
  void exportGDTable(const char* filename);
  void exportGXTable(const char* filename);
  void exportAllTables(const char* output_dir);
};
