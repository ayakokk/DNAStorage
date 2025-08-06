#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <set>

class SLFBAdec {
private:
  int    Ns;            // block length (symbol)
  int    Nu;            // symbol size (bit)
  int    Nb;            // = Ns*Nu (bit)
  int    Q;             // number of inner codewords |C|
  int    Dmax,Dmin;     // channel
  int    Drng;          // = Dmax-Dmin+1
  int    Nu2max,Nu2min; // recv symbol length
  long   Nu2p;          // = pow(2,Nu)
  double Pi,Pd,Ps,Pt;   // channel
  // k-mer依存遷移確率テーブル: [kmer][previous_event] -> [next_event_probabilities]
  std::map<std::string, std::map<int, std::vector<double>>> transition_probs;
  int    k_mer_length;  // k-merの長さ
  std::set<int> genew_loaded_ly;  // 既に読み込み済みのlyを記録
  double **GD;          // [Drng][Drng]: P(d_Nu | d_0)
  double ***GX;         // [Nu2][y][x]:  P(y|x,d0,d1)
  double ****GXNew;      // [Nu2][y][x][e]:  P(y|x,e) 新しい確率分布（イベント条件付き）
  double ****GENew;      // [Nu2][y][x][e]:  P(e_(i+1)v|e_iv, φ0(u'_(i-⌈k/v⌉+1)), φ0(u'_(i-⌈k/v⌉)), φ0(u'_i), φ0(u'_(i+1)), d_iv, d_(i+1)v)
  double **PF, **PB;    // [Ns+1][Drng]: FG forward/backward
  double **PU, **PD;    // [Ns][Q]:      FG up/down
  double **PI, **PO;    // [Ns][Q]:      FG input/output
  unsigned char *Yin;   // [Nb+Dmax]: received word
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
  void SetGXNew();
  void DelGXNew();
  double CalcPyxNew(long y, long x, int ly, int lx, int error_state);  // P(y|x,e)
  // 格子計算ヘルパー関数
  void initializeLattice(std::vector<std::vector<std::vector<double>>>& F, int lx, int ly);
  void calculateLatticeDP(std::vector<std::vector<std::vector<double>>>& F, 
                         const unsigned char* X, const unsigned char* Y, int lx, int ly);
  double getTransitionProbability(int prev_state, int curr_state);
  double GetGXNew(int Nu2, long y, long xi, int error_state);
  void SetGENew();
  void SetGENewMinimal();     // 最小限のテーブル初期化
  void ExpandGENew(int ly);   // 必要な ly のみ拡張
  void DelGENew();
  double CalcPexNew(long y, long x, int ly, int lx);  // P(e_(i+1)v|...) 後方互換性用
  double CalcPexNewWithErrorState(long y, long x, int ly, int lx, int target_error_state);  // エラー状態別確率計算
  std::array<double,4> CalcPexNewWithObservationIntegrated(long y, long x, int ly, int lx);  // calc_prob統合版
  
  // calc_prob統合システム用のヘルパー関数
  std::vector<int> longToBinaryVector(long value, int length);
  std::vector<char> longToDNAVector(long value, int length);
  int findPhi0IndexInCodebook(const std::vector<int>& pattern, int start, int length);
  std::string buildKmerFromPhi0PatternsIntegrated(int phi0_i_idx, int phi0_i1_idx, const std::vector<int>& m_pattern);
  double compute3DLatticeIntegrated(const std::string& kmer, const std::vector<char>& observed_y, int target_event);
  double getTransitionProbability(const std::string& kmer, int event);
  void outputLatticeDebug(const std::vector<std::vector<std::vector<double>>>& lattice,
                         int target_i, int target_j, const std::string& description,
                         std::ofstream& debug_file);
  double GetGENew(int Nu2, long y, long xi);
  // k-mer確率テーブル関連メソッド
  void loadTransitionProbabilities(int k);
  std::string binaryToDNA(const unsigned char* binary_seq, int length);
  std::string extractKmerFromSequence(const unsigned char* seq, int seq_len, int pos, int k);
  void SetFG();
  void DelFG();
  void InitFG(const unsigned char *RW, const double **Pin, int Nb2);
  void CalcPU(int idx);           // PU
  void CalcPO(int idx);           // PO  
  void CalcPF(int idx, int Nb2);  // PF
  void CalcPB(int idx, int Nb2);  // PB
  void CalcPD(int idx, int Nb2);  // PD
public:
  SLFBAdec(class InnerCodebook *_ICB, class ChannelMatrix *_ECM, class IDSchannel *_CH);
  ~SLFBAdec();
  void Decode(double **Pout, const unsigned char *RW, int Nb2, const int *dbgIW, const double **Pin);
  void Decode(double **Pout, const unsigned char *RW, int Nb2, const int *dbgIW);  // (uniform prior)
  void PrintNode(int idx, int iw);
  void PrintNode(const int *dbgIW);
  // 確率テーブル出力メソッド
  void exportTransitionProbabilities(const char* filename);
  void exportGXNewTable(const char* filename);
  void exportGENewTable(const char* filename);
  void exportAllProbabilityTables(const char* output_dir);
  // 格子計算デバッグ用メソッド
  void debugLatticeCalculation(long y, long x, int ly, int lx, const char* filename);
  // 事前計算システム用メソッド
  void PrecomputeGENewForLy(int ly);  // 事前計算専用（制限なし完全計算）
  void SaveGENewToFile(int ly, const char* filename);    // バイナリファイル保存
  bool LoadGENewFromFile(int ly, const char* filename);  // バイナリファイル読み込み
  bool LoadAllPrecomputedGENew(const char* precomputed_dir);  // 事前計算済みデータ一括読み込み
};
