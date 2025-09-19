#include <map>
#include <tuple>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
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
  double Pi,Pd,Ps,Pt;   // channel
  double **GD;          // [Drng][Drng]: P(d_Nu | d_0)
  double ***GX;         // [Nu2][y][x]:  P(y|x,d0,d1)
  double **PF, **PB;    // [Ns+1][Drng]: FG forward/backward
  double **PU, **PD;    // [Ns][Q]:      FG up/down
  double **PI, **PO;    // [Ns][Q]:      FG input/output
  unsigned char *Yin;   // [Nb+Dmax]: received word
  // エラー状態関連データ構造 (3D)
  double ***PFE, ***PBE; // [Ns+1][Drng][NUM_ERROR_STATES]: エラー状態を含む前進/後進確率
  double ***PE;          // [Ns][NUM_ERROR_STATES][NUM_ERROR_STATES]: エラー状態遷移確率
  
  // k-mer依存ラティス関連データ構造（分離アプローチ）
  static const int KMER_LENGTH = 4;           // k-mer長 (k=4固定)
  int num_kmers;                              // 4^KMER_LENGTH = 256 for k=4

  // 【新アプローチ】横方向メッセージ（3次元）- 時間軸での情報伝播用
  double ***FwdMsg;                           // [Ns+1][Drng][NUM_ERROR_STATES]: 前向きメッセージ P(d_t, e_t | y_1..t)
  double ***BwdMsg;                           // [Ns+1][Drng][NUM_ERROR_STATES]: 後ろ向きメッセージ P(y_t+1..T | d_t, e_t)

  // 【互換性用】既存4次元配列（段階的移行用、最終的に削除予定）
  double ****PFE4D, ****PBE4D;               // [Ns+1][Drng][NUM_ERROR_STATES][num_kmers]: 4D前進/後進確率
  
  // 論文式(41): P(e_{t+1} | e_t, η_{t+1}) のテーブル
  // キー: (直前のエラー e_t, 次のk-mer η_{t+1})
  // 値: [P(e_{t+1}=Match), P(e_{t+1}=Insertion), P(e_{t+1}=Deletion), P(e_{t+1}=Substitution)]
  void* kmer_error_probs;                     // std::map<std::pair<int,int>, std::vector<double>>* として使用

  // メモ化キャッシュ：遷移確率のオンデマンド計算（メモリ効率）
  // キー: (e0, k0, xi, xi_next, e1) → 値: P(e1|e0, k0, xi, xi_next)
  std::map<std::tuple<int, int, int, int, int>, double> transition_prob_cache;
  
  // 4値置換確率行列 (論文Table 2対応)
  double SubMatrix[4][4];                     // P(b|a) = aがbに置換される条件付き確率 (A=0, C=1, G=2, T=3)
  class InnerCodebook *ICB;
  class ChannelMatrix *ECM;   // (encoding CM)
  class IDSchannel    *CH;
  //-----
  int    max(int a, int b);
  int    min(int a, int b);
  double Psub(unsigned char a, unsigned char b);              // substitution prob (2値用・レガシー)
  double Psub_quaternary(unsigned char a, unsigned char b);  // substitution prob (4値用・論文Table 2対応)
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
  
  // エラー状態関連関数 (3D)
  void SetFGE();                  // エラー状態拡張FGデータ構造の設定
  void DelFGE();                  // エラー状態拡張FGデータ構造の削除
  void InitFGE(const unsigned char *RW, const double **Pin, int Nb2); // エラー状態拡張FGの初期化
  void CalcPE(int idx);           // エラー状態遷移確率の計算
  void CalcPFE(int idx, int Nb2); // エラー状態を含む前進確率
  void CalcPBE(int idx, int Nb2); // エラー状態を含む後進確率
  void CalcPDE(int idx, int Nb2); // エラー状態を含む事後確率計算
  
  // k-mer依存4次元ラティス関数
  void SetFGE4D();                     // 4次元データ構造の設定
  void DelFGE4D();                     // 4次元データ構造の削除
  void InitFGE4D(const unsigned char *RW, const double **Pin, int Nb2); // 4次元FGの初期化
  int  GetKmerIndex(const unsigned char *kmer_seq);      // k-mer配列をインデックスに変換
  void GetKmerFromIndex(int kmer_idx, unsigned char *kmer_seq); // インデックスからk-mer配列に変換
  int  ComputeNextKmer(int current_kmer, int codeword_xi);     // 決定論的k-mer遷移: k_t + xi → k_{t+1}
  double GetKmerErrorProb(int prev_error, int next_kmer, int next_error); // 論文式(41)の確率を取得
  double ComputeObservationProbabilityFromDNA(int idx, int Nu2, int iL, int xi, int k0, int e0, int Nb2); // DNAチャネル統合観測確率計算
  double CalcPyx_dynamic(long y, long x, int ly, int lx, int prev_error, int kmer, int codeword_xi); // 動的確率を使った観測確率計算(CalcPyxの進化版)
  double Psub_quaternary(unsigned char a, unsigned char b, double dynamic_ps); // 動的ps対応のPsub_quaternaryオーバーロード
  // 【新アプローチ】分離型計算メソッド
  void CalcForwardMsg(int idx, int Nb2);     // 横方向前向きメッセージ計算（k次元を周辺化）
  void CalcBackwardMsg(int idx, int Nb2);    // 横方向後ろ向きメッセージ計算（k次元を周辺化）
  void CalcPosterior(int idx, int Nb2);      // 縦方向事後確率計算（各時刻での内部計算）

  // 【ビームサーチ版】計算量削減メソッド
  void CalcForwardMsgBeam(int idx, int Nb2, double sparse_threshold); // ビームサーチ版前向きメッセージ計算

  

  // ビームサーチパラメータ
  // int beam_width;  
  size_t beam_width; 

  // 【動的枝刈りメモ化】計算量削減用データ構造
  // 1つの「高確率な遷移先」の情報
  struct SparseTransitionInfo {
    int next_e;     // 到達点のエラー状態 e₁
    int next_k;     // 到達点のk-mer状態 k₁
    double prob;    // 遷移確率 P(e₁ | e₀, k₀, xi, xi_next)
  };

  // 動的キャッシュ: 遷移先リストをオンデマンドで計算・保存
  // キー: tuple(e₀, k₀, xi, xi_next) -> 値: 高確率遷移先リスト
  using SparseTransitionCache = std::map<std::tuple<int, int, int, int>, std::vector<SparseTransitionInfo>>;
  SparseTransitionCache sparse_transition_cache;

  // 【互換性用】既存4次元計算メソッド（段階的移行用）
  void CalcPFE4D(int idx, int Nb2);    // 4次元前進確率計算
  void CalcPBE4D(int idx, int Nb2);    // 4次元後進確率計算
  void CalcPDE4D(int idx, int Nb2);    // 4次元事後確率計算
  // 拡張k-mer関連の新しい関数（複数符号語対応）
  int ComputeExtendedKmer(int current_kmer, int xi_current, int xi_next, int idx);
  double GetExtendedKmerErrorProb(int prev_error, int prev_kmer, int xi_current, int xi_next, int next_error);
  double ComputeExtendedObservationProbability(int idx, int Nu2, int iL, int xi_current, int xi_next, int k0, int e0, int Nb2);
  
  // 複数符号語系列を考慮した観測確率計算（現在未使用）
  // double CalcPyx_dynamic_extended(long y, long x, int ly, int lx, int prev_error, int kmer, int xi_current, int xi_next);
  

  // ✅ CalcPxy_dynamicの計算結果をキャッシュするマップ
  // キー: <y, x, ly, lx, prev_e, kmer, xi> の組み合わせ
  // 値: 計算結果のdouble値
  std::map<std::tuple<long, long, int, int, int, int, int>, double> pyx_cache;

  // 【監視システム】Mutex と監視スレッド制御
  mutable std::mutex cache_mutex;          // transition_prob_cache への排他制御
  std::atomic<bool> monitoring_active;     // 監視スレッドの動作制御
  std::thread monitoring_thread;           // 監視スレッドのハンドル
public:
  SLFBAdec(class InnerCodebook *_ICB, class ChannelMatrix *_ECM, class IDSchannel *_CH);
  ~SLFBAdec();
  void Decode(double **Pout, const unsigned char *RW, int Nb2, const int *dbgIW, const double **Pin, double sparse_threshold = 1e-6);
  void Decode(double **Pout, const unsigned char *RW, int Nb2, const int *dbgIW, double sparse_threshold = 1e-6);  // (uniform prior)
  void PrintNode(int idx, int iw);
  void PrintNode(const int *dbgIW);
  // テーブル出力機能
  void exportGDTable(const char* filename);
  void exportGXTable(const char* filename);
  void exportAllTables(const char* output_dir);

  // メモ化遷移確率関数（メモリ効率版） - removed, using GetValidTransitions instead
  void ClearTransitionProbCache();          // キャッシュクリア（メモリ解放）
  void ExportTransitionProbCache(const char* filename) const;          // キャッシュ内容をファイル出力（デバッグ用）
  void LoadKmerErrorProbabilities(const char* dir_path); // DNArSim-mainから P(e_{t+1}|e_t,η_{t+1}) を読み込み（公開API）
  void ExportSparseTransitions(const char* filename) const;

  // 【キャッシュ機能】スパース遷移行列の保存・読み込み
  void SaveSparseTransitions(const char* filename) const;
  bool LoadSparseTransitions(const char* filename);

  // 【監視システム】リアルタイムキャッシュ監視機能
  void StartCacheMonitoring(const char* snapshot_filename = "cache_snapshot.txt", int interval_seconds = 30);
  void StopCacheMonitoring();

  // 【動的枝刈りメモ化】新しい核心関数
  const std::vector<SparseTransitionInfo>& GetValidTransitions(int e0, int k0, int xi, int xi_next, double threshold = 1e-6);

 
};
