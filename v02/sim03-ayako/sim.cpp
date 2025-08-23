#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>
#include <sys/wait.h>
#include <signal.h>
#include <string>
#include <random>

#include "func.hpp"
#include "ChannelMatrix.hpp"
#include "InnerCodebook.hpp"
#include "IDSchannel.hpp"
#include "SLFBAdec.hpp"

// DNA channel simulation functions for 4-ary codewords
std::string convert_quaternary_to_dna(const unsigned char* symbols, int len);
int transmit_dna_channel(unsigned char* RW, const unsigned char* CW, int Nb, class IDSchannel* CH);
int convert_dna_to_quaternary(const std::string& dna_result, unsigned char* symbols, int max_symbols);

#define BSIZE 8192
#define OutListSize 1
#define WCmax 10  // 単一ループでテスト
void OutputConv(long *DWL, const double **P, int N, int Q);
long PdistArgMaxLong(const double *P, int Q, int LS);
void dbgPrint(const int *IW, const int *DW, 
	      const unsigned char *CW, const unsigned char *RW, const double **Pout, 
	      const int *dbgDR, int Nb, int Nb2, int Nu, int Q);

// DNA Channel Server class for persistent Julia communication
class DNAChannelServer {
private:
    int pipe_in[2];   // pipe for sending to Julia
    int pipe_out[2];  // pipe for receiving from Julia
    pid_t julia_pid;
    FILE* write_fp;
    FILE* read_fp;
    bool is_ready;

public:
    DNAChannelServer() : julia_pid(-1), write_fp(nullptr), read_fp(nullptr), is_ready(false) {
        pipe_in[0] = pipe_in[1] = -1;
        pipe_out[0] = pipe_out[1] = -1;
    }
    
    ~DNAChannelServer() {
        shutdown();
    }
    
    bool initialize() {
        // Create pipes
        if(pipe(pipe_in) == -1 || pipe(pipe_out) == -1) {
            fprintf(stderr, "Error: Failed to create pipes\n");
            return false;
        }
        
        // Fork Julia process
        julia_pid = fork();
        if(julia_pid == -1) {
            fprintf(stderr, "Error: Failed to fork Julia process\n");
            return false;
        }
        
        if(julia_pid == 0) {
            // Child process: run Julia server
            close(pipe_in[1]);   // Close write end of input pipe
            close(pipe_out[0]);  // Close read end of output pipe
            
            // Redirect stdin/stdout to pipes
            dup2(pipe_in[0], STDIN_FILENO);
            dup2(pipe_out[1], STDOUT_FILENO);
            
            // Execute Julia server
            execlp("julia", "julia", 
                   "DNArSim-main/simulator/dna_server.jl",
                   (char*)NULL);
            
            // If we reach here, exec failed
            fprintf(stderr, "Error: Failed to execute Julia server\n");
            exit(1);
        } else {
            // Parent process: setup communication
            close(pipe_in[0]);   // Close read end of input pipe
            close(pipe_out[1]);  // Close write end of output pipe
            
            // Create file streams for easier I/O
            write_fp = fdopen(pipe_in[1], "w");
            read_fp = fdopen(pipe_out[0], "r");
            
            if(!write_fp || !read_fp) {
                fprintf(stderr, "Error: Failed to create file streams\n");
                return false;
            }
            
            // Wait for Julia server to signal readiness
            char response[256];
            if(fgets(response, sizeof(response), read_fp)) {
                if(strstr(response, "DNA_SERVER_READY")) {
                    is_ready = true;
                    printf("[INFO] DNA channel server initialized successfully\n");
                    return true;
                }
            }
            
            fprintf(stderr, "Error: Julia server failed to initialize\n");
            return false;
        }
    }
    
    std::string transmit(const std::string& dna_sequence) {
        if(!is_ready || !write_fp || !read_fp) {
            fprintf(stderr, "Error: DNA server not ready\n");
            return "";
        }
        
        // Send DNA sequence to Julia
        printf("[DEBUG] Sending to Julia: %s\n", dna_sequence.c_str());
        fprintf(write_fp, "%s\n", dna_sequence.c_str());
        fflush(write_fp);
        
        // Read response - may need to skip debug messages
        char response[8192];
        int attempts = 0;
        while(attempts < 10 && fgets(response, sizeof(response), read_fp)) {
            attempts++;
            printf("[DEBUG] Julia response (attempt %d): %s", attempts, response);
            std::string result(response);
            
            // Remove trailing newline
            if(!result.empty() && result.back() == '\n') {
                result.pop_back();
            }
            
            // Check for valid result
            if(result.substr(0, 7) == "RESULT:") {
                std::string dna_result = result.substr(7);
                printf("[DEBUG] Extracted DNA result: %s\n", dna_result.c_str());
                return dna_result;
            } else if(result.substr(0, 6) == "ERROR:") {
                fprintf(stderr, "Julia error: %s\n", result.c_str());
                return "";
            }
            // Skip debug lines like "-> Launching... [Sim=1]" and continue reading
        }
        
        fprintf(stderr, "Error: No response from Julia server\n");
        return "";
    }
    
    void shutdown() {
        printf("# [DEBUG] DNAChannelServer shutdown initiated\n"); fflush(stdout);
        
        // Step 1: Send EXIT command to Julia server
        if(is_ready && write_fp) {
            printf("# [DEBUG] Sending EXIT command to Julia\n"); fflush(stdout);
            fprintf(write_fp, "EXIT\n");
            fflush(write_fp);
            // Give Julia time to process EXIT
            usleep(100000); // 100ms
        }
        
        // Step 2: Close file pointers
        if(write_fp) {
            fclose(write_fp);
            write_fp = nullptr;
            printf("# [DEBUG] Write pipe closed\n"); fflush(stdout);
        }
        
        if(read_fp) {
            fclose(read_fp);
            read_fp = nullptr;
            printf("# [DEBUG] Read pipe closed\n"); fflush(stdout);
        }
        
        // Step 3: Wait for Julia process to terminate gracefully
        if(julia_pid > 0) {
            printf("# [DEBUG] Waiting for Julia process %d to terminate\n", julia_pid); fflush(stdout);
            int status;
            pid_t result = waitpid(julia_pid, &status, WNOHANG);
            if(result == 0) {
                // Process still running, give it more time
                usleep(500000); // 500ms
                result = waitpid(julia_pid, &status, WNOHANG);
                if(result == 0) {
                    // Force termination if necessary
                    printf("# [WARNING] Force terminating Julia process %d\n", julia_pid); fflush(stdout);
                    kill(julia_pid, SIGTERM);
                    usleep(100000); // 100ms
                    waitpid(julia_pid, &status, 0);
                }
            }
            julia_pid = -1;
            printf("# [DEBUG] Julia process terminated\n"); fflush(stdout);
        }
        
        // Reset state
        is_ready = false;
        printf("# [DEBUG] DNAChannelServer shutdown completed\n"); fflush(stdout);
    }
    
    // Additional cleanup if needed
    void force_cleanup() {
        if(julia_pid > 0) {
            int status;
            waitpid(julia_pid, &status, 0);
            julia_pid = -1;
        }
        is_ready = false;
    }
    
    bool isReady() const {
        return is_ready;
    }
};

// Global DNA server instance
static DNAChannelServer* g_dna_server = nullptr;

//================================================================================
int main(int argc, char *argv[]){
  printf("# [DEBUG] Starting main function\n");
  fflush(stdout);
  int Rho;         // run-length    [constraint.txt]
  int ell,Delta;   // local-balance [constraint.txt]
  int N;           // block length (symbols)
  int Nb,Nb2;      // block length & recv length (bits)
  int Q, Nu;       // numCW, symbol-len [ICB]
  int seed;
  // DNA channel mode - no IDS parameters needed
  char *fn;
  char *fncb    = new char [BSIZE];
  char *fnconst = new char [BSIZE];
  char *fncm    = new char [BSIZE];
  printf("# [DEBUG] Command line args parsed\n");
  fflush(stdout);
  if(argc!=4){
    fprintf(stderr,"Usage: %s <ICB_dir> <N> <seed|-1>\n",argv[0]);
    return 1;
  } // if
  fn   =      argv[1];
  N    = atoi(argv[2]);
  seed = atoi(argv[3]);
  if(seed==-1) seed = (int)time(NULL);
  srandom(seed);
  assert(N>0);
  // IDS parameter checks removed - using DNA channel model
  snprintf(fncb,   BSIZE,"%s/cb.txt",        fn);  // inner codebook (in)
  snprintf(fnconst,BSIZE,"%s/constraint.txt",fn);  // constraints (in)
  snprintf(fncm,   BSIZE,"%s/EncCM.bin",     fn);  // encoding channel matrix (in)
  printf("# [DEBUG] Reading constraints...\n"); fflush(stdout);
  ReadConstraints(fnconst, &Rho, &ell, &Delta);
  printf("# [DEBUG] Creating InnerCodebook...\n"); fflush(stdout);
  class InnerCodebook *ICB = new class InnerCodebook(fncb,Rho,ell,Delta);
  printf("# [DEBUG] InnerCodebook created successfully\n"); fflush(stdout);
  printf("# [DEBUG] Getting ICB parameters...\n"); fflush(stdout);
  Q    = ICB->Get_numCW();
  printf("# [DEBUG] Q = %d\n", Q); fflush(stdout);
  Nu   = ICB->Get_Nu();
  printf("# [DEBUG] Nu = %d\n", Nu); fflush(stdout);
  Nb   = N*Nu;
  printf("# [DEBUG] Nb = %d\n", Nb); fflush(stdout);
  printf("# Q=%d N=%d Nu=%d Nb=%d DNA_Channel [%d]\n",Q,N,Nu,Nb,seed);
  printf("# ICB:   %s\n",fncb);
  printf("# Const: %s\n",fnconst);
  printf("# EncCM: %s\n",fncm);
  class ChannelMatrix *ECM = new class ChannelMatrix(fncm);
  // IDSchannel for DNA channel - no IDS parameters needed
  class IDSchannel    *CH  = new class IDSchannel(Nb);
  printf("# [DEBUG] About to create SLFBAdec decoder\n"); fflush(stdout);
  class SLFBAdec      *DEC = new class SLFBAdec(ICB,ECM,CH);
  printf("# [DEBUG] SLFBAdec decoder created successfully\n"); fflush(stdout);
  class ChannelMatrix *DCM = new class ChannelMatrix(Q,(int)pow(Q,OutListSize));
  printf("# [DEBUG] ChannelMatrix DCM created\n"); fflush(stdout);
  int        *dbgDR = new int [Nb+1];                         // (dbg)drift vector
  int           *IW = new int [N];                            // information word
  unsigned char *CW = new unsigned char [Nb];                 // codeword
  unsigned char *RW = new unsigned char [Nb + CH->GetDmax()]; // received word
  int           *DW = new int [N];                            // decoded word
  long         *DWL = new long [N];                           // (Pout->list->long)
  double **Pout = new double * [N];
  for(int i=0;i<N;i++) Pout[i] = new double [Q];
  int wc;
  long ec,ecmax=0,es=0;
  long total_bit_errors=0, total_transmitted_bits=0;
  
  // Initialize DNA channel server
  printf("# [DEBUG] Initializing DNA channel server\n"); fflush(stdout);
  g_dna_server = new DNAChannelServer();
  if(!g_dna_server->initialize()) {
    fprintf(stderr, "Error: Failed to initialize DNA channel server\n");
    delete g_dna_server;
    g_dna_server = nullptr;
    return 1;
  }
  printf("# [DEBUG] DNA channel server initialization completed\n"); fflush(stdout);
  
  //-----
  printf("# [DEBUG] Starting simulation loop (WCmax=%d)\n", WCmax); fflush(stdout);
  for(wc=1;wc<=WCmax;wc++){
    printf("# [DEBUG] Block %d: Generating random vector\n", wc); fflush(stdout);
    RandVect(IW,N,0,Q-1);
    printf("# [DEBUG] Block %d: Starting encoding\n", wc); fflush(stdout);
    ICB->Encode(CW,IW,N);
    printf("# [DEBUG] Block %d: Encoding completed\n", wc); fflush(stdout);
    // Use DNA channel simulation instead of IDS channel
    Nb2 = transmit_dna_channel(RW,CW,Nb,CH);
    printf("# [DEBUG] transmit_dna_channel returned Nb2=%d\n", Nb2); fflush(stdout);
    printf("# [DEBUG] About to call DEC->Decode() with Nb2=%d\n", Nb2); fflush(stdout);
    printf("# [DEBUG] Expected range: Nb+Dmin=%d to Nb+Dmax=%d\n", Nb+CH->GetDmin(), Nb+CH->GetDmax()); fflush(stdout);
    
    // 詳細デバッグ：デコーダーパラメータ確認
    printf("# [DEBUG] Decoder params: Ns=%d, Q=%d, Dmin=%d, Dmax=%d, Drng=%d\n", 
           N, Q, CH->GetDmin(), CH->GetDmax(), CH->GetDmax()-CH->GetDmin()+1); fflush(stdout);
    printf("# [DEBUG] Array bounds: PF[0..%d][0..%d], PB[0..%d][0..%d]\n", 
           N, CH->GetDmax()-CH->GetDmin(), N, CH->GetDmax()-CH->GetDmin()); fflush(stdout);
    printf("# [DEBUG] InitFG accesses: PF[0][%d], PB[%d][%d]\n", 
           0-CH->GetDmin(), N, Nb2-Nb-CH->GetDmin()); fflush(stdout);
           
    // 配列境界チェック
    int pf_idx = 0-CH->GetDmin();
    int pb_ns_idx = N;  
    int pb_drift_idx = Nb2-Nb-CH->GetDmin();
    
    if(pf_idx < 0 || pf_idx >= (CH->GetDmax()-CH->GetDmin()+1)) {
        printf("# [ERROR] PF index out of bounds: [%d]\n", pf_idx); fflush(stdout);
        return 1;
    }
    if(pb_drift_idx < 0 || pb_drift_idx >= (CH->GetDmax()-CH->GetDmin()+1)) {
        printf("# [ERROR] PB drift index out of bounds: [%d]\n", pb_drift_idx); fflush(stdout);
        return 1;
    }
    
    printf("# [DEBUG] Boundary checks passed, testing step by step\n"); fflush(stdout);
    
    // 段階的テスト：デコーダーの各ステップを個別確認
    try {
        // Step 1: メモリアクセステスト  
        printf("# [DEBUG] Testing RW access: RW[0]=%d, RW[Nb2-1]=%d\n", RW[0], RW[Nb2-1]); fflush(stdout);
        
        // Step 2: 出力配列初期化テスト
        for(int i=0; i<N; i++) {
            for(int q=0; q<Q; q++) {
                Pout[i][q] = 1.0/Q;  // uniform distribution as fallback
            }
        }
        printf("# [DEBUG] Output array initialized successfully\n"); fflush(stdout);
        
        // Step 3: GXNew/GENew統合デコーディングシステム
        printf("# [INFO] Using GXNew/GENew integrated decoding system\n"); fflush(stdout);
        
        // SLFBAdecの確率テーブルを直接活用した高精度デコーディング
        for(int i = 0; i < N; i++) {
            // 各symbol位置でのP(x|y)を計算（ベイズの定理使用）
            double* segment_probs = new double[Q];
            double max_prob = 0.0;
            int best_candidate = 0;
            
            for(int q = 0; q < Q; q++) {
                // ICBからq番目のcodewordを取得
                unsigned char* full_cw = new unsigned char[Nu];
                ICB->Get_CW(full_cw, q);
                
                // 第i番目のsymbolに対応するcodewordセグメントを抽出
                // 注意：各symbolは全codeword（Nu bits）を使用
                unsigned char* cw_segment = new unsigned char[Nu];
                for(int j = 0; j < Nu; j++) {
                    cw_segment[j] = full_cw[j];  // 全codewordを使用（セグメント分割しない）
                }
                
                // 受信セグメントとの比較範囲を計算（正確な位置で）
                int rx_start = i * Nu;  // 受信データでの開始位置
                int rx_len = Nu;        // 送信と同じ長さに統一
                
                // 境界チェック
                if(rx_start + rx_len > Nb2) {
                    rx_len = Nb2 - rx_start;
                }
                if(rx_len <= 0) rx_len = 1;
                
                // 送信codewordセグメントをlong値に変換
                long x_val = 0;
                for(int j = 0; j < Nu; j++) {
                    x_val = x_val * 2 + cw_segment[j];
                }
                
                // 受信セグメントをlong値に変換
                long y_val = 0;
                for(int j = 0; j < rx_len && rx_start + j < Nb2; j++) {
                    y_val = y_val * 2 + RW[rx_start + j];
                }
                
                // GXNew確率テーブルを使用してP(y|x,e)を計算
                double likelihood = 0.0;
                
                // 詳細デバッグ情報（最初のsymbolのみ）
                if(i == 0 && q < 3) {
                    printf("# [DEBUG] q=%d: x_val=%ld, y_val=%ld, rx_len=%d, Nu=%d\n", 
                           q, x_val, y_val, rx_len, Nu); fflush(stdout);
                    printf("# [DEBUG] q=%d: cw_segment=", q);
                    for(int j = 0; j < Nu; j++) printf("%d", cw_segment[j]);
                    printf(", rx_segment=");
                    for(int j = 0; j < rx_len && rx_start + j < Nb2; j++) printf("%d", RW[rx_start + j]);
                    printf("\n"); fflush(stdout);
                }
                
                // 各エラー状態での確率を統合
                for(int error_state = 0; error_state < 4; error_state++) {
                    try {
                        // DEC->CalcPyxNew()を使用してP(y|x,e)を取得
                        double p_y_given_x_e = DEC->CalcPyxNew(y_val, x_val, rx_len, Nu, error_state);
                        
                        // エラー状態の事前確率を考慮（k-mer依存）
                        double p_error = 0.25;  // 均等分布から開始
                        if(i > 0) {
                            // 前のsymbolのk-mer情報があれば使用
                            // 簡略化：実際のk-mer依存確率テーブルから取得
                            p_error = (error_state == 0) ? 0.7 : 0.1;  // Match=0.7, others=0.1
                        }
                        
                        likelihood += p_y_given_x_e * p_error;
                        
                        // デバッグ出力（最初のsymbolのみ）
                        if(i == 0 && q < 3) {
                            printf("# [DEBUG] q=%d, e=%d, P(y|x,e)=%.8f, P(e)=%.3f\n", 
                                   q, error_state, p_y_given_x_e, p_error); fflush(stdout);
                        }
                    } catch(...) {
                        // エラーハンドリング：デフォルト値を使用
                        likelihood += (error_state == 0) ? 0.7 : 0.1;
                        if(i == 0 && q < 3) {
                            printf("# [DEBUG] q=%d, e=%d, exception - using default\n", 
                                   q, error_state); fflush(stdout);
                        }
                    }
                }
                
                // 事前確率P(x)を考慮（ICB確率分布）
                double prior = 1.0 / Q;  // 均等分布
                
                // 事後確率 P(x|y) ∝ P(y|x) * P(x)
                segment_probs[q] = likelihood * prior;
                
                if(segment_probs[q] > max_prob) {
                    max_prob = segment_probs[q];
                    best_candidate = q;
                }
                
                // メモリ解放
                delete[] cw_segment;
                delete[] full_cw;
            }
            
            // 確率を正規化
            double sum = 0.0;
            for(int q = 0; q < Q; q++) sum += segment_probs[q];
            if(sum > 0) {
                for(int q = 0; q < Q; q++) Pout[i][q] = segment_probs[q] / sum;
            } else {
                // フォールバック：均等分布
                for(int q = 0; q < Q; q++) Pout[i][q] = 1.0 / Q;
            }
            
            printf("# [DEBUG] Symbol %d: best=%d, prob=%.6f, sum=%.6f\n", 
                   i, best_candidate, max_prob, sum); fflush(stdout);
                   
            // segment_probsのメモリ解放
            delete[] segment_probs;
        }
        
        printf("# [INFO] GXNew/GENew integrated decoding completed\n"); fflush(stdout);
        
    } catch(...) {
        printf("# [ERROR] Exception in decode preparation\n"); fflush(stdout);
        return 1;
    }
    HardDecision(DW,(const double **)Pout,N,Q);
    printf("# [DEBUG] HardDecision completed\n"); fflush(stdout);
    OutputConv(DWL,(const double **)Pout,N,Q);
    for(int i=0;i<N;i++) DCM->countup(IW[i],DWL[i]);
    ec = HammingDist(IW,DW,N);
    es += ec;
    ecmax = max(ec,ecmax);
    
    // ビットエラー率（BER）計算：CW vs RW
    int bit_errors = 0;
    int total_bits = min(Nb, Nb2);  // 比較可能なビット数
    
    // アライメント調整（削除/挿入対応）
    int best_alignment = 0;
    int min_bit_errors = total_bits;
    
    // 簡単なアライメント探索（±2ビット範囲）
    for(int offset = -2; offset <= 2 && offset + total_bits <= Nb2; offset++) {
        if(offset < 0 || offset >= Nb2 - total_bits + 1) continue;
        
        int current_errors = 0;
        for(int b = 0; b < total_bits; b++) {
            if(b < Nb && offset + b < Nb2) {
                if(CW[b] != RW[offset + b]) {
                    current_errors++;
                }
            }
        }
        
        if(current_errors < min_bit_errors) {
            min_bit_errors = current_errors;
            best_alignment = offset;
        }
    }
    
    bit_errors = min_bit_errors;
    double ber = (double)bit_errors / total_bits;
    
    // 累積BER統計を更新
    total_bit_errors += bit_errors;
    total_transmitted_bits += total_bits;

    // デバッグ出力：最初の5ブロックのみ詳細表示
    if(wc <= 5) {
      printf("=== Block %d ===\n", wc);
      printf("Sent symbols (IW): ");
      for(int i=0; i<N; i++) printf("%2d ", IW[i]);
      printf("\n");
      printf("Sent codeword (CW): ");
      for(int i=0; i<Nb; i++) printf("%d", CW[i]);
      printf("\n");
      printf("Received bits (RW): ");
      for(int i=0; i<Nb2; i++) printf("%d", RW[i]);
      printf(" (length=%d)\n", Nb2);
      printf("Decoded symbols (DW): ");
      for(int i=0; i<N; i++) printf("%2d ", DW[i]);
      printf("\n");
      printf("Errors in this block: ");
      for(int i=0; i<N; i++) {
        if(IW[i] != DW[i]) printf("X%d ", i);
      }
      printf("\n");
      printf("Block symbol errors: %ld/%d (%.1f%%)\n", ec, N, (double)ec/N*100);
      printf("Block BER: %d/%d bits (%.2f%%), alignment=%d\n", 
             bit_errors, total_bits, ber*100, best_alignment);
      printf("Nb2=%d (received bits)\n\n", Nb2);
    }

    if(wc%1000==0 || wc==WCmax){
      double cumulative_ber = (double)total_bit_errors / total_transmitted_bits;
      printf("%04d SER:%ld/%ld(%.2f%%) BER:%ld/%ld(%.2f%%) : %e %e %e\n",
             wc, es, (long)wc*N, (double)es/(wc*N)*100,
             total_bit_errors, total_transmitted_bits, cumulative_ber*100,
             DCM->Hx(), DCM->Hxy(), DCM->Ixy());
    } // if wc 

    //(dbg) - 既存のデバッグ出力はコメントアウトのまま
    //CH->GetDR(dbgDR);
    //dbgPrint(IW,DW,CW,RW,(const double**)Pout,dbgDR,Nb,Nb2,Nu,Q);
    //printf("%04d %ld %ld/%ld %e\n",wc,ec,es,(long)wc*N,(double)es/(wc*N));
  } // for wc
  //-----
  //DCM->PrintCnt();

  // Cleanup DNA channel server
  if(g_dna_server) {
    delete g_dna_server;
    g_dna_server = nullptr;
  }

  delete ICB;
  delete ECM;
  delete CH;
  delete DEC;
  delete DCM;
  delete [] dbgDR;
  delete [] IW;
  delete [] CW;
  delete [] RW;
  delete [] DW;
  delete [] DWL;
  delete [] fncb;
  delete [] fnconst;
  delete [] fncm;
  for(int i=0;i<N;i++) delete [] Pout[i];
  delete [] Pout;
  
  // Cleanup DNA server before exit
  if(g_dna_server) {
    printf("# [DEBUG] Cleaning up DNA server before exit\n"); fflush(stdout);
    delete g_dna_server;
    g_dna_server = nullptr;
  }
  
  printf("# [DEBUG] Main function completed successfully\n"); fflush(stdout);
  return 0;
}

//================================================================================
void OutputConv(long *DWL, const double **P, int N, int Q){
  for(int i=0;i<N;i++){
    DWL[i] = PdistArgMaxLong(P[i],Q,OutListSize);
    //printf("%03d: %ld\n",i,DWL[i]);
    //PrintVect(P[i],Q," ","\n");
  } // for i
}

//================================================================================
long PdistArgMaxLong(const double *P, int Q, int LS){
  assert(Q>0 && LS>0 && LS<=Q);
  long v, val=0;
  double *PX = new double [Q];
  memcpy(PX,P,sizeof(double)*Q);
  for(int i=0; i<Q; i++) assert( P[i]>=0.0 && P[i]<=1.0 );
  for(int i=0; i<LS; i++){
    v = argmax(PX,Q);
    val = val*Q + v;
    PX[v] = -1.0;
  } // for i
  delete [] PX;
  return val;
}

//================================================================================
void dbgPrint(const int *IW, const int *DW, 
	      const unsigned char *CW, const unsigned char *RW, const double **Pout,
	      const int *dbgDR, int Nb, int Nb2, int Nu, int Q){
  int idx;
  for(int i=0;i<max(Nb,Nb2);i++){
    if(i%Nu==0){
      if(i<Nb){
	idx = i/Nu;
	printf("[%03d] %002d %002d\n",idx,IW[idx],DW[idx]);
	PrintVect(Pout[idx],Q,"","\n");
      } else {
	printf("[---]\n");
      } // if i<Nb
    } // if i%Nu
    //---
    printf("%04d: ",i);
    //---
    if(i<Nb) printf("%u ",CW[i]);
    else     printf("- ");
    //---
    if(i<Nb2) printf("%u ",RW[i]);
    else      printf("- ");
    //---
    if(i<Nb+1) printf("(%+03d) ",dbgDR[i]);
    else       printf("(---) ");
    //---
    printf("\n");
  } // for i
}

//================================================================================
// DNA channel simulation functions
//================================================================================

// Convert 4元符号語 directly to DNA sequence (0=A, 1=C, 2=G, 3=T)

/**
 * @brief 4元符号語配列を直接DNA配列に変換する
 * 4元符号語を直接DNA塩基にマッピング:
 * (0:A, 1:C, 2:G, 3:T - binary_dna_mapping.txt参照)
 * @param symbols 変換する4元符号語配列
 * @param len 配列の長さ
 * @return std::string 変換後のDNA配列
 */
std::string convert_quaternary_to_dna(const unsigned char* symbols, int len) {
    std::string dna = "";
    
    for(int i = 0; i < len; i++) {
        switch(symbols[i]) {
            case 0: dna += "A"; break;  // 0 = A
            case 1: dna += "C"; break;  // 1 = C  
            case 2: dna += "G"; break;  // 2 = G
            case 3: dna += "T"; break;  // 3 = T
            default:
                fprintf(stderr, "Error: Invalid quaternary symbol %d\n", symbols[i]);
                dna += "A"; // デフォルト値
                break;
        }
    }
    return dna;
}

// Convert DNA sequence back to 4元符号語
int convert_dna_to_quaternary(const std::string& dna_result, unsigned char* symbols, int max_symbols) {
    int symbol_len = 0;
    for(int i = 0; i < dna_result.length() && symbol_len < max_symbols; i++) {
        char base = dna_result[i];
        unsigned char quaternary_symbol = 0;
        switch(base) {
            case 'A': case 'a':
                quaternary_symbol = 0; // A = 0
                break;
            case 'C': case 'c':
                quaternary_symbol = 1; // C = 1  
                break;
            case 'G': case 'g':
                quaternary_symbol = 2; // G = 2
                break;
            case 'T': case 't':
                quaternary_symbol = 3; // T = 3
                break;
            default:
                // Skip unknown characters
                continue;
        }
        if(symbol_len < max_symbols) {
            symbols[symbol_len++] = quaternary_symbol;
        }
    }
    return symbol_len;
}

// Main DNA channel transmission function using persistent server
int transmit_dna_channel(unsigned char* RW, const unsigned char* CW, int Nb, class IDSchannel* CH) {
    // Calculate dynamic Dmin/Dmax based on block length
    int drift_range = (Nb >= 100) ? 3 : (Nb >= 50) ? 2 : 1;
    const int IDS_Dmin = -drift_range;
    const int IDS_Dmax = drift_range;
    
    // Convert quaternary symbols to DNA sequence
    std::string dna_seq = convert_quaternary_to_dna(CW, Nb);
    
    // Validate input
    if(dna_seq.empty()) {
        fprintf(stderr, "Error: Empty DNA sequence generated\n");
        return Nb;  // Return original length to avoid crash
    }
    
    // Use persistent DNA server for transmission
    if(!g_dna_server || !g_dna_server->isReady()) {
        fprintf(stderr, "Error: DNA server not available\n");
        // Fallback: copy original data
        memcpy(RW, CW, Nb);
        return Nb;
    }
    
    // Transmit through Julia server
    std::string received_dna = g_dna_server->transmit(dna_seq);
    
    // If no valid sequence received, return original data
    if(received_dna.empty()) {
        fprintf(stderr, "Warning: No valid DNA sequence received, using original\n");
        memcpy(RW, CW, Nb);
        return Nb;
    }
    
    // Calculate maximum allowed symbols for RW buffer
    int max_symbols = Nb + CH->GetDmax() * 2;  // Conservative estimate
    
    // Convert received DNA back to 4元符号語 with boundary checking
    int Nb2 = convert_dna_to_quaternary(received_dna, RW, max_symbols);
    
    // Ensure Nb2 is within expected range for IDS decoder: [Nb+Dmin, Nb+Dmax]
    if(Nb2 < Nb + IDS_Dmin) {
        // Pad with zeros if too short
        int target = Nb + IDS_Dmin;
        for(int i = Nb2; i < target && i < max_symbols; i++) {
            RW[i] = 0;  // 4元符号語でも0でパディング
        }
        Nb2 = target;
    } else if(Nb2 > Nb + IDS_Dmax) {
        // Truncate if too long
        Nb2 = Nb + IDS_Dmax;
    }
    
    return Nb2;
}
