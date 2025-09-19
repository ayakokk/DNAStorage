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
#define OutListSize 3
#define WCmax 1  // Dec3 BER evaluation with 1000 test words

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
                    printf("# [INFO] DNA channel server initialized successfully\n");
                    return true;
                }
            }
            
            fprintf(stderr, "Error: DNA server failed to initialize\n");
            return false;
        }
        return false;
    }
    
    std::string transmit(const std::string& dna_sequence) {
        if(!is_ready || !write_fp || !read_fp) {
            fprintf(stderr, "Error: DNA server not ready\n");
            return "";
        }
        
        // Send DNA sequence to Julia server
        fprintf(write_fp, "%s\n", dna_sequence.c_str());
        fflush(write_fp);
        
        // Read response
        char buffer[8192];
        if(fgets(buffer, sizeof(buffer), read_fp)) {
            // Remove newline
            std::string result(buffer);
            if(!result.empty() && result.back() == '\n') {
                result.pop_back();
            }
            
            // Check if it's an error message
            if(result.find("ERROR:") == 0) {
                fprintf(stderr, "Julia server error: %s\n", result.c_str());
                return "";
            }
            
            // Check if it's a RESULT: response
            if(result.find("RESULT:") == 0) {
                return result.substr(7); // Remove "RESULT:" prefix
            }
            
            // Otherwise return as-is
            return result;
        }
        
        return "";
    }
    
    void shutdown() {
        if(write_fp) {
            fprintf(write_fp, "EXIT\n");
            fflush(write_fp);
            fclose(write_fp);
            write_fp = nullptr;
        }
        
        if(read_fp) {
            fclose(read_fp);
            read_fp = nullptr;
        }
        
        if(julia_pid != -1) {
            // Wait for Julia process to terminate
            int status;
            waitpid(julia_pid, &status, 0);
            julia_pid = -1;
        }
        
        is_ready = false;
    }
    
    bool isReady() const { return is_ready; }
};

// Global DNA server instance
DNAChannelServer* g_dna_server = nullptr;

//================================================================================
int main(int argc, char *argv[]){
  // 標準出力のバッファリングを無効にする（常に即時フラッシュ）
  setvbuf(stdout, NULL, _IONBF, 0); 
  // 標準エラー出力も同様に
  setvbuf(stderr, NULL, _IONBF, 0);
  int Rho;         // run-length    [constraint.txt]
  int ell,Delta;   // local-balance [constraint.txt]
  int N;           // block length (symbols)
  int Nb,Nb2;      // block length & recv length (bits)
  int Q, Nu;       // numCW, symbol-len [ICB]
  int seed;
  double sparse_threshold;
  char *fn;
  char *fncb    = new char [BSIZE];
  char *fnconst = new char [BSIZE];
  char *fncm    = new char [BSIZE];
  if(argc!=5){
    fprintf(stderr,"Usage: %s <ICB_dir> <N> <seed|-1> <sparse_threshold>\n",argv[0]);
    fprintf(stderr,"  sparse_threshold: 事前計算で枝刈りする確率閾値 (e.g., 1e-6)\n");
    return 1;
  } // if
  fn               =      argv[1];
  N                = atoi(argv[2]);
  seed             = atoi(argv[3]);
  sparse_threshold = atof(argv[4]);
  if(seed==-1) seed = (int)time(NULL);
  srandom(seed);
  assert(N>0);
  snprintf(fncb,   BSIZE,"%s/cb.txt",        fn);  // inner codebook (in)
  snprintf(fnconst,BSIZE,"%s/constraint.txt",fn);  // constraints (in)
  snprintf(fncm,   BSIZE,"%s/EncCM.bin",     fn);  // encoding channel matrix (in)
  ReadConstraints(fnconst, &Rho, &ell, &Delta);
  class InnerCodebook *ICB = new class InnerCodebook(fncb,Rho,ell,Delta);
  Q    = ICB->Get_numCW();
  Nu   = ICB->Get_Nu();
  Nb   = N*Nu;
  printf("# Q=%d N=%d Nu=%d Nb=%d DNA_Channel [%d]\n",Q,N,Nu,Nb,seed);
  printf("# ICB:   %s\n",fncb);
  printf("# Const: %s\n",fnconst);
  printf("# EncCM: %s\n",fncm);
  
  // Initialize DNA channel server
  g_dna_server = new DNAChannelServer();
  if(!g_dna_server->initialize()){
    fprintf(stderr, "Error: Failed to initialize DNA channel server\n");
    delete g_dna_server;
    return 1;
  }
  
  class ChannelMatrix *ECM = new class ChannelMatrix(fncm);
  class IDSchannel    *CH  = new class IDSchannel(Nb,0.0,0.0,0.0); // dummy init
  // Initialize Dec3 decoder for BER performance test
  class SLFBAdec      *DEC = new class SLFBAdec(ICB,ECM,CH);
  // [DEBUG] DCM size calculation: Q=%d, OutListSize=%d, pow(Q,OutListSize)=%.0f
  printf("# [DEBUG] DCM size calculation: Q=%d, OutListSize=%d, pow(Q,OutListSize)=%.0f\n", Q, OutListSize, pow(Q,OutListSize));
  printf("# [DEBUG] This would require approximately %.2f GB of memory - DISABLED\n", 
         (double)Q * pow(Q,OutListSize) * sizeof(double) / (1024*1024*1024));
  printf("# [DEBUG] Skipping DCM allocation for Dec3 mode\n");
  class ChannelMatrix *DCM = nullptr;  // Dec3では巨大すぎるDCMを無効化
  int        *dbgDR = new int [Nb+1];                         // (dbg)drift vector
  int           *IW = new int [N];                            // information word
  unsigned char *CW = new unsigned char [Nb]();                 // codeword (zero-initialized)
  unsigned char *RW = new unsigned char [Nb + CH->GetDmax()](); // received word (zero-initialized)
  // Dec3 decoder variables for BER calculation
  int           *DW = new int [N];                            // decoded word
  long         *DWL = new long [N];                           // (Pout->list->long)
  double **Pout = new double * [N];
  for(int i=0;i<N;i++) Pout[i] = new double [Q];
  int wc;
  long ec,ecmax=0,es=0;  // Error counting for BER calculation
  
  // ▼▼▼【ここからコードを挿入】▼▼▼
  const char* decoder_mode = getenv("DECODER_MODE");
  if (decoder_mode == nullptr) {
      decoder_mode = "DEC3"; // デフォルトをDEC3に
  }

  // DEC3モードの場合のみ、スパース遷移行列の準備（キャッシュ機能付き）
  if (strcmp(decoder_mode, "DEC3") == 0) {
      // キャッシュファイル名を生成（閾値に応じてファイル名が変わる）
      char cache_filename[256];
      snprintf(cache_filename, sizeof(cache_filename), "sparse_transitions_thresh_%.0e.cache", sparse_threshold);

      printf("# [INFO] Setting up sparse transition matrices for DEC3 performance optimization...\n");

      // 1. まずキャッシュの読み込みを試す
      if (!DEC->LoadSparseTransitions(cache_filename)) {
          // 2. 失敗した場合のみ、事前計算を実行
          printf("# [INFO] No cache found. Precomputing sparse transition matrices...\n");
          // PrecomputeSparseTransitions removed - using dynamic computation

          // 3. 計算結果を次のために保存
          DEC->SaveSparseTransitions(cache_filename);
          printf("# [INFO] Cache saved for future use.\n");
      } else {
          printf("# [INFO] Using cached sparse transition matrices. Precomputation skipped.\n");
      }

      //【確認用】作成したスパース行列をファイルに出力
      DEC->ExportSparseTransitions("sparse_transitions_result.txt");

      // 【監視システム】リアルタイムキャッシュ監視を開始
      DEC->StartCacheMonitoring("cache_snapshot.txt", 30); // 30秒間隔でスナップショット
  }
// ▲▲▲【ここまでコードを挿入】▲▲▲
  //===== Dec3 BER Performance Test Loop =====
  printf("# Starting Dec3 BER performance evaluation with %d test words\n", WCmax);
  for(wc=1;wc<=WCmax;wc++){
    // ランダムに選択
    // IWを不均一にすれば
    RandVect(IW,N,0,Q-1);
    ICB->Encode(CW,IW,N);
    Nb2 = transmit_dna_channel(RW,CW,Nb,CH);  // Use DNA channel transmission
    if(Nb2 < 0){
      fprintf(stderr, "Error: DNA transmission failed at word %d\n", wc);
      break;
    }
    
    // Print original and received codewords for comparison
    printf("# [COMPARISON] Original CW: ");
    for(int i = 0; i < Nb; i++) printf("%d", CW[i]);
    printf("\n# [COMPARISON] Received RW: ");
    for(int i = 0; i < Nb2; i++) printf("%d", RW[i]);
    printf("\n# [COMPARISON] Length change: %d -> %d\n", Nb, Nb2);
    
    // ✅ Dec3 4D Lattice Decoding (k-mer dependent with error state memory)
    // const char* decoder_mode = getenv("DECODER_MODE");
    const char* decoder_mode = "DEC3"; // Forcing Dec3 mode for this test

    if(decoder_mode && strcmp(decoder_mode, "DEC3") == 0) {
      printf("# [Dec3] Starting 4D lattice decoding (k-mer + error state memory)\n");
    }
    
    DEC->Decode(Pout,RW,Nb2,IW, sparse_threshold);  // Dec3 ultimate decoder
    HardDecision(DW,(const double **)Pout,N,Q);
    OutputConv(DWL,(const double **)Pout,N,Q);
    if(DCM != nullptr) {
      for(int i=0;i<N;i++) DCM->countup(IW[i],DWL[i]);
    }
    ec = HammingDist(IW,DW,N);
    es += ec;
    ecmax = max(ec,ecmax);

    if(wc%1000==0 || wc==WCmax){
      if(DCM != nullptr) {
        printf("%04d %ld/%ld %ld %e : %e %e %e\n",
               wc,es,(long)wc*N,ecmax,(double)es/(wc*N), DCM->Hx(), DCM->Hxy(), DCM->Ixy());
      } else {
        printf("%04d %ld/%ld %ld %e\n",
               wc,es,(long)wc*N,ecmax,(double)es/(wc*N));
      }
      printf("# [Dec3 Progress] Word %d/%d, BER = %.4f%%, Errors = %ld/%ld\n", 
             wc, WCmax, 100.0*(double)es/(wc*N), es, (long)wc*N);
    } // if wc 

    //(dbg)
    //CH->GetDR(dbgDR);
    //dbgPrint(IW,DW,CW,RW,(const double**)Pout,dbgDR,Nb,Nb2,Nu,Q);
    //printf("%04d %ld %ld/%ld %e\n",wc,ec,es,(long)wc*N,(double)es/(wc*N));
  } // for wc
  //-----
  //DCM->PrintCnt();

  // 【監視システム】キャッシュ監視を停止
  if (strcmp(decoder_mode, "DEC3") == 0) {
    DEC->StopCacheMonitoring();
    printf("# Final cache export...\n");
    DEC->ExportTransitionProbCache("final_transition_cache.txt");
  }

  // Shutdown DNA server
  if(g_dna_server){
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
  // delete [] DW;
  // delete [] DWL;
  delete [] fncb;
  delete [] fnconst;
  delete [] fncm;
  for(int i=0;i<N;i++) delete [] Pout[i];
  delete [] Pout;
  delete [] DW;
  delete [] DWL;
  return 0;
}

//================================================================================
// DNA conversion functions for 4-ary codewords
std::string convert_quaternary_to_dna(const unsigned char* symbols, int len){
  std::string dna_sequence;
  dna_sequence.reserve(len);
  
  for(int i = 0; i < len; i++){
    assert(symbols[i] >= 0 && symbols[i] <= 3);
    switch(symbols[i]){
      case 0: dna_sequence += 'A'; break;  // 0 -> A
      case 1: dna_sequence += 'C'; break;  // 1 -> C  
      case 2: dna_sequence += 'G'; break;  // 2 -> G
      case 3: dna_sequence += 'T'; break;  // 3 -> T
      default: assert(false);
    }
  }
  return dna_sequence;
}

int transmit_dna_channel(unsigned char* RW, const unsigned char* CW, int Nb, class IDSchannel* CH){
  if(!g_dna_server || !g_dna_server->isReady()){
    fprintf(stderr, "Error: DNA server not initialized\n");
    return -1;
  }
  
  // Convert 4-ary codeword to DNA sequence
  std::string dna_input = convert_quaternary_to_dna(CW, Nb);
  printf("# [DEBUG] Sending DNA sequence: %s\n", dna_input.c_str());
  
  // Transmit through Julia DNA channel  
  std::string dna_result = g_dna_server->transmit(dna_input);
  
  if(dna_result.empty()){
    fprintf(stderr, "Error: DNA transmission failed\n");
    return -1;
  }
  
  printf("# [DEBUG] Received DNA result: %s\n", dna_result.c_str());
  
  // Convert DNA result back to 4-ary symbols
  int result_length = convert_dna_to_quaternary(dna_result, RW, Nb + CH->GetDmax());
  
  printf("# [DEBUG] Converted to %d quaternary symbols\n", result_length);
  
  return result_length;
}

int convert_dna_to_quaternary(const std::string& dna_result, unsigned char* symbols, int max_symbols){
  int len = dna_result.length();
  if(len > max_symbols){
    len = max_symbols;
  }
  
  for(int i = 0; i < len; i++){
    char nucleotide = dna_result[i];
    switch(nucleotide){
      case 'A': symbols[i] = 0; break;  // A -> 0
      case 'C': symbols[i] = 1; break;  // C -> 1
      case 'G': symbols[i] = 2; break;  // G -> 2
      case 'T': symbols[i] = 3; break;  // T -> 3
      default: 
        fprintf(stderr, "Warning: Invalid nucleotide '%c' at position %d, using A\n", nucleotide, i);
        symbols[i] = 0; // Default to A
        break;
    }
  }
  
  return len;
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
	printf("[%03d] %02d %02d\n",idx,IW[idx],DW[idx]);
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
