#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>
#include <sys/wait.h>
#include <string>
#include <random>

#include "func.hpp"
#include "ChannelMatrix.hpp"
#include "InnerCodebook.hpp"
#include "IDSchannel.hpp"
#include "SLFBAdec.hpp"

// DNA channel simulation functions
std::string convert_bits_to_dna(const unsigned char* bits, int len);
int transmit_dna_channel(unsigned char* RW, const unsigned char* CW, int Nb, class IDSchannel* CH);
int convert_dna_to_bits(const std::string& dna_result, unsigned char* bits, int max_bits);

#define BSIZE 8192
#define OutListSize 3
#define WCmax 1
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
                   "/Users/ayako/research/calc_prob/DNArSim-main/simulator/dna_server.jl",
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
        if(is_ready && write_fp) {
            fprintf(write_fp, "EXIT\n");
            fflush(write_fp);
        }
        
        if(write_fp) {
            fclose(write_fp);
            write_fp = nullptr;
        }
        
        if(read_fp) {
            fclose(read_fp);
            read_fp = nullptr;
        }
        
        if(julia_pid > 0) {
            // Wait for Julia process to terminate
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
  int Rho;         // run-length    [constraint.txt]
  int ell,Delta;   // local-balance [constraint.txt]
  int N;           // block length (symbols)
  int Nb,Nb2;      // block length & recv length (bits)
  int Q, Nu;       // numCW, symbol-len [ICB]
  int seed;
  double Pi,Pd,Ps; // IDS prob
  char *fn;
  char *fncb    = new char [BSIZE];
  char *fnconst = new char [BSIZE];
  char *fncm    = new char [BSIZE];
  if(argc!=7){
    fprintf(stderr,"Usage: %s <ICB_dir> <N> <Pi> <Pd> <Ps> <seed|-1>\n",argv[0]);
    return 1;
  } // if
  fn   =      argv[1];
  N    = atoi(argv[2]);
  Pi   = atof(argv[3]);
  Pd   = atof(argv[4]);
  Ps   = atof(argv[5]);
  seed = atoi(argv[6]);
  if(seed==-1) seed = (int)time(NULL);
  srandom(seed);
  assert(N>0);
  assert(Pi>=0.0 && Pi<0.5);
  assert(Pd>=0.0 && Pd<0.5);
  assert(Ps>=0.0 && Ps<0.5);
  snprintf(fncb,   BSIZE,"%s/cb.txt",        fn);  // inner codebook (in)
  snprintf(fnconst,BSIZE,"%s/constraint.txt",fn);  // constraints (in)
  snprintf(fncm,   BSIZE,"%s/EncCM.bin",     fn);  // encoding channel matrix (in)
  ReadConstraints(fnconst, &Rho, &ell, &Delta);
  class InnerCodebook *ICB = new class InnerCodebook(fncb,Rho,ell,Delta);
  Q    = ICB->Get_numCW();
  Nu   = ICB->Get_Nu();
  Nb   = N*Nu;
  printf("# Q=%d N=%d Nu=%d Nb=%d (Pi,Pd,Ps)=(%e,%e,%e) [%d]\n",Q,N,Nu,Nb,Pi,Pd,Ps,seed);
  printf("# ICB:   %s\n",fncb);
  printf("# Const: %s\n",fnconst);
  printf("# EncCM: %s\n",fncm);
  class ChannelMatrix *ECM = new class ChannelMatrix(fncm);
  class IDSchannel    *CH  = new class IDSchannel(Nb,Pi,Pd,Ps);
  class SLFBAdec      *DEC = new class SLFBAdec(ICB,ECM,CH);
  class ChannelMatrix *DCM = new class ChannelMatrix(Q,(int)pow(Q,OutListSize));
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
  
  // Initialize DNA channel server
  g_dna_server = new DNAChannelServer();
  if(!g_dna_server->initialize()) {
    fprintf(stderr, "Error: Failed to initialize DNA channel server\n");
    fprintf(stderr, "Using fallback mode (copy original data)\n");
  }
  
  //-----
  for(wc=1;wc<=WCmax;wc++){
    RandVect(IW,N,0,Q-1);
    ICB->Encode(CW,IW,N);
    // Use DNA channel simulation instead of IDS channel
    Nb2 = transmit_dna_channel(RW,CW,Nb,CH);
    DEC->Decode(Pout,RW,Nb2,IW);
    HardDecision(DW,(const double **)Pout,N,Q);
    OutputConv(DWL,(const double **)Pout,N,Q);
    for(int i=0;i<N;i++) DCM->countup(IW[i],DWL[i]);
    ec = HammingDist(IW,DW,N);
    es += ec;
    ecmax = max(ec,ecmax);

    if(wc%1000==0 || wc==WCmax){
      printf("%04d %ld/%ld %ld %e : %e %e %e\n",
	     wc,es,(long)wc*N,ecmax,(double)es/(wc*N), DCM->Hx(), DCM->Hxy(), DCM->Ixy());
    } // if wc 

    //(dbg)
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

// Convert binary bits to DNA sequence (2 bits per base: 00=A, 01=C, 10=G, 11=T)

/**
 * @brief ビット列をランダムなビットと組み合わせてDNA配列に変換する
 * * 各入力ビット `bits[i]` を下位ビットとし、ランダムに生成したビットを
 * 上位ビットとして2ビットのペアを生成し、対応するDNA塩基に変換する。
 * (00:A, 01:C, 10:G, 11:T)
 * * @param bits 変換するビット(0または1)を格納したunsigned char型の配列
 * @param len 配列の長さ
 * @return std::string 変換後のDNA配列
 */

std::string convert_bits_to_dna(const unsigned char* bits, int len) {
    std::string dna = "";
    
    //乱数生成器
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1); // 0または1のランダムなビットを生成

    for(int i = 0; i < len; i += 2) {
        //　ランダムなビットを生成
        unsigned char random_bit = dis(gen);

        //入力ビットを取得
        unsigned char input_bit = bits[i];
        // ランダムビットを入力ビットでペアを作成
        unsigned char pair = (random_bit << 1) | input_bit;
        switch(pair) {
            case 0: dna += "A"; break;  // 00
            case 1: dna += "C"; break;  // 01
            case 2: dna += "G"; break;  // 10
            case 3: dna += "T"; break;  // 11
        }
    }
    return dna;
}

// Convert DNA sequence back to binary bits
int convert_dna_to_bits(const std::string& dna_result, unsigned char* bits, int max_bits) {
    int bit_len = 0;
    for(int i = 0; i < dna_result.length() && bit_len + 1 < max_bits; i++) {
        char base = dna_result[i];
        switch(base) {
            case 'A': case 'a':
                if(bit_len + 1 < max_bits) {
                    bits[bit_len++] = 0; bits[bit_len++] = 0;
                }
                break;
            case 'C': case 'c':
                if(bit_len + 1 < max_bits) {
                    bits[bit_len++] = 0; bits[bit_len++] = 1;
                }
                break;
            case 'G': case 'g':
                if(bit_len + 1 < max_bits) {
                    bits[bit_len++] = 1; bits[bit_len++] = 0;
                }
                break;
            case 'T': case 't':
                if(bit_len + 1 < max_bits) {
                    bits[bit_len++] = 1; bits[bit_len++] = 1;
                }
                break;
            default:
                // Skip unknown characters
                break;
        }
    }
    return bit_len;
}

// Main DNA channel transmission function using persistent server
int transmit_dna_channel(unsigned char* RW, const unsigned char* CW, int Nb, class IDSchannel* CH) {
    // Calculate dynamic Dmin/Dmax based on block length
    int drift_range = (Nb >= 100) ? 3 : (Nb >= 50) ? 2 : 1;
    const int IDS_Dmin = -drift_range;
    const int IDS_Dmax = drift_range;
    
    // Convert bits to DNA sequence
    std::string dna_seq = convert_bits_to_dna(CW, Nb);
    
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
    
    // Calculate maximum allowed bits for RW buffer
    int max_bits = Nb + CH->GetDmax() * 2;  // Conservative estimate
    
    // Convert received DNA back to bits with boundary checking
    int Nb2 = convert_dna_to_bits(received_dna, RW, max_bits);
    
    // Ensure Nb2 is within expected range for IDS decoder: [Nb+Dmin, Nb+Dmax]
    if(Nb2 < Nb + IDS_Dmin) {
        // Pad with zeros if too short
        int target = Nb + IDS_Dmin;
        for(int i = Nb2; i < target && i < max_bits; i++) {
            RW[i] = 0;
        }
        Nb2 = target;
    } else if(Nb2 > Nb + IDS_Dmax) {
        // Truncate if too long
        Nb2 = Nb + IDS_Dmax;
    }
    
    return Nb2;
}
