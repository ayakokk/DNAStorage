#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"
#include "ChannelMatrix.hpp"
#include "InnerCodebook.hpp"
#include "IDSchannel.hpp"
#include "SLFBAdec.hpp"

#define BSIZE 8192
#define OutListSize 3
#define WCmax_EXP 20000  // 実験用に短縮

void OutputConv(long *DWL, const double **P, int N, int Q);
long PdistArgMaxLong(const double *P, int Q, int LS);

// 実験タイプの定義
enum ExperimentType {
    IDS_CHANNEL = 1,  // Figure 6: Pi=Pd=Ps
    ID_CHANNEL = 2,   // Figure 7: Pi=Pd, Ps=0
    DEL_CHANNEL = 3   // Figure 8: Pi=Ps=0, Pd only
};

// 重み設定の定義
struct WeightConfig {
    int weights[4];
    const char* name;
    const char* description;
};

//================================================================================
int main(int argc, char *argv[]){
    if(argc != 4) {
        fprintf(stderr,"Usage: %s <ICB_dir> <N> <seed|-1>\n", argv[0]);
        fprintf(stderr,"Example: %s ./ICB/ex01 50 1\n", argv[0]);
        return 1;
    }
    
    char *fn = argv[1];
    int N = atoi(argv[2]);
    int seed = atoi(argv[3]);
    
    if(seed == -1) seed = (int)time(NULL);
    srandom(seed);
    assert(N > 0);
    
    // 複数の重み設定を定義
    WeightConfig weight_configs[] = {
        {{1, 1, 1, 1}, "uniform", "Uniform (Random equivalent)"},
        {{6, 5, 4, 3}, "light", "Light weighting"},
        {{5, 4, 3, 2}, "medium", "Medium weighting"},
        {{4, 3, 2, 1}, "strong", "Strong weighting"},
        {{10, 6, 3, 1}, "extreme", "Extreme weighting"},

        {{3, 4, 5, 6}, "light_reverce", "Light_reverce weighting"},
        {{2, 3, 4, 5}, "medium_reverce", "Medium_reverce weighting"},
        {{1,2,3,4}, "strong_reverce", "Strong_reverce weighting"},
        {{1,3,6,10}, "extreme_reverce", "Extreme_reverce weighting"}
    };
    int num_weight_configs = sizeof(weight_configs) / sizeof(weight_configs[0]);
    
    // ファイルパス設定
    char *fncb    = new char [BSIZE];
    char *fnconst = new char [BSIZE];
    char *fncm    = new char [BSIZE];
    snprintf(fncb,   BSIZE, "%s/cb.txt", fn);
    snprintf(fnconst,BSIZE, "%s/constraint.txt", fn);
    snprintf(fncm,   BSIZE, "%s/EncCM.bin", fn);
    
    // 制約読み込み
    int Rho, ell, Delta;
    ReadConstraints(fnconst, &Rho, &ell, &Delta);
    
    // Inner Codebook初期化
    class InnerCodebook *ICB = new class InnerCodebook(fncb, Rho, ell, Delta);
    int Q = ICB->Get_numCW();
    int Nu = ICB->Get_Nu();
    int Nb = N * Nu;
    
    printf("# ========== Experiment Configuration ==========\n");
    printf("# Codebook: %s\n", fncb);
    printf("# Q=%d N=%d Nu=%d Nb=%d seed=%d\n", Q, N, Nu, Nb, seed);
    printf("# WCmax_EXP=%d\n", WCmax_EXP);
    printf("# Weight configurations:\n");
    for(int i = 0; i < num_weight_configs; i++) {
        printf("#   %s: [%d,%d,%d,%d] - %s\n", 
               weight_configs[i].name,
               weight_configs[i].weights[0], weight_configs[i].weights[1],
               weight_configs[i].weights[2], weight_configs[i].weights[3],
               weight_configs[i].description);
    }
    printf("# ===============================================\n\n");
    
    // エンコーディングチャネルマトリックス
    class ChannelMatrix *ECM = new class ChannelMatrix(fncm);
    
    // エラー率設定（論文に合わせて）
    double error_rates[] = {0.000, 0.001, 0.002, 0.003, 0.004, 0.005, 
                           0.006, 0.007, 0.008, 0.009, 0.010};
    int num_rates = sizeof(error_rates) / sizeof(error_rates[0]);
    
    // 実験タイプのリスト
    ExperimentType experiments[] = {IDS_CHANNEL, ID_CHANNEL, DEL_CHANNEL};
    const char* exp_names[] = {"IDS", "ID", "DEL"};
    int num_experiments = sizeof(experiments) / sizeof(experiments[0]);
    
    // メモリ割り当て
    int *dbgDR = new int [Nb+1];
    int *IW = new int [N];
    unsigned char *CW = new unsigned char [Nb];
    unsigned char *RW = new unsigned char [Nb + 12]; // 最大ドリフト想定
    int *DW = new int [N];
    long *DWL = new long [N];
    double **Pout = new double * [N];
    for(int i = 0; i < N; i++) Pout[i] = new double [Q];
    
    // 各実験タイプについて実行
    for(int exp_idx = 0; exp_idx < num_experiments; exp_idx++) {
        ExperimentType exp_type = experiments[exp_idx];
        
        printf("# ========== %s Channel Experiment ==========\n", exp_names[exp_idx]);
        
        // 各重み設定について実験
        for(int weight_idx = 0; weight_idx < num_weight_configs; weight_idx++) {
            WeightConfig& config = weight_configs[weight_idx];
            
            printf("## Weight Config: %s [%d,%d,%d,%d]\n", 
                   config.name, config.weights[0], config.weights[1], 
                   config.weights[2], config.weights[3]);
            printf("# Error_Rate\tMutual_Info\tEntropy_X\tEntropy_XY\tSER\n");
            
            // 各エラー率について実験
            for(int rate_idx = 0; rate_idx < num_rates; rate_idx++) {
                double error_rate = error_rates[rate_idx];
                double Pi, Pd, Ps;
                
                // チャネルタイプに応じてエラー率設定
                switch(exp_type) {
                    case IDS_CHANNEL:
                        Pi = Pd = Ps = error_rate;
                        break;
                    case ID_CHANNEL:
                        Pi = Pd = error_rate;
                        Ps = 0.0;
                        break;
                    case DEL_CHANNEL:
                        Pi = Ps = 0.0;
                        Pd = error_rate;
                        break;
                }
                
                // チャネルとデコーダの初期化
                class IDSchannel *CH = new class IDSchannel(Nb, Pi, Pd, Ps);
                class SLFBAdec *DEC = new class SLFBAdec(ICB, ECM, CH);
                class ChannelMatrix *DCM = new class ChannelMatrix(Q, (int)pow(Q, OutListSize));
                
                // シミュレーション実行
                long es = 0;
                long ecmax = 0;
                
                for(int wc = 1; wc <= WCmax_EXP; wc++) {
                    // 情報語生成（重み付き）
                    if(weight_idx == 0) {
                        // uniform の場合はランダム
                        RandVect(IW, N, 0, Q-1);
                    } else {
                        // 重み付き選択
                        FixedVect(IW, N, 0, Q-1, config.weights);
                    }
                    
                    // エンコーディング・チャネル・デコーディング
                    ICB->Encode(CW, IW, N);
                    int Nb2 = CH->transmit(RW, CW);
                    DEC->Decode(Pout, RW, Nb2, IW);
                    HardDecision(DW, (const double **)Pout, N, Q);
                    OutputConv(DWL, (const double **)Pout, N, Q);
                    
                    // 統計更新
                    for(int i = 0; i < N; i++) DCM->countup(IW[i], DWL[i]);
                    long ec = HammingDist(IW, DW, N);
                    es += ec;
                    ecmax = max(ec, ecmax);
                }
                
                // 結果出力
                double ser = (double)es / (WCmax_EXP * N);
                printf("%.3f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
                       error_rate, DCM->Ixy(), DCM->Hx(), DCM->Hxy(), ser);
                
                // 進行状況表示
                if(rate_idx % 3 == 0) {
                    fprintf(stderr, "  %s-%s: %.3f completed\n", 
                           exp_names[exp_idx], config.name, error_rate);
                }
                
                // クリーンアップ
                delete CH;
                delete DEC;
                delete DCM;
            }
            printf("\n");
        }
        printf("\n");
    }
    
    // 最終クリーンアップ
    delete ICB;
    delete ECM;
    delete [] dbgDR;
    delete [] IW;
    delete [] CW;
    delete [] RW;
    delete [] DW;
    delete [] DWL;
    delete [] fncb;
    delete [] fnconst;
    delete [] fncm;
    for(int i = 0; i < N; i++) delete [] Pout[i];
    delete [] Pout;
    
    printf("# ========== Experiment Completed ==========\n");
    return 0;
}

//================================================================================
void OutputConv(long *DWL, const double **P, int N, int Q){
    for(int i = 0; i < N; i++){
        DWL[i] = PdistArgMaxLong(P[i], Q, OutListSize);
    }
}

//================================================================================
long PdistArgMaxLong(const double *P, int Q, int LS){
    assert(Q > 0 && LS > 0 && LS <= Q);
    long v, val = 0;
    double *PX = new double [Q];
    memcpy(PX, P, sizeof(double) * Q);
    for(int i = 0; i < Q; i++) assert( P[i] >= 0.0 && P[i] <= 1.0 );
    for(int i = 0; i < LS; i++){
        v = argmax(PX, Q);
        val = val * Q + v;
        PX[v] = -1.0;
    }
    delete [] PX;
    return val;
}