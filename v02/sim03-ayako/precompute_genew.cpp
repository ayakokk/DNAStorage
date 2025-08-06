#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "InnerCodebook.hpp"
#include "ChannelMatrix.hpp"
#include "IDSchannel.hpp"
#include "SLFBAdec.hpp"

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s <ICB_dir> <N> <output_dir>\n", argv[0]);
        printf("Example: %s ICB/ex01/ 3 precomputed_genew/\n", argv[0]);
        printf("Note: Pi,Pd,Ps parameters are automatically set to optimal values (0.05)\n");
        return 1;
    }

    printf("====================================================\n");
    printf("GENew Precomputation Tool\n");
    printf("====================================================\n");

    // パラメータ解析
    const char *icb_dir = argv[1];
    int N = atoi(argv[2]);
    const char *output_dir = argv[3];
    
    // 最適値を自動設定（IDSchannelと同じ値）
    double Pi = 0.05;
    double Pd = 0.05;
    double Ps = 0.05;

    printf("Parameters:\n");
    printf("  ICB directory: %s\n", icb_dir);
    printf("  Block length N: %d\n", N);
    printf("  Channel: Pi=%.6f, Pd=%.6f, Ps=%.6f (auto-optimized)\n", Pi, Pd, Ps);
    printf("  Output directory: %s\n", output_dir);

    auto total_start = std::chrono::high_resolution_clock::now();

    try {
        printf("\n--- 1. Component Initialization ---\n");
        
        // InnerCodebook初期化
        std::string cb_path = std::string(icb_dir) + "/cb.txt";
        printf("Loading InnerCodebook: %s\n", cb_path.c_str());
        InnerCodebook *ICB = new InnerCodebook(cb_path.c_str(), N, N*6, 2);
        printf("Loaded: Nu=%d, Q=%d\n", ICB->Get_Nu(), ICB->Get_numCW());

        // ChannelMatrix初期化
        printf("Initializing ChannelMatrix...\n");
        ChannelMatrix *ECM = new ChannelMatrix(ICB->Get_numCW(), ICB->Get_numCW());

        // IDSchannel初期化
        printf("Initializing IDSchannel...\n");
        IDSchannel *CH = new IDSchannel(N*ICB->Get_Nu());

        printf("\n--- 2. SLFBAdec Initialization ---\n");
        SLFBAdec *decoder = new SLFBAdec(ICB, ECM, CH);

        printf("\n--- 3. Precomputation Process ---\n");
        
        // 出力ディレクトリ作成
        std::string mkdir_cmd = "mkdir -p " + std::string(output_dir);
        system(mkdir_cmd.c_str());

        // 全てのlyについて事前計算実行
        int Nu = ICB->Get_Nu();
        int Nu2min = Nu - (int)ceil((double)Nu * Pd) - 2;
        int Nu2max = Nu + (int)ceil((double)Nu * Pi) + 2;
        Nu2max = std::min(Nu2max, Nu * 2);
        Nu2min = std::max(Nu2min, 0);

        printf("Computing range: ly=%d to %d\n", Nu2min, Nu2max);

        for (int ly = Nu2min; ly <= Nu2max; ly++) {
            auto ly_start = std::chrono::high_resolution_clock::now();
            
            printf("\n--- Computing ly=%d ---\n", ly);
            long ly2p = (long)pow(2, ly);
            int Q = ICB->Get_numCW();
            
            printf("Entries to compute: %ld y-values × %d codewords × 4 error states = %ld total\n", 
                   ly2p, Q, ly2p * Q * 4);

            // 事前計算実行
            decoder->PrecomputeGENewForLy(ly);

            // バイナリファイルに保存
            std::string filename = std::string(output_dir) + "/GENew_ly" + std::to_string(ly) + ".bin";
            decoder->SaveGENewToFile(ly, filename.c_str());

            auto ly_end = std::chrono::high_resolution_clock::now();
            auto ly_duration = std::chrono::duration_cast<std::chrono::seconds>(ly_end - ly_start);
            printf("ly=%d completed in %lld seconds\n", ly, ly_duration.count());
        }

        printf("\n--- 4. Creating Metadata ---\n");
        
        // メタデータファイル作成
        std::string metadata_file = std::string(output_dir) + "/metadata.txt";
        std::ofstream meta(metadata_file);
        if (meta.is_open()) {
            meta << "# GENew Precomputed Data Metadata\n";
            meta << "# Generated on: " << __DATE__ << " " << __TIME__ << "\n";
            meta << "Nu=" << Nu << "\n";
            meta << "Q=" << ICB->Get_numCW() << "\n";
            meta << "Pi=" << std::scientific << std::setprecision(6) << Pi << "\n";
            meta << "Pd=" << std::scientific << std::setprecision(6) << Pd << "\n";
            meta << "Ps=" << std::scientific << std::setprecision(6) << Ps << "\n";
            meta << "Nu2min=" << Nu2min << "\n";
            meta << "Nu2max=" << Nu2max << "\n";
            meta << "k_mer_length=4\n";
            
            for (int ly = Nu2min; ly <= Nu2max; ly++) {
                std::string filename = "GENew_ly" + std::to_string(ly) + ".bin";
                meta << "file_ly" << ly << "=" << filename << "\n";
            }
            meta.close();
            printf("Metadata saved: %s\n", metadata_file.c_str());
        }

        printf("\n--- 5. Cleanup ---\n");
        delete decoder;
        delete CH;
        delete ECM;
        delete ICB;

        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::minutes>(total_end - total_start);

        printf("\n====================================================\n");
        printf("Precomputation Completed Successfully!\n");
        printf("Total time: %ld minutes\n", total_duration.count());
        printf("Output directory: %s\n", output_dir);
        printf("====================================================\n");

    } catch (const std::exception& e) {
        printf("ERROR: %s\n", e.what());
        return 1;
    }

    return 0;
}