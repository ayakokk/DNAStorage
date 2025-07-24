#include <iostream>
#include <string>
#include <vector>
#include "func.hpp"
#include "ChannelMatrix.hpp"
#include "InnerCodebook.hpp"
#include "IDSchannel.hpp"
#include "SLFBAdec.hpp"

int main() {
    std::cout << "====================================================" << std::endl;
    std::cout << "格子計算デバッグテスト" << std::endl;
    std::cout << "特定のケース: 6 110000 110010" << std::endl;
    std::cout << "====================================================" << std::endl;

    try {
        // --- 1. コンポーネント初期化 ---
        std::cout << "\n--- 1. コンポーネント初期化 ---" << std::endl;
        
        // InnerCodebook初期化
        std::cout << "InnerCodebook読み込み中..." << std::endl;
        InnerCodebook *ICB = new InnerCodebook("ICB/ex01/cb.txt", 3, 8, 2);
        std::cout << "読み込み完了: Nu=" << ICB->Get_Nu() << ", Q=" << ICB->Get_numCW() << std::endl;
        
        // ChannelMatrix初期化
        std::cout << "ChannelMatrix初期化中..." << std::endl;
        ChannelMatrix *ECM = new ChannelMatrix(ICB->Get_numCW(), ICB->Get_numCW());
        std::cout << "ChannelMatrix初期化完了" << std::endl;
        
        // IDSchannel初期化
        std::cout << "IDSchannel初期化中..." << std::endl;
        IDSchannel *CH = new IDSchannel(12, 0.05, 0.05, 0.05);  // N=12, Pi=Pd=Ps=0.05
        std::cout << "IDSchannel初期化完了: Pi=" << CH->GetPi() << ", Pd=" << CH->GetPd() << ", Ps=" << CH->GetPs() << std::endl;

        // --- 2. SLFBAdec初期化 ---
        std::cout << "\n--- 2. SLFBAdec初期化 ---" << std::endl;
        SLFBAdec *decoder = new SLFBAdec(ICB, ECM, CH);
        std::cout << "SLFBAdec初期化完了" << std::endl;

        // --- 3. 特定ケースの格子計算デバッグ ---
        std::cout << "\n--- 3. 格子計算デバッグ ---" << std::endl;
        
        // 対象パラメータ: 6 110000 110010
        int Nu2 = 6;
        long y = 0b110000;   // 110000 (binary) = 48 (decimal)
        long xi = 18;        // 110010 (binary) = 50 (decimal), but we need to find the actual codeword index
        
        // まず、110010が何番目のコードワードかを確認
        std::cout << "コードワード検索..." << std::endl;
        unsigned char target_cw[6] = {1,1,0,0,1,0}; // 110010
        int target_xi = -1;
        
        for(int i = 0; i < ICB->Get_numCW(); i++){
            unsigned char cw[6];
            ICB->Get_CW(cw, i);
            bool match = true;
            for(int j = 0; j < 6; j++){
                if(cw[j] != target_cw[j]){
                    match = false;
                    break;
                }
            }
            if(match){
                target_xi = i;
                break;
            }
        }
        
        if(target_xi == -1){
            std::cout << "エラー: コードワード110010が見つかりません" << std::endl;
            return 1;
        }
        
        std::cout << "コードワード110010はインデックス " << target_xi << " です" << std::endl;
        
        // 送信コードワードの取得
        unsigned char X[6];
        ICB->Get_CW(X, target_xi);
        long x = VectToLong(X, 6);
        
        std::cout << "送信: ";
        for(int i = 0; i < 6; i++) std::cout << (int)X[i];
        std::cout << " (binary), " << x << " (decimal)" << std::endl;
        std::cout << "受信: ";
        unsigned char Y[6];
        LongToVect(Y, y, Nu2);
        for(int i = 0; i < Nu2; i++) std::cout << (int)Y[i];
        std::cout << " (binary), " << y << " (decimal)" << std::endl;
        
        // 格子計算デバッグの実行
        std::string filename = "lattice_debug_6_110000_110010.txt";
        std::cout << "格子計算を実行中..." << std::endl;
        decoder->debugLatticeCalculation(y, x, Nu2, 6, filename.c_str());
        
        // --- 4. 結果の確認 ---
        std::cout << "\n--- 4. 結果確認 ---" << std::endl;
        std::cout << "格子計算の詳細は " << filename << " に保存されました" << std::endl;
        
        // 実際のGXNew値との比較
        std::cout << "\n実際のGXNewテーブルとの比較:" << std::endl;
        for(int error_state = 0; error_state < 4; error_state++){
            // Note: GetGXNewはprivateなので、CalcPyxNewを直接呼び出し
            // double prob = decoder->CalcPyxNew(y, x, Nu2, 6, error_state);
            std::cout << "Error state " << error_state << ": (CalcPyxNewで計算)" << std::endl;
        }

        // --- 5. クリーンアップ ---
        delete decoder;
        delete CH;
        delete ECM;
        delete ICB;
        
        std::cout << "\n====================================================" << std::endl;
        std::cout << "テスト完了" << std::endl;
        std::cout << "====================================================" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "エラー: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}