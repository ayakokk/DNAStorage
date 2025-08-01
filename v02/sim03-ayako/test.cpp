#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <algorithm>
#include <stdexcept>

#include "InnerCodebook.hpp"
#include "ChannelMatrix.hpp"
#include "IDSchannel.hpp"
#include "SLFBAdec.hpp"

int main(){
  printf("====================================================\n");
  printf("SetGXNew と SetGENew の動作テスト\n");
  printf("====================================================\n");

  // パラメータ設定
  double Pi = 0.05;  // 挿入確率
  double Pd = 0.05;  // 削除確率 
  double Ps = 0.05;  // 置換確率
  int N = 12;        // ブロック長

  try {
    printf("\n--- 1. コンポーネント初期化 ---\n");
    
    // InnerCodebook初期化
    printf("InnerCodebook読み込み中...\n");
    InnerCodebook *ICB = new InnerCodebook("ICB/ex01/cb.txt", 3, 8, 2);
    printf("読み込み完了: Nu=%d, Q=%d\n", ICB->Get_Nu(), ICB->Get_numCW());

    // ChannelMatrix初期化
    printf("ChannelMatrix初期化中...\n");
    ChannelMatrix *ECM = new ChannelMatrix(ICB->Get_numCW(), ICB->Get_numCW());
    printf("ChannelMatrix初期化完了\n");

    // IDSchannel初期化 (DNA channelモード)
    printf("IDSchannel初期化中...\n");
    IDSchannel *CH = new IDSchannel(N);
    printf("IDSchannel初期化完了 (DNA channel mode)\n");

    printf("\n--- 2. SLFBAdec初期化（SetGXNew, SetGENew含む） ---\n");
    SLFBAdec *decoder = new SLFBAdec(ICB, ECM, CH);
    printf("SLFBAdec初期化完了\n");

    printf("\n--- 3. 基本パラメータの確認 ---\n");
    int Nu = ICB->Get_Nu();
    int Q = ICB->Get_numCW();
    printf("Nu (ブロック長): %d\n", Nu);
    printf("Q (コードワード数): %d\n", Q);
    printf("デコーダーの初期化が完了しました\n");
    
    printf("注意: GetGXNewとGetGENewはprivateメソッドのため、直接テストできません\n");

    printf("\n--- 4. 確率テーブル出力テスト ---\n");
    printf("確率テーブルをファイルに出力します...\n");
    
    // 確率テーブルを出力
    decoder->exportAllProbabilityTables("probability_tables");
    
    printf("確率テーブル出力完了\n");
    printf("出力先: probability_tables/ ディレクトリ\n");

    printf("\n--- 5. テスト完了 ---\n");
    printf("SLFBAdeckの初期化が正常に完了しました。\n");

    // クリーンアップ
    delete decoder;
    delete CH;
    delete ECM;
    delete ICB;

  } catch (const std::exception& e) {
    printf("エラーが発生しました: %s\n", e.what());
    return 1;
  } catch (...) {
    printf("未知のエラーが発生しました\n");
    return 1;
  }

  printf("\n====================================================\n");
  printf("テスト終了\n");
  printf("====================================================\n");
  
  return 0;
}