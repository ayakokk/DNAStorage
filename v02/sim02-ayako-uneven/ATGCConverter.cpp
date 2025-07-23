#include "ATGCConverter.hpp"
#include <cassert>
#include <cstring>

// 単体変換関数（前回と同じ）
char binary_to_atgc(unsigned char b1, unsigned char b2) {
    if (b1 == 0 && b2 == 0) return 'A';  // 00 → A
    if (b1 == 0 && b2 == 1) return 'T';  // 01 → T  
    if (b1 == 1 && b2 == 0) return 'G';  // 10 → G
    if (b1 == 1 && b2 == 1) return 'C';  // 11 → C
    return 'A'; // エラー時のデフォルト
}

void atgc_to_binary(char atgc, unsigned char &b1, unsigned char &b2) {
    switch(atgc) {
        case 'A': b1 = 0; b2 = 0; break;  // A → 00
        case 'T': b1 = 0; b2 = 1; break;  // T → 01
        case 'G': b1 = 1; b2 = 0; break;  // G → 10
        case 'C': b1 = 1; b2 = 1; break;  // C → 11
        default:  b1 = 0; b2 = 0; break;  // エラー時
    }
}

// 配列変換関数
void binary_array_to_atgc(char *atgc_output, const unsigned char *binary_input, int binary_len) {
    assert(binary_len % 2 == 0);  // 偶数長であることを確認
    
    int atgc_len = binary_len / 2;
    for(int i = 0; i < atgc_len; i++) {
        atgc_output[i] = binary_to_atgc(binary_input[i*2], binary_input[i*2+1]);
    }
    atgc_output[atgc_len] = '\0';  // 終端文字を追加
}

int atgc_array_to_binary(unsigned char *binary_output, const char *atgc_input, int atgc_len) {
    for(int i = 0; i < atgc_len; i++) {
        unsigned char b1, b2;
        atgc_to_binary(atgc_input[i], b1, b2);
        binary_output[i*2] = b1;
        binary_output[i*2+1] = b2;
    }
    return atgc_len * 2;  // バイナリ配列の長さを返す
}