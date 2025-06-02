#ifndef JULIACHANNEL_HPP
#define JULIACHANNEL_HPP

#include "ATGCConverter.hpp"
#include <iostream>

class JuliaChannel {
public:
    JuliaChannel(int k = 6);
    ~JuliaChannel();
    
    int transmit(unsigned char *RW, const unsigned char *CW, int Nb);
    
    // IDSchannelとの互換性のためのメソッド
    int GetN() const { return Nb; }
    int GetDmax() const { return max_expansion; }
    int GetDmin() const { return -max_contraction; }
    void print_statistics() const {
        std::cout << "=== JuliaChannel Statistics ===" << std::endl;
        std::cout << "Total transmissions: " << stats.total_transmissions << std::endl;
        std::cout << "Average insertions: " << (double)stats.total_insertions / stats.total_transmissions << std::endl;
        std::cout << "Average deletions: " << (double)stats.total_deletions / stats.total_transmissions << std::endl;
        std::cout << "Max expansion: " << stats.max_expansion << std::endl;
        std::cout << "Max contraction: " << stats.max_contraction << std::endl;
    }

private:
    int Nb;
    int k;  // チャネルメモリ長
    int max_expansion;
    int max_contraction;
    
    // Julia通信路を呼び出す関数
    char* call_julia_channel(const char* input_seq, int input_len, int* output_len);
    struct Statistics {
        int total_transmissions = 0;
        int total_insertions = 0;
        int total_deletions = 0;
        int max_expansion = 0;
        int max_contraction = 0;
    } stats;

};

#endif