// JuliaChannel.hpp
#ifndef JULIACHANNEL_HPP
#define JULIACHANNEL_HPP

#include "ATGCConverter.hpp"
#include <julia.h>

class JuliaChannel {
public:
    JuliaChannel(int k = 6);
    ~JuliaChannel();
    
    int transmit(unsigned char *RW, const unsigned char *CW, int Nb);
    
    // IDSchannelとの互換性のためのメソッド
    int GetN() const { return Nb; }
    int GetDmax() const { return max_expansion; }
    int GetDmin() const { return -max_contraction; }

private:
    int Nb;
    int k;
    int max_expansion;
    int max_contraction;
    
    // Julia関数への参照
    jl_function_t *channel_func;
    jl_module_t *main_module;
    
    // Julia初期化フラグ（静的メンバ）
    static bool julia_initialized;
    static int julia_ref_count;
    
    // Julia環境の初期化
    void init_julia();
    void cleanup_julia();
    
    // Juliaチャネル関数を呼び出す
    char* call_julia_channel(const char* input_seq, int input_len, int* output_len);
};

#endif

// JuliaChannel.cpp
#include "JuliaChannel.hpp"
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <sstream>

// 静的メンバの初期化
bool JuliaChannel::julia_initialized = false;
int JuliaChannel::julia_ref_count = 0;

JuliaChannel::JuliaChannel(int _k) : k(_k), channel_func(nullptr), main_module(nullptr) {
    max_expansion = 30;
    max_contraction = 100;
    
    init_julia();
    
    std::cout << "# JuliaChannel: k=" << k 
              << " expected range: [" << GetDmin() << ", " << GetDmax() << "]" << std::endl;
}

JuliaChannel::~JuliaChannel() {
    cleanup_julia();
    std::cout << "# JuliaChannel: destroyed" << std::endl;
}

void JuliaChannel::init_julia() {
    julia_ref_count++;
    
    if (!julia_initialized) {
        // Julia初期化
        jl_init();
        julia_initialized = true;
        
        std::cout << "# Julia initialized" << std::endl;
    }
    
    // Juliaコードをロード
    std::string julia_code = R"(
        # パスを追加
        push!(LOAD_PATH, "simulator")
        
        # 必要なモジュールをロード
        using DelimitedFiles
        using Random
        
        # グローバル変数とロード処理
        global k_global = 0
        global loaded = false
        
        # チャネル関数の定義
        function setup_channel(k_val::Int)
            global k_global = k_val
            global loaded = false
            
            # ここで必要なファイルをinclude
            include("simulator/loadProb.jl")
            include("simulator/functions.jl")
            include("simulator/channel.jl")
            
            global loaded = true
            return true
        end
        
        # シンプルなチャネル関数
        function simulate_channel(seq::String, k_val::Int)
            if !loaded || k_global != k_val
                setup_channel(k_val)
            end
            
            # channel関数を呼び出し
            result = channel(1, k_val, [">sequence", seq])
            
            # 結果を文字列に変換
            if length(result) > 0 && length(result[1]) > 0
                return join(result[1])
            else
                return ""
            end
        end
    )";
    
    // Juliaコードを評価
    jl_eval_string(julia_code.c_str());
    if (jl_exception_occurred()) {
        std::cerr << "Julia error during initialization: "
                  << jl_typeof_str(jl_exception_occurred()) << std::endl;
        // Base.stdout を取得して showerror に渡す
        jl_value_t* stdout_obj = jl_get_global(jl_base_module, jl_symbol("stdout"));
        jl_call2(jl_get_function(jl_base_module, "showerror"),
                 stdout_obj, jl_exception_occurred());
        jl_printf(jl_stdout_stream(), "\n");
    }
    
    // モジュールと関数を取得
    main_module = jl_main_module;
    channel_func = jl_get_function(main_module, "simulate_channel");
    
    if (!channel_func) {
        std::cerr << "Error: Could not find simulate_channel function" << std::endl;
    }
}

void JuliaChannel::cleanup_julia() {
    julia_ref_count--;
    
    if (julia_ref_count == 0 && julia_initialized) {
        jl_atexit_hook(0);
        julia_initialized = false;
        std::cout << "# Julia cleaned up" << std::endl;
    }
}

char* JuliaChannel::call_julia_channel(const char* input_seq, int input_len, int* output_len) {
    if (!channel_func) {
        std::cerr << "Error: Julia channel function not initialized" << std::endl;
        return nullptr;
    }
    
    // Julia文字列を作成
    jl_value_t *seq_str = jl_cstr_to_string(input_seq);
    jl_value_t *k_val = jl_box_int64(k);
    
    // GCルート保護
    JL_GC_PUSH2(&seq_str, &k_val);
    
    // 関数を呼び出し
    jl_value_t *result = jl_call2(channel_func, seq_str, k_val);
    
    if (jl_exception_occurred()) {
        std::cerr << "Julia error during channel simulation: " 
                  << jl_typeof_str(jl_exception_occurred()) << std::endl;
        // 修正: Base.stdout を取得して showerror に渡す
        jl_value_t* stdout_obj = jl_get_global(jl_base_module, jl_symbol("stdout"));
        jl_call2(jl_get_function(jl_base_module, "showerror"), 
                 stdout_obj, jl_exception_occurred());
        jl_printf(jl_stdout_stream(), "\n");
        JL_GC_POP();
        return nullptr;
    }
    
    // 結果を取得
    char* output = nullptr;
    if (jl_is_string(result)) {
        const char* julia_str = jl_string_ptr(result);
        *output_len = strlen(julia_str);
        output = new char[*output_len + 1];
        strcpy(output, julia_str);
    } else {
        std::cerr << "Error: Julia function did not return a string" << std::endl;
        *output_len = 0;
    }
    
    JL_GC_POP();
    
    return output;
}

int JuliaChannel::transmit(unsigned char *RW, const unsigned char *CW, int Nb) {
    assert(Nb % 2 == 0);
    
    // 1. バイナリ→ATGC変換
    int atgc_len = Nb / 2;
    char *atgc_input = new char[atgc_len + 1];
    binary_array_to_atgc(atgc_input, CW, Nb);
    
    std::cout << "# Julia Input: length=" << atgc_len << " bases" << std::endl;
    
    // 2. Julia通信路を呼び出し
    int output_atgc_len;
    char* atgc_output = call_julia_channel(atgc_input, atgc_len, &output_atgc_len);

    
    if (atgc_output == nullptr) {
        std::cerr << "Error: Julia channel simulation failed" << std::endl;
        delete[] atgc_input;
        return -1;
    }
    
    std::cout << "# Julia Output: length=" << output_atgc_len << " bases" << std::endl;
    
    // 3. ATGC→バイナリ逆変換
    int output_binary_len = atgc_array_to_binary(RW, atgc_output, output_atgc_len);
    
    std::cout << "# Binary conversion: " << Nb << " bits -> " << output_binary_len << " bits" << std::endl;
    std::cout << "#   Change: " << (output_binary_len - Nb) << " bits" << std::endl;
    
    // メモリ解放
    delete[] atgc_input;
    delete[] atgc_output;
    
    return output_binary_len;
}