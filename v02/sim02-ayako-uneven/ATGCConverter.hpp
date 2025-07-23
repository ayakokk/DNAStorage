#ifndef ATGCCONVERTER_HPP
#define ATGCCONVERTER_HPP

// 単体変換関数
char binary_to_atgc(unsigned char b1, unsigned char b2);
void atgc_to_binary(char atgc, unsigned char &b1, unsigned char &b2);

// 配列変換関数
void binary_array_to_atgc(char *atgc_output, const unsigned char *binary_input, int binary_len);
int atgc_array_to_binary(unsigned char *binary_output, const char *atgc_input, int atgc_len);

#endif