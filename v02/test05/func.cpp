#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "func.hpp"

//================================================================================
void PrintVectX(const int *V, int len, const char *pre, const char *post){
  int w=20;
  printf("%s",pre);
  for(int i=0;i<len;i++){
    if(i%w==0) printf("%04d: ",i);
    printf("%04X ",V[i]);
    if(i%w==w-1) printf("\n");
  } // for i
  printf("%s",post);
}

//================================================================================
void PrintVectB(const unsigned char *V, int len, const char *pre, const char *post){
  int w=100;
  printf("%s",pre);
  for(int i=0;i<len;i++){
    if(i%w==0) printf("%04d: ",i);
    printf("%u",V[i]);
    if(i%w==w-1) printf("\n");
  } // for i
  printf("%s",post);
}

//================================================================================
void RandVect(int *V, int len, int Vmin, int Vmax){
  assert(len>0);
  assert(Vmin<Vmax);
  int A = Vmax-Vmin+1;
  for(int i=0;i<len;i++){
    V[i] = random()%A;
    V[i] += Vmin;
  } // for i
}

//================================================================================
long VectToLong(const unsigned char *V, int len){
  assert(len>0);
  long val = 0;
  for(int i=0;i<len;i++){
    assert(V[i]==0 || V[i]==1);
    val <<= 1;
    if(V[i]==1) val |= 0x1;
  } // for i
  return val;
}

//================================================================================
void LongToVect(unsigned char *V, long val, int len){
  assert(len>0 && val>=0);
  long mask = 0x1 << (len-1);
  for(int i=0;i<len;i++){
    V[i] = ( (val & mask)==0 )? 0 : 1;
    mask >>= 1;
  } // for i
}
