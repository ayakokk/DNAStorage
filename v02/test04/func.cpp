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
