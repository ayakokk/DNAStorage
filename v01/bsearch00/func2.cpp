#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "func.hpp"

//============================================================================
void PrintVect(const unsigned char *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d ",V[i]);
  printf("%s",post);
}

//============================================================================
void PrintVect(const int *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d ",V[i]);
  printf("%s",post);
}

//============================================================================
int max(const int *V, int len){
  int x=V[0];
  for(int i=1;i<len;i++)
    if(V[i]>x) x=V[i];
  return x;
}

//============================================================================
int binom(int n, int k){
  assert(n>=k);
  if(n==k || k==0) return 1;
  if(2*k>n) k=n-k; 
  return binom(n-1,k) + binom(n-1,k-1);
}

//============================================================================
void ConvUintVect(unsigned char *V, unsigned int u, unsigned int len){
  assert(len<=(unsigned int)sizeof(unsigned int)*8);
  for(unsigned int i=0;i<len;i++){
    V[i] = 0x00000001 & u;
    u >>= 1;
  } // for i
}

//============================================================================
unsigned int ConvVectUint(const unsigned char *V,   unsigned int len){
  assert(len<=(unsigned int)sizeof(unsigned int)*8);
  unsigned int u=0;
  for(int i=len-1;i>=0;i--){
    assert(V[i]==0 || V[i]==1);
    u<<=1;
    u += V[i];
  } // for i
  return u;
}
