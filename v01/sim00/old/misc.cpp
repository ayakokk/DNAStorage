#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "misc.hpp"


//=================================================================================
void PrintVect(const unsigned long *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%lu,",V[i]);
  printf("%s",post);
}

//=================================================================================
void PrintVectRatio(const unsigned long *V, int len){
  unsigned long sum=0;
  for(int i=0;i<len;i++) sum+=V[i];
  assert(sum>0);
  for(int i=0;i<len;i++) printf("%d:%.2e ",i,(double)V[i]/sum);
}

//=================================================================================
void PrintArray(const double **V, int l0, int l1){
  int w=16;
  for(int i=0;i<l0;i++){
    printf("[%04d]\n",i);
    for(int j=0;j<l1;j++){
      printf("%02X:%.2e ",j,V[i][j]);
      if(j%w==w-1) printf("\n");
    } // for j
    if(l1%w!=0) printf("\n");
  } // for i
}

//=================================================================================
void GenRndVect(int *V, int mod, int len){
  assert(mod>=2);
  for(int i=0;i<len;i++) V[i] = random()%mod;
}

//=================================================================================
double ErrEntropy(const int *X, const double **P, int len, int Q){
  int    x;
  double p,e=0;
  for(int i=0;i<len;i++){
    x = X[i];
    assert(x>=0 && x<Q);
    p = P[i][x];
    assert(p>0);
    e += log2(p);
  } // for i
  return -e/(double)len;
}

//=================================================================================
/*
int HammingDist(const int *V0, const int *V1, int len){
  int d=0;
  for(int i=0;i<len;i++) if(V0[i]!=V1[i]) d++;
  return d;
}
*/
