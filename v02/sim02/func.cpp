#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <vector>

#include "func.hpp"

#define BSIZE 8192

//================================================================================
int max(int a, int b){return (a>b)? a : b;}
int min(int a, int b){return (a<b)? a : b;}

//================================================================================
int argmax(const double *V, int len){
  assert(len>0);
  int ret=0;
  for(int i=1;i<len;i++){
    if(V[i]>V[ret]) ret=i;
  } // for i
  return ret;
}

//================================================================================
int argmin(const double *V, int len){
  assert(len>0);
  int ret=0;
  for(int i=1;i<len;i++){
    if(V[i]<V[ret]) ret=i;
  } // for i
  return ret;
}

//================================================================================
double max(const double *V, int len){return V[ argmax(V,len) ]; }
double min(const double *V, int len){return V[ argmin(V,len) ]; }

//================================================================================
void PrintVect(const double *V, int len, const char *pre, const char *post){
  int w=20;
  printf("%s",pre);
  for(int i=0;i<len;i++){
    printf("%02d:%.2e ",i,V[i]);
    if(i%w==w-1) printf("\n");
  } // for i
  printf("%s",post);
}

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
    if(i%w==0 && len>w) printf("%04d: ",i);
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
void FixedVect(int *V, int len, int Vmin, int Vmax, int *group_weights) {
  assert(len > 0);
  assert(Vmin < Vmax);
  assert(group_weights != NULL);
  
  int A = Vmax - Vmin + 1;
  
  // 重み付けを設定（値が大きいほど選ばれやすい）
  int *weights = new int[A];
  for (int i = 0; i < A; i++) {
    weights[i] = 1; // デフォルト重み
  }
  
  // グループ1: 最も高い確率（基準の2倍）
  int group1[] = {2, 8, 9, 15};
  for (int i = 0; i < sizeof(group1)/sizeof(group1[0]); i++) {
    int idx = group1[i] - Vmin;
    if (idx >= 0 && idx < A) {
      weights[idx] = group_weights[0];
    }
  }
  
  // グループ2: 2番目の確率（基準の1.5倍）
  int group2[] = {0, 1, 3, 6, 11, 14, 16, 17};
  for (int i = 0; i < sizeof(group2)/sizeof(group2[0]); i++) {
    int idx = group2[i] - Vmin;
    if (idx >= 0 && idx < A) {
      weights[idx] = group_weights[1];
    }
  }
  
  // グループ3: 3番目の確率（基準と同等）
  int group3[] = {5, 7, 10, 12};
  for (int i = 0; i < sizeof(group3)/sizeof(group3[0]); i++) {
    int idx = group3[i] - Vmin;
    if (idx >= 0 && idx < A) {
      weights[idx] = group_weights[2];
    }
  }
  
  // グループ4: 最も低い確率（基準の半分）
  int group4[] = {4, 13};
  for (int i = 0; i < sizeof(group4)/sizeof(group4[0]); i++) {
    int idx = group4[i] - Vmin;
    if (idx >= 0 && idx < A) {
      weights[idx] = group_weights[3];
    }
  }
  
  // 重み付きランダム選択
  int total_weight = 0;
  for (int i = 0; i < A; i++) {
    total_weight += weights[i];
  }
  
  for (int i = 0; i < len; i++) {
    int r = random() % total_weight;
    int cumulative = 0;
    for (int j = 0; j < A; j++) {
      cumulative += weights[j];
      if (r < cumulative) {
        V[i] = j + Vmin;
        break;
      }
    }
  }
  delete[] weights;
}
// ================================================================================
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

//================================================================================
void ReadConstraints(const char *fn, int *Rho, int *ell, int *Delta){
  FILE *fp;
  char *buf = new char [BSIZE];
  if((fp=fopen(fn,"r"))==NULL){
    fprintf(stderr,"Cannot open %s\n",fn);
    exit(1);
  } // if
  assert( fgets(buf,BSIZE,fp)!=NULL );
  (*Rho)   = atoi(strtok(buf, " \t"));
  (*ell)   = atoi(strtok(NULL," \t"));
  (*Delta) = atoi(strtok(NULL," \t\n"));
  fclose(fp);
  assert( (*Rho)  >0 );
  assert( (*ell)  >0 );
  assert( (*Delta)>0 );
}

//================================================================================
void HardDecision(int *V, const double **P, int N, int Q){
  assert( N>0 && Q>0 );
  for(int i=0;i<N;i++) V[i] = argmax( P[i], Q );
}

//================================================================================
int  HammingWeight(const unsigned char *V, int len){
  assert(len>0);
  int cnt=0;
  for(int i=0;i<len;i++){
    if(V[i]!=0) cnt++;
  } // for i
  return cnt;
}

//================================================================================
int  HammingWeight(const bool *V, int len){
  assert(len>0);
  int cnt=0;
  for(int i=0;i<len;i++){
    if(V[i]) cnt++;
  } // for i
  return cnt;
}

//================================================================================
int HammingDist(const int *V0, const int *V1, int len){
  assert(len>0);
  int cnt=0;
  for(int i=0;i<len;i++){
    if(V0[i]!=V1[i]) cnt++;
  } // for i
  return cnt;
}

//================================================================================
int RunLength(const unsigned char *V, int pos, int len){
  int L;
  int posL,posR;
  assert(len>0);
  assert(pos>=0 && pos<len);
  for(posL=pos-1;posL>=0;posL--){
    if(V[pos]!=V[posL]) break;
  } // for posL
  for(posR=pos+1;posR<len;posR++){
    if(V[pos]!=V[posR]) break;
  } // for posR
  L = posR-posL-1;
  assert(L>=0 && L<=len);
  return L;
}

//================================================================================
int MaxRunLength(const unsigned char *V, int len){
  int L,Lmax=1;
  int posL,posR;
  assert(len>0);
  for(posL=0;posL<len;posL++){
    for(posR=posL+1;posR<len;posR++){
      if(V[posR]!=V[posL]) break;
    } // for posR
    L = posR-posL;
    assert(L>=1 && L<=len);
    if(L>Lmax) Lmax=L;
    posL += (L-1); // jump
  } // for posL
  return Lmax;
}
