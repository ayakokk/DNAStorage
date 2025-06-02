#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"


//================================================================================
void PrintVect(int *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d ",V[i]);
  printf("%s",post);
}

//================================================================================
void PrintVect2(unsigned char *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%u",V[i]);
  printf("%s",post);
}

//================================================================================
void PrintVect2(bool *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d",V[i]);
  printf("%s",post);
}

//=================================================================================
void PrintSP(int len){
  for(int i=0;i<len;i++) printf(" ");
}

//=================================================================================
void SetVal(bool *V, bool val, int len){
  for(int i=0;i<len;i++) V[i]=val;
}

//=================================================================================
int MaxRunLength(const unsigned char *V, int len){
  int RL=1,RLmax=0;
  for(int i=1;i<len;i++){
    if(V[i]==V[i-1]){
      RL++;
      if(RL>RLmax) RLmax=RL;
    } else {
      RL=1;
    } // if
  } // for i
  return RLmax;
}

//=================================================================================
int HammingWeight(const unsigned char *V, int len){
  int w=0;
  for(int i=0;i<len;i++)
    if(V[i]!=0) w++;
  return w;
}

//=================================================================================
int HammingWeight(const bool *V, int len){
  int w=0;
  for(int i=0;i<len;i++)
    if(V[i]) w++;
  return w;
}

//=================================================================================
int HammingDist(const unsigned char *V0, const unsigned char *V1, int len){
  int d=0;
  for(int i=0;i<len;i++)
    if(V0[i]!=V1[i]) d++;
  return d;
}

//=================================================================================
int HammingDist(const bool *V0, const bool *V1, int len){
  int d=0;
  for(int i=0;i<len;i++)
    if(V0[i]!=V1[i]) d++;
  return d;
}

//=================================================================================
//=================================================================================
//=================================================================================

//=================================================================================
void ConvIntV4(unsigned char *V, unsigned int u, int len){
  assert(len>=0 && len<=(int)sizeof(unsigned int)*4);
  for(int i=0;i<len;i++){
    V[i] = u%4;
    u /= 4;
  } // for i
}

//=================================================================================
unsigned int ConvV4Int(const unsigned char *V, int len){
  assert(len>=0 && len<=(int)sizeof(unsigned int)*4);
  unsigned int u=0;
  for(int i=len-1;i>=0;i--){
    assert(V[i]<4);
    u *= 4;
    u += V[i];
  } // for i
  return u;
}

//=================================================================================
void ConvIntV2(unsigned char *V, unsigned int x, int len){
  assert(len>=0 && x<pow(2,len));
  unsigned int x2=x;
  for(int i=0;i<len;i++){
    V[i] = x%2;
    x>>=1;
  } // for i
  assert(ConvV2Int(V,len)==x2);
}

//=================================================================================
unsigned int ConvV2Int(const unsigned char *V, int len){
  assert(len>=0);
  unsigned int x=0;
  for(int i=len-1;i>=0;i--){
    x<<=1;
    if(V[i]==1) x|=0x00000001;
  } // for i
  return x;
}

//=================================================================================
//=================================================================================
//=================================================================================

//=================================================================================
int max(int v0, int v1){return (v0>=v1)? v0:v1;}
int min(int v0, int v1){return (v0<=v1)? v0:v1;}
