#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "bmatrix.hpp"
#include "func.hpp"

//=================================================================================
void ClearVect(unsigned long *V, int len){
  for(int i=0;i<len;i++) V[i]=0;
}

//=================================================================================
void PrintVect(const unsigned char *V, int len){
  for(int i=0;i<len;i++) printf("%u ",V[i]);
}

//=================================================================================
void PrintVect(const unsigned long *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%lu,",V[i]);
  printf("%s",post);
}

//=================================================================================
void PrintVect(const int *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d,",V[i]);
  printf("%s",post);
}

//=================================================================================
void PrintVect4(const unsigned char *V, int len){
  char c;
  for(int i=0;i<len;i++){
    switch(V[i]){
    case 0: c='A'; break;
    case 1: c='T'; break;
    case 2: c='G'; break;
    case 3: c='C'; break;
    default: assert(false);
    }
    printf("%c",c);
  } // for i
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
void PrintArray(const int **V, int l0, int l1){
  int w=32;
  for(int i=0;i<l0;i++){
    printf("[%04d] ",i);
    for(int j=0;j<l1;j++){
      printf("%d ",V[i][j]);
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
void ConvIntVect(int u, unsigned char *V, int len){
  assert(len>=0 && len<=(int)sizeof(unsigned int)*4);
  for(int i=0;i<len;i++){
    V[i] = u%4;
    u /= 4;
  } // for i
}

//=================================================================================
int ConvVectInt(const unsigned char *V, int len){
  assert(len>=0 && len<=(int)sizeof(unsigned int)*4);
  int u=0;
  for(int i=len-1;i>=0;i--){
    assert(V[i]<4);
    u *= 4;
    u += V[i];
  } // for i
  return u;
}

//=================================================================================
void ConvToBinVect(unsigned char *V, unsigned int x, int len){
  assert(len>=0 && x<pow(2,len));
  unsigned int x2=x;
  for(int i=0;i<len;i++){
    V[i] = x%2;
    x>>=1;
  } // for i
  assert(ConvFromBinVect(V,len)==x2);
}

//=================================================================================
unsigned int ConvFromBinVect(const unsigned char *V, int len){
  assert(len>=0);
  unsigned int x=0;
  for(int i=len-1;i>=0;i--){
    x<<=1;
    if(V[i]==1) x|=0x00000001;
  } // for i
  // check
  /*
  unsigned char *V2 = new unsigned char [len];
  ConvToBinVect(V2,x,len);
  assert(HammingDist(V,V2,len)==0);
  delete [] V2;
  */
  return x;
}

//=================================================================================
int SymbolCnt(const unsigned char *V, unsigned char val, int len){
  int c=0;
  for(int i=0;i<len;i++)
    if(V[i]==val) c++;
  return c;
}

//=================================================================================
int GCbalance(const unsigned char *V, int len){
  for(int i=0;i<len;i++) assert(V[i]>=0 && V[i]<=3);
  return
    SymbolCnt(V,2,len)+ // G
    SymbolCnt(V,3,len)- // C
    SymbolCnt(V,0,len)- // A
    SymbolCnt(V,1,len); // T
}

//=================================================================================
double GCbalanceWin(const unsigned char *V, int len, int ws){
  int sum=0;
  assert(ws>=2 && ws<=len);
  for(int i=0;i<=len-ws;i++) sum += abs(GCbalance(&V[i],ws));
  return (double)sum/(len-ws+1);
}

//=================================================================================
void GCbalanceCnt(const unsigned char *V, int len, int ws, unsigned long *BLC){
  int b;
  for(int i=0;i<=len-ws;i++){
    b = GCbalance(&V[i],ws);
    assert(b>=-ws && b<=ws);
    BLC[b+ws]++;
  } // for i
}

//=================================================================================
int MaxRunLength(const unsigned char *V, int len){
  int L=1,Lmax=1;
  for(int i=1;i<len;i++){
    if(V[i]==V[i-1]){
      L++;
      if(L>Lmax)Lmax=L;
    } else {
      L=1;
    } // if
  } // for i
  return Lmax;
}

//=================================================================================
double AveRunLength(const unsigned char *V, int len){
  int L=1,Lsum=0,Lcnt=1;
  for(int i=1;i<len;i++){
    if(V[i]==V[i-1]){
      L++;
    } else {
      Lsum+=L;
      Lcnt++;
      L=1;
    } // if
  } // for i
  return (double)Lsum/Lcnt;
}

//=================================================================================
void CntRunLength(const unsigned char *V, int len, unsigned long *RLC, int RLmax){
  int L=1;
  for(int i=1;i<len;i++){
    if(V[i]==V[i-1]){
      L++;
    } else {
      RLC[min(L,RLmax)]++;
      L=1;
    } // if
  } // for i
}

//=================================================================================
void RandomSel(class bmatrix *CW, int num, class bmatrix *F){
  int M = CW->getM();
  int W = F->HammingWeight();
  assert(CW->getN()==1 && F->getN()==1);
  assert(F->getM()==M);
  assert(M-W>=num);

  // vector of valid indices
  int *V = new int [M-W];
  int pos=0;
  for(int i=0;i<M;i++){
    if(F->getV(i,0)==0){
      V[pos]=i;
      pos++;
    } // if
  } // for i
  assert(pos==M-W);

  // randomize
  for(int i=0;i<M-W;i++){
    pos = random()%(M-W);
    Swap(&V[i],&V[pos]);
  }  // for i

  // set CW
  CW->clear();
  for(int i=0;i<num;i++) CW->setV(V[i],0,1);

  // check
  for(int i=0;i<M;i++)
    if(F->getV(i,0)==1) assert(CW->getV(i,0)==0);
  assert(CW->HammingWeight()==num);

  delete [] V;
}

//=================================================================================
void Swap(int *a, int *b){
  int x;
  x=*a;
  *a=*b;
  *b=x;
}

//=================================================================================
int min(int a, int b){
  return (a<=b)? a : b;
}

//=================================================================================
int max(int a, int b){
  return (a>=b)? a : b;
}

//=================================================================================
int min(const int *V, int len){
  int x = V[0];
  for(int i=1;i<len;i++)
    if(V[i]<x) x=V[i];
  return x;
}

//=================================================================================
int max(const int *V, int len){
  int x = V[0];
  for(int i=1;i<len;i++)
    if(V[i]>x) x=V[i];
  return x;
}

//=================================================================================
int HammingDist(const int *V0, const int *V1, int len){
  int d=0;
  for(int i=0;i<len;i++)
    if(V0[i]!=V1[i]) d++;
  return d;
}

//=================================================================================
int HammingDist(const unsigned char *V0, const unsigned char *V1, int len){
  int d=0;
  for(int i=0;i<len;i++)
    if(V0[i]!=V1[i]) d++;
  return d;
}

//=================================================================================
double Normalize(double *V, int len){
  double s=0.0;
  for(int i=0;i<len;i++){
    assert(V[i]>=0.0);
    s+=V[i];
  } // for i
  assert(s>0.0);
  for(int i=0;i<len;i++) V[i]/=s;
  return s;
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
