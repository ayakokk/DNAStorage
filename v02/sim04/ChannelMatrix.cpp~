#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "MutualInfo.hpp"

//================================================================================
void MutualInfo::PrintVect(const unsigned long *V, int len, const char *pre, const char *post){
  int w=16;
  printf("%s",pre);
  for(int i=0;i<len;i++){
    if(i%w==0) printf("%04d: ",i);
    printf("%06lu ",V[i]);
    if(i%w==w-1) printf("\n");
  } // for i
  if(len%w!=0) printf("\n");
  printf("%s",post);
}

//================================================================================
unsigned long MutualInfo::sum(const unsigned long *V, int len){
  unsigned long s=0;
  for(int i=0;i<len;i++) s+=V[i];
  return s;
}

//================================================================================
double MutualInfo::entropy(const unsigned long *V, int len){
  unsigned long s = sum(V,len);
  double p, h=0;
  assert(s>0);
  for(int i=0;i<len;i++){
    p = (double) V[i]/s;
    if(p>0) h -= p*log2(p);
  } // for i
  return h;
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
MutualInfo::MutualInfo(int _M, int _N){
  M = _M;
  N = _N;
  printf("# Mutual Info: M=%d N=%d\n",M,N);
  assert(M>0 && N>0);
  //---
  cntX = new unsigned long [M];
  cntY = new unsigned long [N];
  cnt  = new unsigned long * [M];
  for(int i=0;i<M;i++) cnt[i] = new unsigned long [N];
  clear();
}

//================================================================================
MutualInfo::~MutualInfo(){
  for(int i=0;i<M;i++) delete [] cnt[i];
  delete [] cnt;
  delete [] cntX;
  delete [] cntY;
  printf("# Mutual Info: deleted\n");
}

//================================================================================
void MutualInfo::clear(){
  for(int i=0;i<M;i++){
    for(int j=0;j<N;j++) cnt[i][j]=0;
    cntX[i] = 0;
  } // for i
  for(int j=0;j<N;j++) cntY[j] = 0;
}

//================================================================================
void MutualInfo::countup(int x, int y){
  assert(x>=0 && x<M);
  assert(y>=0 && y<N);
  cnt[x][y]++;
  cntX[x]++;
  cntY[y]++;
}

//================================================================================
void MutualInfo::PrintCnt(){
  printf("cntX:\n");
  PrintVect(cntX,M,"","");
  printf("cntY:\n");
  PrintVect(cntY,N,"","");
  for(int i=0;i<M;i++){
    printf("cnt[%d]\n",i);
    PrintVect(cnt[i],N,"","");
  } // for i
}

//================================================================================
double MutualInfo::Hx(){
  return entropy(cntX,M);
}

//================================================================================
double MutualInfo::Hxy(){
  double h, hxy = 0.0;
  unsigned long s = sum(cntY,N);
  unsigned long *V = new unsigned long [M];
  for(int j=0;j<N;j++){
    if(cntY[j]>0){
      for(int i=0;i<M;i++) V[i] = cnt[i][j];
      h = entropy(V,M);
      hxy += ((double)cntY[j]/s)*h;
    } // for j
  } // for j
  delete [] V;
  return hxy;
}

//================================================================================
double MutualInfo::Ixy(){
  return Hx()-Hxy();
}

//================================================================================
void MutualInfo::GetPy(double *Py){
  unsigned long s = sum(cntY,N);
  assert(s>0);
  for(int j=0;j<N;j++) Py[j] = (double)cntY[j]/s;
}
