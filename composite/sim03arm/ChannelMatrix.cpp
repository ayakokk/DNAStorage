#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <fstream>
#include <iostream>

#include "ChannelMatrix.hpp"

//================================================================================
void ChannelMatrix::PrintVect(const unsigned long *V, int len, const char *pre, const char *post){
  int w=20;
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
unsigned long ChannelMatrix::sum(const unsigned long *V, int len){
  unsigned long s=0;
  for(int i=0;i<len;i++) s+=V[i];
  return s;
}

//================================================================================
double ChannelMatrix::entropy(const unsigned long *V, int len){
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
ChannelMatrix::ChannelMatrix(int _M, int _N){
  M = _M;
  N = _N;
  printf("# ChannelMatrix: M=%d N=%d\n",M,N);
  assert(M>0 && N>0);
  //---
  cntX = new unsigned long [M];
  cntY = new unsigned long [N];
  cnt  = new unsigned long * [M];
  Pxy  = new double * [M];
  for(int i=0;i<M;i++){
    cnt[i] = new unsigned long [N];
    Pxy[i] = new double [N];
  } // for i
  clear();
}

//================================================================================
ChannelMatrix::ChannelMatrix(const char *fnPxy){
  double val;
  std::ifstream fin;
  fin.open(fnPxy, std::ios::in | std::ios::binary);
  if(!fin){
    std::cout << "Cannot open " << fnPxy << std::endl;
    exit(1);
  } // if
  fin.read((char *)&M, sizeof(int));
  fin.read((char *)&N, sizeof(int));
  //---
  cntX = new unsigned long [M];
  cntY = new unsigned long [N];
  cnt  = new unsigned long * [M];
  Pxy  = new double * [M];
  for(int i=0;i<M;i++){
    cnt[i] = new unsigned long [N];
    Pxy[i] = new double [N];
  } // for i
  clear();
  //---
  for(int i=0;i<M;i++){
    for(int j=0;j<N;j++){
      fin.read((char *)&val, sizeof(double));
      assert(val>=0.0 && val<=1.0);
      Pxy[i][j] = val;
    } // for j
  } // for i
  fin.close();
  printf("# ChannelMatrix: %s (%dx%d)\n",fnPxy,M,N);
}

//================================================================================
ChannelMatrix::~ChannelMatrix(){
  for(int i=0;i<M;i++){
    delete [] cnt[i];
    delete [] Pxy[i];
  } // for i
  delete [] cnt;
  delete [] Pxy;
  delete [] cntX;
  delete [] cntY;
  printf("# ChannelMatrix: deleted\n");
}

//================================================================================
void ChannelMatrix::clear(){
  for(int i=0;i<M;i++){
    for(int j=0;j<N;j++){
      cnt[i][j] = 0;
      Pxy[i][j] = 0.0;
    } // for j
    cntX[i] = 0;
  } // for i
  for(int j=0;j<N;j++) cntY[j] = 0;
  cntAll=0;
}

//================================================================================
void ChannelMatrix::countup(int x, int y){
  assert(x>=0 && x<M);
  assert(y>=0 && y<N);
  cnt[x][y]++;
  cntX[x]++;
  cntY[y]++;
  cntAll++;
  Pxy[x][y] = (double)cnt[x][y]/cntAll;
}

//================================================================================
void ChannelMatrix::PrintCnt(){
  printf("cntX:\n");
  PrintVect(cntX,M,"","");
  printf("cntY:\n");
  PrintVect(cntY,N,"","");
  for(int i=0;i<M;i++){
    printf("cnt[%02d]: ",i);
    for(int j=0;j<N;j++) printf("%06lu ",cnt[i][j]);
    printf("\n");
    //PrintVect(cnt[i],N,"","");
  } // for i
}

//================================================================================
void ChannelMatrix::PrintPxy(){
  for(int i=0;i<M;i++){
    printf("Pxy[%02d]: ",i);
    for(int j=0;j<N;j++) printf("%.2e ",Pxy[i][j]);
    printf("\n");
    //PrintVect(cnt[i],N,"","");
  } // for i
}

//================================================================================
void ChannelMatrix::WritePxy(const char *fn){
  double val;
  std::ofstream fout;
  fout.open(fn, std::ios::out | std::ios::binary | std::ios::trunc );
  if(!fout){
    std::cout << "Cannot open " << fn << std::endl;
    exit(1);
  } // if
  fout.write((char *)&M,sizeof(int));
  fout.write((char *)&N,sizeof(int));
  for(int i=0;i<M;i++){
    for(int j=0;j<N;j++){
      val = Pxy[i][j];
      fout.write((char *)&val,sizeof(double));
    } // for j
  } // for i
  fout.close();
  // ---- verify
  class ChannelMatrix *CMV = new class ChannelMatrix(fn);
  assert( ChannelMatrix_IsEqual(this,CMV) );
  delete CMV;
}
  
//================================================================================
double ChannelMatrix::Hx(){
  return entropy(cntX,M);
}

//================================================================================
double ChannelMatrix::Hxy(){
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
double ChannelMatrix::Ixy(){
  return Hx()-Hxy();
}

//================================================================================
void ChannelMatrix::GetPy(double *Py){
  unsigned long s = sum(cntY,N);
  assert(s>0);
  for(int j=0;j<N;j++) Py[j] = (double)cntY[j]/s;
}

//================================================================================
int ChannelMatrix::GetM(){return M;}
int ChannelMatrix::GetN(){return N;}

//================================================================================
double ChannelMatrix::GetPxy(int x, int y){
  assert( x>=0 && x<M );
  assert( y>=0 && y<N );
  return Pxy[x][y];
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
bool ChannelMatrix_IsEqual(class ChannelMatrix *CM0, class ChannelMatrix *CM1){
  int M = CM0->GetM();
  int N = CM0->GetN();
  if( CM1->GetM() != M ) return false;
  if( CM1->GetN() != N ) return false;
  for(int i=0;i<M;i++){
    for(int j=0;j<N;j++){
      if( CM0->GetPxy(i,j) != CM1->GetPxy(i,j) ) return false;
    } // for j
  } // for i
  return true;
}
