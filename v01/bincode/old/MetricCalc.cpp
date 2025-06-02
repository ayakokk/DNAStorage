#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <iostream>

#include "MetricCalc.hpp"


//============================================================================
int MetricCalc::min(int a, int b){ return (a<=b)? a : b; }

//============================================================================
void MetricCalc::PrintVect(unsigned char *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%u",V[i]);
  printf("%s",post);
}

//============================================================================
bool MetricCalc::CheckConv(){
  unsigned char *V = new unsigned char [beta];
  unsigned int u,u2;
  for(u=0;u<(unsigned int)beta2p;u++){
    ConvUintVect(V,u,beta);
    u2 = ConvVectUint(V,beta);
    // printf("%04X ",u);
    // PrintVect(V,beta,""," ");
    // printf("%04X\n",u2);
    if(u!=u2)return false;
  } // for u
  delete [] V;
  return true;
}

//============================================================================
void MetricCalc::CalcLval(int s, int t, const unsigned char *V0, const unsigned char *V1){
  assert(s>=0 && s<=beta);
  assert(t>=0 && t<=beta);
  assert(s>0 || t>0);  // (s,t)!=(0,0)
  int v;
  
  L[s][t]=INT_MAX;
  if(s>0){        // del -----
    v = L[s-1][t]+1;
    //if(t<beta) v++;
    L[s][t] = min(L[s][t],v);
  } // if s
  if(t>0){        // ins -----
    v = L[s][t-1]+1;
    //if(s<beta) v++;
    L[s][t] = min(L[s][t],v);
  } // if t
  if(s>0 && t>0){ // sub -----
    v = L[s-1][t-1];
    if(V0[s-1]!=V1[t-1]) v++;
    L[s][t] = min(L[s][t],v);
  } // if s,t>0
}

//============================================================================
void MetricCalc::SetMT(){
  unsigned int u0,u1;
  unsigned char *V0 = new unsigned char [beta];
  unsigned char *V1 = new unsigned char [beta];
  int t;
  
  //TMP
  u0 = 3;
  u1 = 13;

  for(u0=0;u0<(unsigned int)beta2p;u0++){
    for(u1=0;u1<(unsigned int)beta2p;u1++){
      ConvUintVect(V0,u0,beta);
      ConvUintVect(V1,u1,beta);
      L[0][0]=0;
      for(int i=1;i<=2*beta;i++){
	for(int s=i;s>=0;s--){
	  t = i-s;
	  if(s>=0 && s<=beta && t>=0 && t<=beta){
	    CalcLval(s,t,V0,V1);
	    //printf("%d %d\n",s,t);
	  } // if s,t
	} // for s
      } // for i
      MT[u0][u1] = L[beta][beta];
      
      //(dbg)
      // printf("%04X %04X: ",u0,u1);
      // PrintVect(V0,beta,""," ");
      // PrintVect(V1,beta,""," ");
      // printf("[%d]\n",MT[u0][u1]);
      // DumpL();
    } // for u1
  } // for u0
  
  delete [] V0;
  delete [] V1;
}

//============================================================================
void MetricCalc::DumpL(){
  printf("L:\n");
  for(int i=0;i<=beta;i++){
    for(int j=0;j<=beta;j++) printf("%02d ",L[i][j]);
    printf("\n");
  } // for i
}

//============================================================================
void MetricCalc::DumpMT(){
  for(int i=0;i<beta2p;i++){
    printf("%04X: ",i);
    for(int j=0;j<beta2p;j++) printf("%d ",MT[i][j]);
    printf("\n");
  } // for i
}

//============================================================================
//============================================================================
//============================================================================

//============================================================================
MetricCalc::MetricCalc(int _beta){
  beta  = _beta;
  beta2p= (int)pow(2,beta);
  printf("# MetricCalc: beta=%d beta2p=%d\n",beta,beta2p);
  assert(beta>=1 && beta<=15);
  L = new int * [beta+1];
  MT= new int * [beta2p];
  for(int i=0;i<beta+1;i++) L[i] = new int [beta+1];
  for(int i=0;i<beta2p;i++) MT[i]= new int [beta2p];
  
  assert(CheckConv());
  SetMT();
  //DumpMT();
  printf("# MetricCalc: generated MT\n");
}

//============================================================================
MetricCalc::~MetricCalc(){
  for(int i=0;i<beta+1;i++) delete [] L[i];
  for(int i=0;i<beta2p;i++) delete [] MT[i];
  delete [] L;
  delete [] MT;
  printf("# MetricCalc: deleted\n");
}

//============================================================================
int MetricCalc::calc(const unsigned char *V0, const unsigned char *V1){
  unsigned int u0,u1;
  u0 = ConvVectUint(V0,beta);
  u1 = ConvVectUint(V1,beta);
  assert(u0>=0 && u0<(unsigned int)beta2p);
  assert(u1>=0 && u1<(unsigned int)beta2p);
  return MT[u0][u1];
}

//============================================================================
void MetricCalc::ConvUintVect(unsigned char *V, unsigned int u, unsigned int len){
  assert(len<=(unsigned int)sizeof(unsigned int)*8);
  for(unsigned int i=0;i<len;i++){
    V[i] = 0x00000001 & u;
    u >>= 1;
  } // for i
}

//============================================================================
unsigned int MetricCalc::ConvVectUint(const unsigned char *V,   unsigned int len){
  assert(len<=(unsigned int)sizeof(unsigned int)*8);
  unsigned int u=0;
  for(int i=len-1;i>=0;i--){
    assert(V[i]==0 || V[i]==1);
    u<<=1;
    u += V[i];
  } // for i
  return u;
}

