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
int MetricCalc::HammingDist(unsigned char v0, unsigned char v1){ return (v0==v1)? 0 : 1; }

//============================================================================
int MetricCalc::HammingDist(const unsigned char *V0, const unsigned char *V1, int len){
  int d=0;
  for(int i=0;i<len;i++)
    if(V0[i]!=V1[i]) d++;
  return d;
}

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
  assert(s>=1 && s<=beta);
  assert(t>=0 && t<=beta);
  assert(s>0 || t>0);  // (s,t)!=(0,0)
  //int v;
  
  L[s][t] = L[s-1][t]*Qd[0];  // del
  if(t>=1) L[s][t]+= L[s-1][t-1]*Qs[ HammingDist(V0[s-1],V1[t-1]) ]; // sub
  if(t>=2) L[s][t]+= L[s-1][t-2]*Qi[ HammingDist(V0[s-1],V1[t-2])+HammingDist(V0[s-1],V1[t-1]) ]; //ins
  
  /*
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
  */
}

//============================================================================
void MetricCalc::SetD(){
  // ins
  Qi[0] = Pi*(1.0-Ps)*(1.0-Ps);
  Qi[1] = Pi*Ps      *(1.0-Ps);
  Qi[2] = Pi*Ps      *Ps      ;
  // del
  Qd[0] = Pd;
  // sub
  Qs[0] = Pt*(1.0-Ps);
  Qs[1] = Pt*Ps;
}

//============================================================================
void MetricCalc::SetMT(){
  unsigned int u0,u1;
  unsigned char *V0 = new unsigned char [beta];
  unsigned char *V1 = new unsigned char [beta];
  
  //TMP
  //u0 = 3;
  //u1 = 23;

  for(u0=0;u0<(unsigned int)beta2p;u0++){
    for(u1=0;u1<(unsigned int)beta2p;u1++){
      ConvUintVect(V0,u0,beta);
      ConvUintVect(V1,u1,beta);
      L[0][0]=1.0;
      for(int j=1;j<=beta;j++) L[0][j]=0.0;
      for(int i=1;i<=beta;i++){
	for(int j=0;j<=beta;j++) CalcLval(i,j,V0,V1);
      } // for i
      MT[u0][u1] = L[beta][beta];
      
      //(dbg)
      // printf("%04X %04X: ",u0,u1);
      // PrintVect(V0,beta,""," ");
      // PrintVect(V1,beta,""," ");
      // printf("[%e]\n",MT[u0][u1]);
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
    for(int j=0;j<=beta;j++) printf("%.2e ",L[i][j]);
    printf("\n");
  } // for i
}

//============================================================================
void MetricCalc::DumpMT(){
  for(int i=0;i<beta2p;i++){
    printf("%04X: ",i);
    for(int j=0;j<beta2p;j++) printf("%.2e ",MT[i][j]);
    printf("\n");
  } // for i
}

//============================================================================
//============================================================================
//============================================================================

//============================================================================
MetricCalc::MetricCalc(int _beta, double _Pi, double _Pd, double _Ps){
  beta  = _beta;
  beta2p= (int)pow(2,beta);
  Pi    = _Pi;
  Pd    = _Pd;
  Ps    = _Ps;
  Pt    = 1-Pi-Pd;
  
  printf("# MetricCalc: beta=%d beta2p=%d\n",beta,beta2p);
  printf("# MetricCalc: (Pi,Pd,Ps;Pt)=(%e,%e,%e;%e)\n",Pi,Pd,Ps,Pt);
  assert(beta>=1 && beta<=15);
  assert(Pi>=0 && Pi<0.5);
  assert(Pd>=0 && Pd<0.5);
  assert(Ps>=0 && Ps<0.5);
  assert(Pt>=0 && Pt<=1.0);

  Qi = new double [3];  // 0,1,2
  Qd = new double [1];  // 0
  Qs = new double [2];  // 0,1
  
  L = new double * [beta+1];
  MT= new double * [beta2p];
  for(int i=0;i<beta+1;i++) L[i] = new double [beta+1];
  for(int i=0;i<beta2p;i++) MT[i]= new double [beta2p];
  
  assert(CheckConv());
  SetD();
  printf("# MetricCalc: Qi=(%.2e,%.2e,%.2e) Qd=(%.2e) Qs=(%.2e,%.2e)\n",Qi[0],Qi[1],Qi[2],Qd[0],Qs[0],Qs[1]);

  SetMT();
  //DumpMT();
  printf("# MetricCalc: generated MT\n");
}

//============================================================================
MetricCalc::~MetricCalc(){
  delete [] Qi;
  delete [] Qd;
  delete [] Qs;
  for(int i=0;i<beta+1;i++) delete [] L[i];
  for(int i=0;i<beta2p;i++) delete [] MT[i];
  delete [] L;
  delete [] MT;
  printf("# MetricCalc: deleted\n");
}

//============================================================================
double MetricCalc::calc(const unsigned char *V0, const unsigned char *V1){
  double p=0.0;
  unsigned int u0,u1;
  unsigned char **Vx = new unsigned char * [4];  // 0:L0 1:L1 2:R0 3:R1
  for(int i=0;i<4;i++) Vx[i] = new unsigned char [beta];

  for(int i=0;i<beta-1;i++){
    Vx[0][i]=V1[i+1];
    Vx[1][i]=V1[i+1];
  } // for i
  for(int i=1;i<beta;i++){
    Vx[2][i]=V1[i-1];
    Vx[3][i]=V1[i-1];
  } // for i
  Vx[0][beta-1]=0;
  Vx[1][beta-1]=1;
  Vx[2][0]     =0;
  Vx[3][0]     =1;

  u0 = ConvVectUint(V0,beta);
  assert(u0>=0 && u0<(unsigned int)beta2p);

  for(int i=0;i<4;i++){
    u1 = ConvVectUint(Vx[i],beta);
    assert(u1>=0 && u1<(unsigned int)beta2p);
    if(MT[u0][u1]>p) p=MT[u0][u1];
    //printf("  %02d %02d: %e\n",u0,u1,MT[u0][u1]);
  } // for i
  //printf("---\n");
  
  for(int i=0;i<4;i++) delete [] Vx[i];
  delete [] Vx;
  return p;
}

//============================================================================
double MetricCalc::calc1(const unsigned char *V0, const unsigned char *V1){
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

