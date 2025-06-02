#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "CIDSchannel.hpp"

#define Drange 2.0 // Dmin=-Drange*N*Pd Dmax=Drange*N*Pi  

//================================================================================
int CIDSchannel::max(int a, int b){return (a>b)? a : b; }
int CIDSchannel::min(int a, int b){return (a<b)? a : b; }  

//================================================================================
int CIDSchannel::sum(const int *V, int len){
  int s=0;
  for(int i=0;i<len;i++) s+=V[i];
  return s;
}

//================================================================================
int CIDSchannel::HammingDist(const unsigned char *V0, const unsigned char *V1, int len){
  int d=0;
  for(int i=0; i<len; i++){
    if(V0[i] != V1[i]) d++;
  } // for i
  return d;
}

//================================================================================
unsigned char CIDSchannel::inv(unsigned char a){
  assert(a==0 || a==1);
  return (a==0)? 1 : 0 ;
}

//================================================================================
unsigned char CIDSchannel::CompSel(const int *dist){
  int    v = 0;
  double p = ((double)random()/RAND_MAX) * (double)Kc;
  double pp= dist[0];
  assert( sum(dist,4) == Kc );
  while(pp<p){
    v++;
    assert(v<=3);
    pp += dist[v];
  } // while
  assert(v>=0 && v<=3);
  return (unsigned char)v;
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
CIDSchannel::CIDSchannel(int _N, double _Pi, double _Pd, double _Ps, int _Dmin, int _Dmax, int _Kc, int _Nseq){
  N  = _N;
  Pi = _Pi;
  Pd = _Pd;
  Ps = _Ps;
  Kc = _Kc;
  Nseq = _Nseq;
  Dmin = _Dmin;
  Dmax = _Dmax;
  // Dmin = (int)floor(-(double)Drange*N*Pd);
  // Dmax = (int)ceil ( (double)Drange*N*Pi);
  printf("# CIDSchannel: N=%d (Pi,Pd,Ps)=(%e,%e,%e) (Dmin,Dmax)=(%d,%d)\n",N,Pi,Pd,Ps,Dmin,Dmax);
  printf("# CIDSchannel: Kc=%d Nseq=%d\n",Kc,Nseq);
  assert( N>0 );
  assert( Pi>=0.0 && Pi<0.5 );
  assert( Pd>=0.0 && Pd<0.5 );
  assert( Ps>=0.0 && Ps<0.5 );
  assert( Kc>0 && Nseq>0);
  XV = new unsigned char * [Nseq];
  DV = new int * [Nseq];
  for(int i=0;i<Nseq;i++){
    XV[i] = new unsigned char [N];
    DV[i] = new int [N+1];
  } // for i
}

//================================================================================
CIDSchannel::~CIDSchannel(){
  for(int i=0;i<Nseq;i++){
    delete [] XV[i];
    delete [] DV[i];
  } // for i
  delete [] XV;
  delete [] DV;
  printf("# CIDSchannel: deleted\n");
}

//================================================================================
void CIDSchannel::transmit(unsigned char **YV, int *N2, const int **XC){
  //----- Set XV from XC
  for(int i=0;i<Nseq;i++){
    for(int j=0;j<N;j++){
      XV[i][j] = CompSel( XC[j] );
    } // for j
  } // for i
  //----- transmit XV->YV
  for(int i=0;i<Nseq;i++){
    N2[i] = transmitOne(YV[i],DV[i],XV[i]);
  } // for i
}

//================================================================================
int CIDSchannel::transmitOne(unsigned char *Y, int *DD, const unsigned char *X){
  double p;
  int v;
  //----- set drift vector
  DD[0] = 0;
  for(int i=0;i<N;i++){
    p = (double) random()/RAND_MAX;
    if(DD[i]==Dmin){
      //--- Trans or Ins
      if(p<Pi) DD[i+1] = DD[i]+1; // ins
      else     DD[i+1] = DD[i];   // trans
    } else if(DD[i]==Dmax){
      //--- Trans or Del
      if(p<Pd) DD[i+1] = DD[i]-1; // del
      else     DD[i+1] = DD[i];   // trans
    } else {
      //--- Trans/Ins/Del
      if(p<Pi)           DD[i+1] = DD[i]+1; // ins
      else if(p < Pi+Pd) DD[i+1] = DD[i]-1; // del
      else               DD[i+1] = DD[i];   // trans
    } // if DD[i]
    //printf("%03d %+02d %e\n",i,DD[i],p);
  } // for i
  //----- set received word
  for(int i=0;i<N;i++){
    for(int j=DD[i];j<=DD[i+1];j++){
      p = (double) random()/RAND_MAX;
      if(p<Ps){
	v = X[i];
	v = ( v+1+(random()%3) ) %4;  // symmetric
	Y[i+j] = v;
      } else {
	Y[i+j] = X[i];
      } // if
      //Y[i+j] = (p<Ps)? inv(X[i]) : X[i];
    } // for j
  } // for i
  //----- (check)
  if(Pi==0 && Pd==0 && Ps==0){
    assert( DD[N]==0 && HammingDist(X, Y, N)==0 ); 
  } // if
  return N+DD[N]; //TMP
}

//================================================================================
void CIDSchannel::funcFd(unsigned char **YV2, const unsigned char **YV){
  for(int i=0; i<Nseq; i++){
    for(int j=0; j<N+Dmax; j++){
      YV2[i][j] = (YV[i][j]<=1)? 0 : 1;
    } // for j
  } // for i
}

//================================================================================
int    CIDSchannel::GetN(){   return N;}
int    CIDSchannel::GetDmin(){return Dmin;}
int    CIDSchannel::GetDmax(){return Dmax;}
double CIDSchannel::GetPi(){  return Pi;}
double CIDSchannel::GetPd(){  return Pd;}
double CIDSchannel::GetPs(){  return Ps;}

//================================================================================
void CIDSchannel::GetDV(int **dbgDV){
  for(int i=0;i<Nseq;i++){
    memcpy(dbgDV[i],DV[i],sizeof(int)*(N+1));
  } // for i
}

//================================================================================
void CIDSchannel::GetXV(unsigned char **dbgXV){
  for(int i=0;i<Nseq;i++){
    memcpy(dbgXV[i],XV[i],sizeof(unsigned char)*N);
  } // for i
}
