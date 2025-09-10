#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "IDSchannel.hpp"

#define Drange 2.0 // Dmin=-Drange*N*Pd Dmax=Drange*N*Pi  

//================================================================================
int IDSchannel::max(int a, int b){return (a>b)? a : b; }
int IDSchannel::min(int a, int b){return (a<b)? a : b; }  

//================================================================================
unsigned char IDSchannel::inv(unsigned char a){
  assert(a==0 || a==1);
  return (a==0)? 1 : 0 ;
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
IDSchannel::IDSchannel(int _N, double _Pi, double _Pd, double _Ps){
  N  = _N;
  Pi = _Pi;
  Pd = _Pd;
  Ps = _Ps;
  // Dmin = (int)floor(-(double)Drange*N*Pd);
  // Dmax = (int)ceil ( (double)Drange*N*Pi);
  Dmin = -N;
  Dmax = N * 2;
  printf("# IDSchannel: N=%d (Pi,Pd,Ps)=(%e,%e,%e) (Dmin,Dmax)=(%d,%d)\n",N,Pi,Pd,Ps,Dmin,Dmax);
  assert( N>0 );
  assert( Pi>=0.0 && Pi<0.5 );
  assert( Pd>=0.0 && Pd<0.5 );
  assert( Ps>=0.0 && Ps<0.5 );
  DR = new int [N+1];
}

//================================================================================
IDSchannel::~IDSchannel(){
  delete [] DR;
  printf("# IDSchannel: deleted\n");
}

//================================================================================
int IDSchannel::transmit(unsigned char *Y, const unsigned char *X){
  double p;
  //----- set drift vector
  DR[0] = 0;
  for(int i=0;i<N;i++){
    p = (double) random()/RAND_MAX;
    if(DR[i]==Dmin){
      //--- Trans or Ins
      if(p<Pi) DR[i+1] = DR[i]+1; // ins
      else     DR[i+1] = DR[i];   // trans
    } else if(DR[i]==Dmax){
      //--- Trans or Del
      if(p<Pd) DR[i+1] = DR[i]-1; // del
      else     DR[i+1] = DR[i];   // trans
    } else {
      //--- Trans/Ins/Del
      if(p<Pi)           DR[i+1] = DR[i]+1; // ins
      else if(p < Pi+Pd) DR[i+1] = DR[i]-1; // del
      else               DR[i+1] = DR[i];   // trans
    } // if DR[i]
    //printf("%03d %+02d %e\n",i,DR[i],p);
  } // for i
  //----- set received word
  for(int i=0;i<N;i++){
    for(int j=DR[i];j<=DR[i+1];j++){
      p = (double) random()/RAND_MAX;
      Y[i+j] = (p<Ps)? inv(X[i]) : X[i];
    } // for j
  } // for i
  return N+DR[N]; //TMP
}

//================================================================================
int    IDSchannel::GetN(){   return N;}
int    IDSchannel::GetDmin(){return Dmin;}
int    IDSchannel::GetDmax(){return Dmax;}
double IDSchannel::GetPi(){  return Pi;}
double IDSchannel::GetPd(){  return Pd;}
double IDSchannel::GetPs(){  return Ps;}

//================================================================================
void IDSchannel::GetDR(int *dbgDR){
  memcpy(dbgDR,DR,sizeof(int)*(N+1));
}
