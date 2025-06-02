#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"
#include "InnerCodec.hpp"

#define _Dcoef 300  // Dmin=-_Dcoef*Pd Dmax=_Dcoef*Pi

void dbgPrint(const int *IW, const int *CW, const int *RW, const double **Px, int N, int N2);

//================================================================================
int main(int argc, char *argv[]){
  int RunK, Delta, N, Nmax, N2, seed;
  int Dmin, Dmax;
  double Pi,Pd,Ps;
  if(argc!=8){
    fprintf(stderr,"Usage: %s <Run:K> <LB:delta> <Pi> <Pd> <Ps> <N> <seed|-1>\n",argv[0]);
    return 1;
  } // if
  RunK  = atoi(argv[1]);
  Delta = atoi(argv[2]);
  Pi    = atof(argv[3]);
  Pd    = atof(argv[4]);
  Ps    = atof(argv[5]);
  N     = atoi(argv[6]);
  seed  = atoi(argv[7]);
  if(seed==-1) seed=(int)time(NULL);
  srandom(seed);
  Dmin = -(int)round(Pd*_Dcoef);
  Dmax =  (int)round(Pi*_Dcoef);
  Nmax = N + Dmax;
  printf("# RunK=%d Delta=%d (Pi,Pd,Ps)=(%e,%e,%e) (N,Nmax)=(%d,%d) (Dmin,Dmax)=(%d,%d) [%d]\n",
	 RunK,Delta,Pi,Pd,Ps,N,Nmax,Dmin,Dmax,seed);
  assert(RunK>0 && Delta>0);
  assert(N>0);
  
  class InnerCodec *ICD = new class InnerCodec(RunK,Delta,Pi,Pd,Ps,N,Dmin,Dmax);
  double **Px   = new double * [N];
  double **Pout = new double * [N];
  for(int i=0;i<N;i++){
    Px[i]   = new double [4];
    Pout[i] = new double [4];
  } // for i
  int *IW = new int [N];
  int *CW = new int [N];
  int *RW = new int [Nmax];
  
  //-----
  RandVect(IW,N,0,3);
  ICD->Encode(CW,IW);
  memcpy(RW,CW,sizeof(int)*N); //dummy
  N2 = N; // dummy
  PxUnif(Px,N,4);
  ICD->Decode(Pout,RW,(const double **)Px,N2,CW);  // CW for dbg
  
  // (dbg)
  // PrintVect(IW,N,"IW\n","");
  // PrintVect(CW,N,"CW\n","");
  // dbgPrint(IW,CW,RW,(const double **)Px,N,N2);
  //-----

  delete ICD;
  for(int i=0;i<N;i++){
    delete [] Px[i];
    delete [] Pout[i];
  } // for i
  delete [] Px;
  delete [] Pout;
  delete [] IW;
  delete [] CW;
  delete [] RW;
  return 0;
}

//================================================================================
void dbgPrint(const int *IW, const int *CW, const int *RW, const double **Px, int N, int N2){
  for(int i=0;i<max(N,N2);i++){
    printf("%04d: ",i);
    if(i<N) printf("%d %d ",IW[i],CW[i]);
    else    printf("* * ");
    //---
    if(i<N2) printf("%d ",RW[i]);
    else     printf("* ");
    //---
    if(i<N) PrintVect1(Px[i],4,"[","]");
    printf("\n");
  } // for i
}
