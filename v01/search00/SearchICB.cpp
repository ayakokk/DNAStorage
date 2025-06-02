#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <assert.h>
#include <iostream>

#include "func.hpp"
#include "bmatrix.hpp"
#include "InnerCodebook.hpp"

#define BSIZE 8192

void   RemoveUnbalancedWords(class InnerCodebook *ICB);
void   RemoveSkewedWords(class InnerCodebook *ICB, int delta);
bool   CheckNumCW(class InnerCodebook *ICB, const int *B2p);
void   PruneCW(class InnerCodebook *ICB, int cb, int b2p);
double CalcTmpEntropy(const int **HST, int v, int beta);
void   GetICBhist(class InnerCodebook *ICB, int cb, int **HST);

//============================================================================
int main(int argc, char *argv[]){
  int  beta, nu, delta, seed;
  int  *B, *B2p, *NumCW;
  char *fn;
  if(argc<7){
    fprintf(stderr,"Usage: %s <beta> <nu> <delta> <ICBout.bin> <b0> ... <bnu-1> <seed|-1>\n",argv[0]);
    return 1;
  }
  beta = atoi(argv[1]);
  nu   = atoi(argv[2]);
  delta= atoi(argv[3]);
  fn   =      argv[4];
  assert(beta>0 && nu>0 && delta>=0);
  assert(argc==6+nu);
  B    = new int [nu];
  B2p  = new int [nu];
  NumCW= new int [nu]; 
  for(int i=0;i<nu;i++){
    B[i]  = atoi(argv[5+i]);
    B2p[i]= (int)pow(2,B[i]);
    assert(B[i]<=2*beta);
  } // for i
  seed = atoi(argv[argc-1]);
  if(seed==-1) seed = (int)time(NULL);
  srandom(seed);
  printf("# beta=%d nu=%d delta=%d [%d]\n",beta,nu,delta,seed);
  printf("# output: %s\n",fn);
  PrintVect(B,  nu,"# B:   ","\n");
  PrintVect(B2p,nu,"# B2p: ","\n");

  class InnerCodebook *ICB = new class InnerCodebook(beta,B,nu);
  ICB->CWMset();  // select all

  // (2) balance --------------------
  RemoveUnbalancedWords(ICB);    
  ICB->GetNumCW(NumCW);
  PrintVect(NumCW,nu,"(2) ","\n");

  // (3) skew --------------------
  RemoveSkewedWords(ICB,delta); 
  ICB->GetNumCW(NumCW);
  PrintVect(NumCW,nu,"(3) ","\n");

  // (4) check num CW--------------------
  assert(CheckNumCW(ICB,B2p));

  // (5) prune --------------------
  for(int cb=0;cb<nu;cb++){
    PruneCW(ICB,cb,B2p[cb]);
  } // for cb

  ICB->GenMap();
  ICB->CWMwrite(fn);
  //ICB->PrintCWM();
  //ICB->dump();
  
  delete [] B;
  delete [] B2p;
  delete [] NumCW;
  delete ICB;
  return 0;
}

//============================================================================
void RemoveUnbalancedWords(class InnerCodebook *ICB){
  int beta   = ICB->Get_beta();
  int beta4p = ICB->Get_beta4p();
  int nu     = ICB->Get_nu();
  int w0 = (int)floor((double)beta/2.0);
  int w1 = (int)ceil( (double)beta/2.0);
  int w;

  unsigned char *V = new unsigned char [beta];
  for(int v=0;v<beta4p;v++){
    ConvIntVect(v,V,beta);
    w = GCweight(V,beta);
    for(int cb=0;cb<nu;cb+=2)
      if(w!=w0) ICB->CWMset(v,cb,0);
    for(int cb=1;cb<nu;cb+=2)
      if(w!=w1) ICB->CWMset(v,cb,0);
  } // for v
  delete [] V;
}

//============================================================================
void RemoveSkewedWords(class InnerCodebook *ICB, int delta){
  int beta   = ICB->Get_beta();
  int beta4p = ICB->Get_beta4p();
  int nu     = ICB->Get_nu();
  int ell = (int)floor((double)beta/2.0);
  int wL,wR;
  
  unsigned char *V = new unsigned char [beta];
  for(int v=0;v<beta4p;v++){
    ConvIntVect(v,V,beta);
    wL = GCweight(V,           ell);
    wR = GCweight(&V[beta-ell],ell);
    if(abs(wL-wR)>delta){
      for(int cb=0;cb<nu;cb++) ICB->CWMset(v,cb,0);
    } // if 
  } // for v
  delete [] V;
}

//============================================================================
bool CheckNumCW(class InnerCodebook *ICB, const int *B2p){
  bool ret=true;
  int  nu = ICB->Get_nu();
  int  *NumCW = new int [nu]; 
  ICB->GetNumCW(NumCW);

  for(int cb=0;cb<nu;cb++){
    if(NumCW[cb]<B2p[cb]){
      printf("failed\n");
      ret=false;
      break;
    } // if
  } // for cb

  delete [] NumCW;
  return ret;
}

//============================================================================
void PruneCW(class InnerCodebook *ICB, int cb, int b2p){
  int    beta   = ICB->Get_beta();
  int    beta4p = ICB->Get_beta4p();
  int    nu = ICB->Get_nu();
  int    nc = ICB->GetNumCW(cb);
  int    v,   vmin;
  double ent, emin;
  assert(cb>=0 && cb<nu);
  assert(nc>=b2p);
  int **hst = new int * [beta];
  for(int i=0;i<beta;i++) hst[i] = new int [4];

  for(int cnt=0;cnt<nc-b2p;cnt++){
    GetICBhist(ICB,cb,hst);
    //PrintArray((const int **)hst,beta,4);
    vmin = -1;
    emin = beta*2;
    for(v=0;v<beta4p;v++){
      if(ICB->CWMget(v,cb)==1){
	ent = CalcTmpEntropy((const int **)hst,v,beta);
	if(ent<emin){
	  emin = ent;
	  vmin = v;
	} // if ent<emin
      } // if ICB
    } // for v
    assert(vmin>=0 && ICB->CWMget(vmin,cb)==1);
    ICB->CWMset(vmin,cb,0);
    //printf("vmin=%d emin=%e\n",vmin,emin);
  } // for cnt
  
  for(int i=0;i<beta;i++) delete [] hst[i];
  delete [] hst;
}

//============================================================================
double CalcTmpEntropy(const int **HST, int v, int beta){
  double ent=0.0;
  unsigned char *V = new unsigned char [beta];
  int *hst0= new int [4];
  ConvIntVect(v,V,beta);

  for(int i=0;i<beta;i++){
    memcpy(hst0,HST[i],sizeof(int)*4);
    assert(V[i]>=0 && V[i]<4);
    assert(hst0[V[i]]>0);
    hst0[V[i]]--;
    ent += Entropy((const int *)hst0,4);
  } // for i
  
  //PrintVect4((const unsigned char *)V,beta);
  //printf(": %e\n",ent);
  delete [] hst0;
  delete [] V;
  return ent;//TMP
}


//============================================================================
void GetICBhist(class InnerCodebook *ICB, int cb, int **HST){
  int beta  = ICB->Get_beta();
  int beta4p= ICB->Get_beta4p();
  int nu    = ICB->Get_nu();
  unsigned char *V = new unsigned char [beta];
  assert(cb>=0 && cb<nu);
  
  ClearArray(HST,beta,4);
  for(int v=0;v<beta4p;v++){
    if(ICB->CWMget(v,cb)==1){
      ConvIntVect(v,V,beta);
      for(int i=0;i<beta;i++) HST[i][V[i]]++;
    } // if ICB->CWMget()
  } // for v

  delete [] V;
}
