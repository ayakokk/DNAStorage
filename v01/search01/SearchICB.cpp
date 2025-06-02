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
void   PruneCW(class InnerCodebook *ICB, int *B2p);
double CalcTmpEntropy(const int **HST, int v, int beta);
void   GetICBhist(class InnerCodebook *ICB, int cb, int **HST);
int    CalcDiff(const int ***Hist, int nu, int beta, int cb, int v);

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
  PruneCW(ICB,B2p);

  ICB->GenMap();
  ICB->check();
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
void PruneCW(class InnerCodebook *ICB, int *B2p){
  int  beta   = ICB->Get_beta();
  int  beta4p = ICB->Get_beta4p();
  int  nu     = ICB->Get_nu();
  int  v,cb,diff,maxval,maxpos;
  int  *NC = new int [nu];         // number of codewords
  int  ***Hist = new int ** [nu];  // [nu][beta][4]
  for(int cb=0;cb<nu;cb++){
    Hist[cb] = new int * [beta];
    for(int i=0;i<beta;i++) Hist[cb][i] = new int [4];
  } // for cb
  
  //-----
  while(1){
    maxval = INT_MIN;
    maxpos = -1;
    ICB->GetNumCW(NC);
    for(cb=0;cb<nu;cb++) GetICBhist(ICB,cb,Hist[cb]);
    
    for(cb=0;cb<nu;cb++){
      if(NC[cb]>B2p[cb]){
	for(v=0;v<beta4p;v++){
	  if(ICB->CWMget(v,cb)==1){
	    diff = CalcDiff((const int ***)Hist,nu,beta,cb,v);
	    if(diff>maxval){
	      maxval = diff;
	      maxpos = cb*beta4p + v;
	    } // if diff
	  } // if ICB
	} // for v
      } // if NC[]>B2p[]
    } // for cb
    if(maxpos==-1) break; //-----
    cb = (int)floor((double)maxpos/beta4p);
    v  = maxpos%beta4p;
    assert(ICB->CWMget(v,cb)==1);
    ICB->CWMset(v,cb,0);

    //(dbg)
    
    printf("maxval=%d pos=(%d,%04X) ",maxval,cb,v);
    PrintVect((const int *)NC,nu,"NC:","\n");
    for(int i=0;i<nu;i++){
      printf("cb=%d\n",i);
      PrintArray((const int **)Hist[i],beta,4);
    }
    
  } // while(1)
  //-----
  for(int i=0;i<nu;i++){
    printf("cb=%d\n",i);
    PrintArray((const int **)Hist[i],beta,4);
  }
  
  // delete
  for(cb=0;cb<nu;cb++){
    for(int i=0;i<beta;i++) delete [] Hist[cb][i];
    delete [] Hist[cb];
  } // for cb
  delete [] Hist;
  delete [] NC;  
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

//============================================================================
int CalcDiff(const int ***Hist, int nu, int beta, int cb, int v){
  int cbp = (cb-1+nu)%nu;
  int cbn = (cb+1)%nu;
  int h0d=0,h1d=0;
  unsigned char *V = new unsigned char [beta];
  int **H0 = new int * [beta+2];
  int **H1 = new int * [beta+2];
  for(int i=0;i<beta+2;i++){
    H0[i] = new int [4];
    H1[i] = new int [4];
  } // for i
  assert(cb>=0 && cb<nu);

  // Set H0
  for(int j=0;j<4;j++) H0[0     ][j] = Hist[cbp][beta-1][j]; // last of previous CB
  for(int j=0;j<4;j++) H0[beta+1][j] = Hist[cbn][0     ][j]; // head of next CB
  for(int i=0;i<beta;i++)
    for(int j=0;j<4;j++) H0[i+1][j] = Hist[cb][i][j];

  // Set H1
  ConvIntVect(v,V,beta);
  for(int i=0;i<beta+2;i++) memcpy(H1[i],H0[i],sizeof(int)*4);
  for(int i=0;i<beta;i++){
    assert(H1[i+1][V[i]]>0);
    H1[i+1][V[i]]--;
  } // for i

  // diff
  for(int i=0;i<beta+1;i++){
    for(int j=0;j<4;j++){
      //h0d += (H0[i+1][j]-H0[i][j])*(H0[i+1][j]-H0[i][j]);
      //h1d += (H1[i+1][j]-H1[i][j])*(H1[i+1][j]-H1[i][j]);
      h0d += (int)abs((int)pow(H0[i+1][j]-H0[i][j],0.5));
      h1d += (int)abs((int)pow(H1[i+1][j]-H1[i][j],0.5));
    } // for j
  } // for i
  
  //(dbg)
  /*
  printf("nu=%d beta=%d cb=%d v=%04X: h0d=%d h1d=%d\n",nu,beta,cb,v,h0d,h1d);
  PrintVect4(V,beta); printf("\n");
  printf("H0\n");
  PrintArray((const int **)H0,beta+2,4);
  printf("H1\n");
  PrintArray((const int **)H1,beta+2,4);
  */
  
  for(int i=0;i<beta+2;i++){
    delete [] H0[i];
    delete [] H1[i];    
  } // for i
  delete [] H0;
  delete [] H1;
  delete [] V;
  return h1d-h0d;
}
