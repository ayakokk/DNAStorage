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
void   GetICBhist(class InnerCodebook *ICB, int cb, int **HST);
int    CntDist0(class InnerCodebook *ICB, int cb, int v0, int offset, int len);

bool   RndSel;

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
  RndSel = (seed!=0);
  printf("# beta=%d nu=%d delta=%d [%d:%d]\n",beta,nu,delta,seed,RndSel);
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
  ICB->PrintCWM();
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
  int  v,cb,cntL,cntR,xL,xR,cntx,maxval,maxpos;
  int  v0,cb0, v_rnd, cb_rnd;
  int  *NC = new int [nu];         // number of codewords
  int  ***Hist = new int ** [nu];  // [nu][beta][4]
  for(cb=0;cb<nu;cb++){
    Hist[cb] = new int * [beta];
    for(int i=0;i<beta;i++) Hist[cb][i] = new int [4];
  } // for cb
  unsigned char *V = new unsigned char [beta];

  //-----
  while(1){
    maxval = 0;
    maxpos = -1;
    ICB->GetNumCW(NC);
    for(cb=0;cb<nu;cb++) GetICBhist(ICB,cb,Hist[cb]);
    cb_rnd = (RndSel)? random() : 0;
    v_rnd  = (RndSel)? random() : 0;
    
    for(cb0=0;cb0<nu;cb0++){
      cb = (cb0+cb_rnd)%nu;
      if(NC[cb]>B2p[cb]){
	for(v0=0;v0<beta4p;v0++){
	  v = (v0+v_rnd)%beta4p;
	  if(ICB->CWMget(v,cb)==1){
	    ConvIntVect(v,V,beta);
	    cntL = CntDist0(ICB,cb,v,-1,beta-1);
	    cntR = CntDist0(ICB,cb,v,+1,beta-1);
	    xL = Hist[(cb-1+nu)%nu][beta-1][V[0]     ];
	    xR = Hist[(cb+1   )%nu][0     ][V[beta-1]];
	    cntx = cntL*xL + cntR*xR;
	    if(cntx>maxval){
	      maxval = cntx;
	      maxpos = cb*beta4p + v;
	      //(dbg)
	      //PrintVect4(V,beta);
	      //printf(" cb=%d cntL=%04d cntR=%04d xL=%d xR=%d cntx=%d\n",cb,cntL,cntR,xL,xR,cntx);
	    } //if
	    
	    //(dbg)
	    // PrintVect(V,beta,""," ");
	    // printf("cb=%d cntL=%04d cntR=%04d xL=%d xR=%d\n",cb,cntL,cntR,xL,xR);
	  } // if ICB
	} // for v0
      } // if NC[]>B2p[]
    } // for cb0
    if(maxpos==-1) break; //-----
    cb = (int)floor((double)maxpos/beta4p);
    v  = maxpos%beta4p;
    assert(ICB->CWMget(v,cb)==1);
    ICB->CWMset(v,cb,0);
    
    //(dbg)
    ConvIntVect(v,V,beta);
    printf("del: ");
    PrintVect4(V,beta);
    PrintVect((const int *)NC,nu," NC:"," ");
    printf(": cb=%d maxval=%d\n",cb,maxval);
  } // while(1)
  //-----

  // delete
  for(cb=0;cb<nu;cb++){
    for(int i=0;i<beta;i++) delete [] Hist[cb][i];
    delete [] Hist[cb];
  } // for cb
  delete [] Hist;
  delete [] NC;  
  delete [] V;
  
  /*
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
  */
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
int CntDist0(class InnerCodebook *ICB, int cb, int v0, int offset, int len){
  int cnt=0, i0,i1,v1,dist;
  int beta  = ICB->Get_beta();
  int beta4p= ICB->Get_beta4p();
  int nu    = ICB->Get_nu();
  unsigned char *V0 = new unsigned char [beta];
  unsigned char *V1 = new unsigned char [beta];
  assert(cb>=0 && cb<nu);
  if(offset<=0){
    i0 = -offset;
    i1 = 0;
  } else {
    i0 = 0;
    i1 = offset;
  } // if offset
  assert(len>=0 && len+i1<=beta);

  ConvIntVect(v0,V0,beta);
  for(v1=0;v1<beta4p;v1++){
    if(ICB->CWMget(v1,cb)==1){
      ConvIntVect(v1,V1,beta);
      dist = HammingDist(&V0[i0],&V1[i1],len); 
      if(dist==0) cnt++;

      //(dbg)
      // PrintVect(V0,beta," "," ");
      // PrintVect(V1,beta," "," ");
      // printf("dist=%d: cb=%d offset=%d len=%d\n",dist,cb,offset,len);
    } // if ICB
  } // for v1
  
  delete [] V0;
  delete [] V1;
  return cnt;
}
