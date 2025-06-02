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
#define LCmax 1000

void   RemoveUnbalancedWords(class InnerCodebook *ICB);
void   RemoveSkewedWords(class InnerCodebook *ICB, int delta);
bool   CheckNumCW(class InnerCodebook *ICB, const int *B2p);
int    PruneCW(class InnerCodebook *ICB, int *B2p);
void   GetICBhist(class InnerCodebook *ICB, int cb, int **HST);
int    CntIns(class InnerCodebook *ICB, int cb, int v, const int *HistL);
int    CntDel(class InnerCodebook *ICB, int cb, int v, const int *HistR);

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
  class InnerCodebook *ICBtmp=NULL;
  class InnerCodebook *ICBmin=NULL;
  int LC = (RndSel)? LCmax : 1;
  int val,minval=INT_MAX;
  for(int i=0;i<LC;i++){
    printf("[%d/%d]\n",i+1,LC);
    if(ICBtmp!=NULL) delete ICBtmp;
    ICBtmp = new class InnerCodebook(ICB);
    val = PruneCW(ICBtmp,B2p); //-----
    if(val<minval){
      minval = val;
      if(ICBmin!=NULL) delete ICBmin;
      ICBmin = new class InnerCodebook(ICBtmp);
    } // if
    printf("val=%d minval=%d\n",val,minval);
  } // for i
  assert(ICBtmp!=NULL && ICBmin!=NULL);
  
  ICBmin->GenMap();
  ICBmin->check();
  ICBmin->CWMwrite(fn);
  //ICBmin->PrintCWM();
  //ICBmin->dump();
  
  delete [] B;
  delete [] B2p;
  delete [] NumCW;
  delete ICB;
  delete ICBtmp;
  delete ICBmin;
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
int PruneCW(class InnerCodebook *ICB, int *B2p){
  int  beta   = ICB->Get_beta();
  int  beta4p = ICB->Get_beta4p();
  int  nu     = ICB->Get_nu();
  int  v,cb,maxval,maxpos,maxval_pre=-1;
  int  cntI,cntD,cnt,dcnt;
  int  v0,cb0, v_rnd, cb_rnd;
  int  *NC = new int [nu];         // number of codewords
  int  ***Hist = new int ** [nu];  // [nu][beta][4]
  for(cb=0;cb<nu;cb++){
    Hist[cb] = new int * [beta];
    for(int i=0;i<beta;i++) Hist[cb][i] = new int [4];
  } // for cb
  unsigned char *V = new unsigned char [beta];
  
  ICB->GetNumCW(NC);
  dcnt = Sum(NC,nu);
  printf(">dcnt=%d\n",dcnt);
  //-----
  while(1){
    maxval = -1;
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
	    cntI = CntIns(ICB,cb,v,Hist[(cb-1+nu)%nu][beta-1]);
	    cntD = CntDel(ICB,cb,v,Hist[(cb+1   )%nu][0     ]);
	    cnt  = cntI+cntD;
	    if(cnt>maxval){
	      //printf("  %d %d %d\n",cb0,v0,cnt);
	      maxval = cnt;
	      maxpos = cb*beta4p + v;
	      //(dbg)
	      // ConvIntVect(v,V,beta);
	      // PrintVect4(V,beta);
	      // printf(" cb=%d cntI=%04d cntD=%04d cnt=%04d[up]\n",cb,cntI,cntD,cnt);
	    } //if
	    //(dbg)
	    // ConvIntVect(v,V,beta);
	    // PrintVect4(V,beta);
	    // printf(" cb=%d cntI=%04d cntD=%04d cnt=%04d\n",cb,cntI,cntD,cnt);

	  } // if ICB
	} // for v0
      } // if NC[]>B2p[]
    } // for cb0
    if(maxpos==-1) break; //-----
    cb = (int)floor((double)maxpos/beta4p);
    v  = maxpos%beta4p;
    assert(ICB->CWMget(v,cb)==1);
    ICB->CWMset(v,cb,0);
    dcnt--;
    maxval_pre = maxval;
    
    //(dbg)
    // printf("---\n");
    // ConvIntVect(v,V,beta);
    // printf("del: ");
    // PrintVect4(V,beta);
    // PrintVect((const int *)NC,nu," NC:"," ");
    // printf(": cb=%d maxval=%d\n",cb,maxval);
  } // while(1)
  //-----
  assert(dcnt==Sum(B2p,nu));
  
  //(dbg)
  printf(">dcnt=%d maxval_pre=%d\n",dcnt,maxval_pre);
  for(cb=0;cb<nu;cb++){
    GetICBhist(ICB,cb,Hist[cb]);
    printf(">cb=%d\n",cb);
    PrintArray((const int **)Hist[cb],beta,4);
  }
  
  // delete
  for(cb=0;cb<nu;cb++){
    for(int i=0;i<beta;i++) delete [] Hist[cb][i];
    delete [] Hist[cb];
  } // for cb
  delete [] Hist;
  delete [] NC;  
  delete [] V;
  return maxval_pre;
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
int CntIns(class InnerCodebook *ICB, int cb, int v, const int *HistL){
  int cnt=0, vx, c;
  int beta  = ICB->Get_beta();
  int nu    = ICB->Get_nu();
  unsigned char *V  = new unsigned char [beta];
  unsigned char *Vx = new unsigned char [beta];
  assert(cb>=0 && cb<nu);
  ConvIntVect(v,V,beta);
  //PrintVect(V,beta,"I:V:","\n");
  
  // pos=-1 -----
  memcpy(&Vx[1],V,sizeof(unsigned char)*(beta-1));
  for(int u=0;u<4;u++){
    Vx[0] = u;
    vx = ConvVectInt(Vx,beta);
    c = HistL[u] * ICB->CWMget(vx,cb);
    cnt+=c;
    //(dbg)
    // PrintVect(Vx,beta," I0:Vx:"," ");
    // printf("c=%d\n",c);
  } // for u
  
  // pos=0...beta-2 -----
  for(int pos=0;pos<beta-1;pos++){
    memcpy(Vx,V,sizeof(unsigned char)*(pos+1));
    Vx[pos+1] = V[pos];
    if(pos<beta-2) memcpy(&Vx[pos+2],&V[pos+1],sizeof(unsigned char)*(beta-pos-2));
    vx = ConvVectInt(Vx,beta);
    c = ICB->CWMget(vx,cb);
    cnt+=c;
    //(dbg)
    // PrintVect(Vx,beta," I1:Vx:"," ");
    // printf("c=%d\n",c);    
  } // for pos
  
  delete [] V;
  delete [] Vx;
  return cnt;
}

//============================================================================
int CntDel(class InnerCodebook *ICB, int cb, int v, const int *HistR){
  int cnt=0, vx, c;
  int beta  = ICB->Get_beta();
  int nu    = ICB->Get_nu();
  unsigned char *V  = new unsigned char [beta];
  unsigned char *Vx = new unsigned char [beta];
  assert(cb>=0 && cb<nu);
  ConvIntVect(v,V,beta);
  //PrintVect(V,beta,"D:V:","\n");

  //----
  for(int pos=0;pos<beta;pos++){
    if(pos>0) memcpy(Vx,V,sizeof(unsigned char)*pos);
    memcpy(&Vx[pos],&V[pos+1],sizeof(unsigned char)*(beta-pos-1));
    for(int u=0;u<4;u++){
      Vx[beta-1] = u;
      vx = ConvVectInt(Vx,beta);
      c = HistR[u] * ICB->CWMget(vx,cb);
      cnt+=c;
      //(dbg)
      // PrintVect(Vx,beta," D:Vx:"," ");
      // printf("c=%d\n",c);    
    } // for u
  } // for pos

  delete [] V;
  delete [] Vx;
  return cnt;
}
