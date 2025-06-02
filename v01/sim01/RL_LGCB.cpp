#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "bmatrix.hpp"
#include "func.hpp"
#include "InnerCodebook.hpp"

#define _Wmax 20

int nu,beta,beta4p,Q;
int Wmin,Wmax;
class InnerCodebook *ICB;

void CalcRunLength();
int  RunLenL(const unsigned char *X, unsigned char v);
int  RunLenR(const unsigned char *X, unsigned char v);
int  RunLenC(const unsigned char *X, unsigned char v);
void MaxRunLenSeg(int cb, unsigned char v, int *lmaxL, int *lmaxC, int *lmaxR);   
int  MaxRunLen(unsigned char v, const int **runL, const int **runC, const int **runR);

void CalcLocalGCB();
int  CalcWeight(const unsigned char *X, int len);
int  CalcWeightMax(int posL, int len, const int **maxL, const int **maxR);
int  CalcWeightMin(int posL, int len, const int **minL, const int **minR);

//=================================================================================
int main(int argc, char *argv[]){
  char *fn;
  Q=4; // fixed
  if(argc!=2){
    fprintf(stderr,"Usage: %s <ICB.bin>\n",argv[0]);
    return 1;
  }
  fn = argv[1];
  printf("# ICB: %s\n",fn);
  ICB   = new class InnerCodebook(fn);
  nu    = ICB->Get_nu();
  beta  = ICB->Get_beta();
  beta4p= ICB->Get_beta4p();
  Wmin  = beta;
  Wmax  = _Wmax;
  printf("# nu=%d beta=%d(%d) Q=%d Wmin=%d Wmax=%d\n",nu,beta,beta4p,Q,Wmin,Wmax);

  CalcRunLength();
  CalcLocalGCB();
  
  delete ICB;
  
  return 0;
}

//=================================================================================
void CalcRunLength(){
  int lmax;
  int *len = new int [Q]; 
  int **runL = new int * [nu];
  int **runC = new int * [nu];
  int **runR = new int * [nu];
  for(int i=0;i<nu;i++){
    runL[i] = new int [Q];
    runC[i] = new int [Q];
    runR[i] = new int [Q];
  } // for i

  //--- set
  for(int cb=0;cb<nu;cb++){
    for(unsigned char v=0;v<Q;v++){
      MaxRunLenSeg(cb,v,&runL[cb][v],&runC[cb][v],&runR[cb][v]);
    } // for v
  } // for cb

  //--- calc
  for(unsigned char v=0;v<Q;v++){
    len[v] = MaxRunLen(v,(const int **)runL, (const int **)runC, (const int **)runR);
  } // for v
  lmax = max(len,Q);
  printf(">>> max run length: %d\n",lmax);
  
  //(dbg)
  /*
  printf("runL\n");
  PrintArray((const int**)runL,nu,Q);
  printf("runC\n");
  PrintArray((const int**)runC,nu,Q);
  printf("runR\n");
  PrintArray((const int**)runR,nu,Q);
  PrintVect(len,Q,"len:","\n");
  */
  
  //---
  delete [] len;
  for(int i=0;i<nu;i++){
    delete [] runL[i];
    delete [] runC[i];
    delete [] runR[i];
  } // for i
  delete [] runL;
  delete [] runC;
  delete [] runR;
}

//=================================================================================
int RunLenL(const unsigned char *X, unsigned char v){
  assert(v<Q);
  int len=0;
  for(int i=0;i<beta;i++){
    if(X[i]==v) len++;
    else break;
  } // for 
  return len;
}

//=================================================================================
int RunLenR(const unsigned char *X, unsigned char v){
  assert(v<Q);
  int len=0;
  for(int i=beta-1;i>=0;i--){
    if(X[i]==v) len++;
    else break;
  } // for 
  return len;
}

//=================================================================================
int RunLenC(const unsigned char *X, unsigned char v){
  assert(v<Q);
  int len, lmax=0;
  len = (X[0]==v)? 1:0;
  lmax= len;
  for(int i=1;i<beta;i++){
    if(X[i]==v && X[i-1]==v){
      len++;
      if(len>lmax) lmax=len;
    } else if(X[i]==v){
      len=1;
      if(len>lmax) lmax=len;
    } else {
      len=0;
    } // if
  } // for i
  return lmax;
}

//=================================================================================
void MaxRunLenSeg(int cb, unsigned char v, int *lmaxL, int *lmaxC, int *lmaxR){
  assert(cb>=0 && cb<nu);
  assert(v<Q);
  int lenL,lenC,lenR;
  unsigned char *X = new unsigned char [beta];
  *lmaxL=0;
  *lmaxC=0;
  *lmaxR=0;
  
  //printf("[cb:%d v:%u]\n",cb,v);
  for(int x=0;x<beta4p;x++){
    if(ICB->CWMget(x,cb)==1){
      ConvIntVect(x,X,beta);
      lenL = RunLenL(X,v);
      lenC = RunLenC(X,v);
      lenR = RunLenR(X,v);
      if(lenL>(*lmaxL)) (*lmaxL)=lenL;
      if(lenC>(*lmaxC)) (*lmaxC)=lenC;
      if(lenR>(*lmaxR)) (*lmaxR)=lenR;
      //PrintVect(X,beta); printf(": %d %d\n",len,lmax);
    } // if ICB
  } // for x
  
  delete [] X;
}

//=================================================================================
int MaxRunLen(unsigned char v, const int **runL, const int **runC, const int **runR){
  assert(v<Q);
  int len,lmax=0;
  for(int L=0;L<nu;L++){
    if(runC[L][v]>lmax) lmax=runC[L][v];
    len = runR[L][v];
    if(len>0){
      for(int R=L+1;R<2*nu;R++){
	len+= runL[R%nu][v];
	if(runL[R%nu][v]<beta) break;
      } // for R
      if(len>lmax) lmax=len;
    } // if len
  } // for L
  return lmax;
}

//=================================================================================
//=================================================================================
//=================================================================================

//=================================================================================
void CalcLocalGCB(){
  int **maxL, **minL, **maxR, **minR;
  int wL,wR;
  unsigned char *X = new unsigned char [beta];
  int Nseg = 2*(int)ceil((double)Wmax/beta);
  maxL = new int * [Nseg];
  minL = new int * [Nseg];
  maxR = new int * [Nseg];
  minR = new int * [Nseg];
  for(int i=0;i<Nseg;i++){
    maxL[i] = new int [beta];
    minL[i] = new int [beta];
    maxR[i] = new int [beta];
    minR[i] = new int [beta];
  } // for i
  
  //--- set: maxL,minL,maxR,minR
  for(int idx=0;idx<Nseg;idx++){
    //printf("idx=%d\n",idx);
    // init
    for(int j=0;j<beta;j++){
      maxL[idx][j]=0;
      minL[idx][j]=beta;
      maxR[idx][j]=0;
      minR[idx][j]=beta;
    } // for j
    // set 
    for(int x=0;x<beta4p;x++){
      if(ICB->CWMget(x,idx%nu)==1){
	ConvIntVect(x,X,beta);
	//PrintVect(X,beta); printf("\n"); 
	for(int j=0;j<beta;j++){
	  wL = CalcWeight(X,j+1);
	  wR = CalcWeight(&X[j],beta-j);
	  if(wL>maxL[idx][j]) maxL[idx][j] = wL;
	  if(wL<minL[idx][j]) minL[idx][j] = wL;
	  if(wR>maxR[idx][j]) maxR[idx][j] = wR;
	  if(wR<minR[idx][j]) minR[idx][j] = wR;
	  //printf("(%d %d) ",wL,wR);
	} // for j
	//printf("\n");
      } // if ICB
    } // for x
  
  } // for idx

  //--- calc
  int Xmax, Xmin;
  int d0,d1,dmax;
  for(int ws=Wmin;ws<=Wmax;ws++){
    dmax = 0;
    for(int pos=0;pos<beta*Nseg/2;pos++){
      Xmax = CalcWeightMax(pos,ws,(const int **)maxL, (const int **)maxR);
      Xmin = CalcWeightMin(pos,ws,(const int **)minL, (const int **)minR);
      d0 = abs( Xmax-(ws-Xmax) );
      d1 = abs( Xmin-(ws-Xmin) );
      if(d0>dmax) dmax=d0;
      if(d1>dmax) dmax=d1;
      //printf("ws=%d pos=%d max=%d min=%d d0=%d d1=%d dmax=%d\n",ws,pos,Xmax,Xmin,d0,d1,dmax);
    } // for pos
    printf(">>> ws=%02d dmax=%d\n",ws,dmax);
  }  // for ws

  
  //(dbg)
  /*
  printf("Wmax=%d Nseg=%d\n",Wmax,Nseg);
  printf("maxL\n");
  PrintArray((const int **)maxL,Nseg,beta);
  printf("minL\n");
  PrintArray((const int **)minL,Nseg,beta);
  printf("maxR\n");
  PrintArray((const int **)maxR,Nseg,beta);
  printf("minR\n");
  PrintArray((const int **)minR,Nseg,beta);
  */
  
  //---
  for(int i=0;i<Nseg;i++){
    delete maxL[i];
    delete minL[i];
    delete maxR[i];
    delete minR[i];
  } // for i
  delete [] maxL;
  delete [] minL;
  delete [] maxR;
  delete [] minR;
  delete [] X;
}

//=================================================================================
int CalcWeight(const unsigned char *X, int len){
  assert(len>=0);
  int w=0;
  for(int i=0;i<len;i++){
    assert(X[i]>=0 && X[i]<4);
    if(X[i]==2 || X[i]==3) w++;
  } // for i
  return w;
}

//=================================================================================
int CalcWeightMax(int posL, int len, const int **maxL, const int **maxR){
  assert(posL>=0 && len>=beta);
  int w;
  int posR = posL+len-1;
  int idxL = (int)floor((double)posL/beta);
  int idxR = (int)floor((double)posR/beta);

  if(idxL==idxR){
    w = maxR[idxL][0];
  } else {
    w = maxR[idxL][posL%beta] + maxL[idxR][posR%beta];
    for(int i=idxL+1;i<idxR;i++) w+= maxR[i][0];
  } // if idxL==idxR
  
  return w;
}

//=================================================================================
int CalcWeightMin(int posL, int len, const int **minL, const int **minR){
  assert(posL>=0 && len>=beta);
  int w;
  int posR = posL+len-1;
  int idxL = (int)floor((double)posL/beta);
  int idxR = (int)floor((double)posR/beta);

  if(idxL==idxR){
    w = minR[idxL][0];
  } else {
    w = minR[idxL][posL%beta] + minL[idxR][posR%beta];
    for(int i=idxL+1;i<idxR;i++) w+= minR[i][0];
  } // if idxL==idxR
  
  return w;
}
