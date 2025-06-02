#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <iostream>

#include "func2.hpp"
#include "combin.hpp"
#include "bmatrix.hpp"

#define BSIZE 8192

void SelectVect(class bmatrix *V, class bmatrix *U, const unsigned char *sel);  // V<-U
void WriteBaseFile(const char *fn, const unsigned char* sel, int beta, int B, class bmatrix *U);

//============================================================================
int main(int argc, char *argv[]){
  int    beta, B, T;
  int    NumVect, B2p;
  int    cnt;
  unsigned long NumSelPtn;
  char   *fn = new char [BSIZE];
  char   *od;
  
  if(argc!=5){
    fprintf(stderr,"Usage: %s <beta> <b> <t> <output.dir>\n",argv[0]);
    return 1;
  }
  beta = atoi(argv[1]);
  B    = atoi(argv[2]);
  T    = atoi(argv[3]);
  od   =      argv[4];
  printf("# beta=%d B=%d T=%d\n",beta,B,T);
  printf("# output dir: %s\n",od);
  assert(beta>0);
  assert(B>=0 && B<=beta);
  assert(T>=0 && T<=beta);
  NumVect = binom(beta,T);
  B2p     = (int)pow(2,B);
  printf("# NumVect=%d B2p=%d\n",NumVect,B2p);
  assert(B2p<=NumVect);

  class bmatrix *U = new class bmatrix(NumVect,beta);   // set of all vectors of weight t
  class combin *CMBgen = new class combin(beta,T);      // generate weight-t vector  
  class combin *CMBsel = new class combin(NumVect,B2p); // select B2p vectors
  unsigned char *Ur = new unsigned char [beta];
  unsigned char *sel= new unsigned char [NumVect];
  
  assert((int)CMBgen->getNum()==NumVect);
  NumSelPtn = CMBsel->getNum();
  printf("# NumSelPtn: %lu\n",NumSelPtn);
  
  // set U
  cnt=0;
  do {
    CMBgen->getV(Ur);
    U->setRVunp(cnt,Ur);
    cnt++;
  } while(CMBgen->nextV());
  printf("[U]\n");
  U->print();

  // generate
  cnt = 0;
  do {
    CMBsel->getV(sel);
    snprintf(fn,BSIZE,"%s/%08d.txt",od,cnt);
    printf("> %s\n",fn);
    WriteBaseFile(fn,sel,beta,B,U);
    cnt++;
  } while(CMBsel->nextV());
  
  delete [] Ur;
  delete [] sel;
  delete [] fn;
  delete U;
  delete CMBgen;
  delete CMBsel;
  return 0;
}

//============================================================================
void SelectVect(class bmatrix *V, class bmatrix *U, const unsigned char *sel){
  int B2p    = V->getM();
  int NumVect= U->getM();
  int beta   = V->getN();
  int cnt=0;
  assert(U->getN()==beta);
  for(int i=0;i<NumVect;i++){
    if(sel[i]==1){
      V->setV(U, cnt,0,i,0, 1,beta);
      cnt++;
    } // if sel
  } // for i
  assert(cnt==B2p);
}

//============================================================================
void WriteBaseFile(const char *fn, const unsigned char* sel, int beta, int B, class bmatrix *U){
  int  B2p = (int)pow(2,B);
  unsigned int x;
  unsigned char *X = new unsigned char [beta];
  FILE *fp;
  class bmatrix *V = new class bmatrix(B2p,beta);

  if((fp=fopen(fn,"w"))==NULL){
    fprintf(stderr,"Cannot open %s\n",fn);
    exit(1);
  }
  // L1: beta nu
  fprintf(fp,"%d 1\n",beta);

  // L2: b
  fprintf(fp,"%d\n",beta+B);

  // L3: cw
  SelectVect(V,U,sel);
  for(int i=0;i<B2p;i++){
    V->getRVunp(X,i);
    x = ConvVectUint(X,beta);
    fprintf(fp,"%u ",x);
  } // for i
  fprintf(fp,"\n");
  
  fclose(fp);
  delete [] X;
  delete V;
}
