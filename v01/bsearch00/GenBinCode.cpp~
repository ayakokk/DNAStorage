#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <iostream>
#include <list>
#include <algorithm>

#include "func.hpp"
#include "combin.hpp"
#include "MetricCalc.hpp"
#include "bmatrix.hpp"

double Pi,Pd,Ps;

void   SelectVect(class bmatrix *V, class bmatrix *U, const unsigned char *sel);  // V<-U
double CalcMaxPrb(class bmatrix *V, class MetricCalc *MC);
void   PrintList(std::list<unsigned char*> LS, int len);
void   WriteBaseFile(const char *fn, std::list<unsigned char*> LS, int beta, int B, double pmin,
		   class bmatrix *U, class MetricCalc *MC);

//============================================================================
int main(int argc, char *argv[]){
  int    beta, B, T;
  int    NumVect, B2p;
  int    cnt;
  double prob,pmin;
  unsigned long NumSelPtn;
  char   *fn;
  
  if(argc!=8){
    fprintf(stderr,"Usage: %s <beta> <b> <t> <Pi> <Pd> <Ps> <output.txt>\n",argv[0]);
    return 1;
  }
  beta = atoi(argv[1]);
  B    = atoi(argv[2]);
  T    = atoi(argv[3]);
  Pi   = atof(argv[4]);
  Pd   = atof(argv[4]);
  Ps   = atof(argv[4]);
  fn   =      argv[7];
  printf("# beta=%d B=%d T=%d (Pi,Pd,Ps)=(%e,%e,%e)\n",beta,B,T,Pi,Pd,Ps);
  printf("# output: %s\n",fn);
  assert(beta>0);
  assert(B>=0 && B<=beta);
  assert(T>=0 && T<=beta);
  assert(Pi>=0 && Pi<0.5);
  assert(Pd>=0 && Pd<0.5);
  assert(Ps>=0 && Ps<0.5);
  NumVect = binom(beta,T);
  B2p     = (int)pow(2,B);
  printf("# NumVect=%d B2p=%d\n",NumVect,B2p);
  assert(B2p<=NumVect);

  class bmatrix *U = new class bmatrix(NumVect,beta);   // set of all vectors of weight t
  class bmatrix *V = new class bmatrix(B2p,    beta);   // selected B2p vectors
  class combin *CMBgen = new class combin(beta,T);      // generate weight-t vector  
  class combin *CMBsel = new class combin(NumVect,B2p); // select B2p vectors
  class MetricCalc *MC = new class MetricCalc(beta,Pi,Pd,Ps);
  unsigned char *Ur = new unsigned char [beta];
  unsigned char *sel= new unsigned char [NumVect];
  unsigned char *Lins;
  std::list<unsigned char *> LS;
  
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

  // select
  pmin=1.0;
  LS.clear();
  do {
    CMBsel->getV(sel);
    SelectVect(V,U,sel);
    prob = CalcMaxPrb(V,MC);
    if(prob<pmin){
      pmin = prob;
      for(auto it=LS.begin(); it!=LS.end(); ++it) delete [] (*it);
      LS.clear();
    } // if prob>pmin
    if(prob==pmin){
      Lins = new unsigned char [NumVect];
      memcpy(Lins,sel,sizeof(unsigned char)*NumVect);
      LS.push_back(Lins);
    } // if prob==pmin
    
    //PrintVect(sel,NumVect,""," ");
    //printf("prob=%e pmin=%e %lu\n",prob,pmin,LS.size());
    //V->print();
    //PrintList(LS,NumVect);
  } while(CMBsel->nextV());
  
  printf("LS\n");
  PrintList(LS,NumVect);
  printf("MinProb=%e ListSize=%lu\n",pmin,LS.size());
  WriteBaseFile(fn,LS,beta,B,pmin,U,MC);
  
  delete [] Ur;
  delete [] sel;
  delete U;
  delete V;
  delete CMBgen;
  delete CMBsel;
  delete MC;
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
double CalcMaxPrb(class bmatrix *V, class MetricCalc *MC){
  int    beta= V->getN();
  int    B2p = V->getM();
  double prob,pmax=0.0;
  unsigned char *X0 = new unsigned char [beta];
  unsigned char *X1 = new unsigned char [beta];
  
  // calc
  for(int i=0;i<B2p;i++){
    V->getRVunp(X0,i);
    //PrintVect(X0,beta,"","\n");
    for(int j=i+1;j<B2p;j++){
      V->getRVunp(X1,j);
      prob = MC->calc1(X0,X1);
      if(prob>pmax) pmax=prob;
      //PrintVect(X1,beta," ","");
      //printf(" [%e]\n",prob);
    } // for j
  } // for i
  
  delete [] X0;
  delete [] X1;
  return pmax;
}

//============================================================================
void PrintList(std::list<unsigned char*> LS, int len){
  for(auto it=LS.begin(); it!=LS.end(); ++it){
    for(int i=0;i<len;i++) printf("%u",(*it)[i]);
    printf("\n");
  } // for it
}

//============================================================================
void WriteBaseFile(const char *fn, std::list<unsigned char*> LS, int beta, int B, double pmin,
		   class bmatrix *U, class MetricCalc *MC){
  int nu  = LS.size();
  int B2p = (int)pow(2,B);
  unsigned int x;
  unsigned char *X = new unsigned char [beta];
  FILE *fp;
  class bmatrix *V = new class bmatrix(B2p,beta);

  if((fp=fopen(fn,"w"))==NULL){
    fprintf(stderr,"Cannot open %s\n",fn);
    exit(1);
  }
  // L1: beta nu
  fprintf(fp,"%d %d\n",beta,nu);
  // L2: b0 b1 ...
  for(int i=0;i<nu;i++) fprintf(fp,"%d ",beta+B);
  fprintf(fp,"\n");

  // L3+: cw0 cw1 ...
  for(auto it=LS.begin();it!=LS.end();++it){
    SelectVect(V,U,(*it));
    assert(CalcMaxPrb(V,MC)==pmin);
    for(int i=0;i<B2p;i++){
      V->getRVunp(X,i);
      x = MC->ConvVectUint(X,beta);
      fprintf(fp,"%d ",x);
    } // for i
    fprintf(fp,"\n");
  } // for it

  // last
  fprintf(fp,"# pmin=%e (Pi,Pd,Ps)=(%e,%e,%e)\n",pmin,Pi,Pd,Ps);
  
  fclose(fp);
  delete [] X;
  delete V;
}
