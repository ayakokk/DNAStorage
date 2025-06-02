#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "bmatrix.hpp"
#include "func.hpp"
#include "InnerCodebook.hpp"


//=================================================================================
int InnerCodebook::GetNumCW(int cb){
  assert(cb>=0 && cb<nu);
  int c=0;
  for(int i=0;i<beta4p;i++) c += CWM->getV(i,cb);
  return c;
}

//=================================================================================
void InnerCodebook::FuncTest(){
  int x;
  unsigned char *V = new unsigned char [beta];
  for(int u=0;u<beta4p;u++){
    ConvIntVect(u,V,beta);
    x = ConvVectInt(V,beta);
    assert(x==u);
    // printf("%08d: ",u);
    // PrintVect(V,beta);
    // printf("> %08d\n",x);
  } // for u
  delete [] V;
}

//=================================================================================
//=================================================================================
//=================================================================================

//=================================================================================
InnerCodebook::InnerCodebook(int _beta, int _B, int _nu){
  beta  = _beta;
  B     = _B;
  nu    = _nu;
  beta4p= (int)pow(4,beta);
  B2p   = (int)pow(2,B);
  printf("# InnerCodebook: beta=%d(%d) B=%d(%d) nu=%d\n",beta,beta4p,B,B2p,nu);
  assert(beta>0 && beta<(int)sizeof(int)*4);
  assert(B>0 && B<beta*2);
  assert(nu>0);
  CWM = new class bmatrix(beta4p,nu);
  CWM->clear();
  EncMap=NULL;
  DecMap=NULL;
  FuncTest();
}

//=================================================================================
InnerCodebook::InnerCodebook(const char *fn){
  CWM = new class bmatrix(fn,true);
  beta4p = CWM->getM();
  nu     = CWM->getN();
  B2p    = GetNumCW(0);
  B      = (int)log2(B2p);  
  beta   = (int)log2(beta4p)/2;
  printf("# InnerCodebook: %s\n",fn);
  printf("# InnerCodebook: beta=%d(%d) B=%d(%d) nu=%d\n",beta,beta4p,B,B2p,nu);
  assert((int)pow(2,B)   ==B2p);
  assert((int)pow(4,beta)==beta4p);
  assert(B>0 && B<beta*2);
  assert(nu>0);
  GenMap();
  check();
  FuncTest();
}

//=================================================================================
InnerCodebook::~InnerCodebook(){
  delete CWM;
  if(EncMap!=NULL){
    for(int i=0;i<B2p;i++) delete [] EncMap[i];
    delete [] EncMap;
  } // if EncMap
  if(DecMap!=NULL){
    for(int i=0;i<beta4p;i++) delete [] DecMap[i];
    delete [] DecMap;
  } // if DecMap
  printf("# InnerCodebook: deleted\n");
}

//=================================================================================
int InnerCodebook::Get_beta(){  return beta;  }
int InnerCodebook::Get_B(){     return B;     }
int InnerCodebook::Get_nu(){    return nu;    }
int InnerCodebook::Get_beta4p(){return beta4p;}
int InnerCodebook::Get_B2p(){   return B2p;   }

//=================================================================================
void InnerCodebook::Encode(unsigned char *ICout, const int *ICin, int Nseg){
  int c;
  for(int i=0;i<Nseg;i++){
    assert(ICin[i]>=0 && ICin[i]<B2p);
    c = EncMap[ICin[i]][i%nu];
    assert(c>=0 && c<beta4p);
    ConvIntVect(c,&ICout[i*beta],beta);
  } // for i
  // (check)
  int *IDout = new int [Nseg];
  Decode(IDout,ICout,Nseg);
  assert(HammingDist(ICin,IDout,Nseg)==0);
  delete [] IDout;
}

//=================================================================================
void InnerCodebook::Decode(int *IDout, const unsigned char *IDin, int Nseg){
  int c,d;
  for(int i=0;i<Nseg;i++){
    d = ConvVectInt(&IDin[i*beta],beta);
    c = DecMap[d][i%nu];
    assert(c>=0 && c<B2p);
    IDout[i] = c;
  } // for i
}

//=================================================================================
void InnerCodebook::GetPrior(double **Px, int Nseg){
  for(int j=0;j<Nseg;j++){
    for(int i=0;i<beta4p;i++){
      if(CWM->getV(i,j%nu)==1) Px[j][i]=(double)1.0/B2p;
      else                     Px[j][i]=0.0;
    } // for i
  } // for j
}

//=================================================================================
int InnerCodebook::CWMget(int idx, int cb){
  assert(idx>=0 && idx<beta4p);
  assert(cb>=0  && cb<nu);
  return CWM->getV(idx,cb);
}

//=================================================================================
void InnerCodebook::CWMset(int idx, int cb, int val){
  assert(idx>=0 && idx<beta4p);
  assert(cb>=0  && cb<nu);
  CWM->setV(idx,cb,val);
}

//=================================================================================
void InnerCodebook::CWMset(class bmatrix *CW, int cb){
  assert(cb>=0  && cb<nu);
  CWM->setV(CW, 0,cb, 0,0, beta4p,1);
}

//=================================================================================
void InnerCodebook::CWMwrite(const char *fn){
  CWM->write(fn);
}

//=================================================================================
void InnerCodebook::GenMap(){
  int cnt;
  
  // EncMap
  EncMap = new int * [B2p];
  for(int i=0;i<B2p;i++) EncMap[i] = new int [nu];
  for(int j=0;j<nu;j++){
    cnt = 0;
    for(int i=0;i<beta4p;i++){
      if(CWM->getV(i,j)==1){
	EncMap[cnt][j] = i;
	cnt++;
	assert(cnt<=B2p);
      } // if
    } // for i
  } // for j
  
  // DecMap
  DecMap = new int * [beta4p];
  for(int i=0;i<beta4p;i++) DecMap[i] = new int [nu];
  for(int j=0;j<nu;j++){
    cnt = 0;
    for(int i=0;i<beta4p;i++){
      if(CWM->getV(i,j)==1){
	DecMap[i][j] = cnt;
	cnt++;
	assert(cnt<=B2p);
      } else {
	DecMap[i][j] = -1;
      } // if
    } // for i
  } // for j
}

//=================================================================================
bool InnerCodebook::check(){
  int c;
  // CWM -----
  for(int j=0;j<nu;j++){
    if(GetNumCW(j)!=B2p) return false;
  } // for j
  // EncMap, DecMap -----
  for(int j=0;j<nu;j++){
    for(int i=0;i<B2p;i++){
      c = EncMap[i][j];
      if(DecMap[c][j]!=i) return false;
    } // for i
  } // for j
  printf("# InnerCodebook: checked CWM, EncMap, DecMap\n");
  return true;
}

//=================================================================================
void InnerCodebook::dump(){
  int fb,fr;
  unsigned char *V = new unsigned char [beta];
  //---
  printf("CWM | DecMAP\n");
  for(int i=0;i<beta4p;i++){
    ConvIntVect(i,V,beta);
    fb = GCbalance(V,beta);
    fr = MaxRunLength(V,beta);
    printf("%04X: ",i);
    PrintVect4(V,beta);
    printf(" [%+d %d] ",fb,fr);
    for(int j=0;j<nu;j++) printf("%d",CWM->getV(i,j));
    if(DecMap!=NULL){
      printf(" | ");
      for(int j=0;j<nu;j++) printf("%+04d ",DecMap[i][j]);
    } // if DecMap
    printf("\n");    
  } // for i
  //---
  printf("EncMAP\n");
  for(int i=0;i<B2p;i++){
    printf("%02X: ",i);
    for(int j=0;j<nu;j++) printf("%04d ",EncMap[i][j]);
    printf("\n");
  } // for i
  
}
