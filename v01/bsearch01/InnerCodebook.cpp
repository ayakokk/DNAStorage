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
void InnerCodebook::PrintParam(){
  printf("# InnerCodebook: beta=%d(%d) nu=%d\n",beta,beta4p,nu);
  printf("# InnerCodebook: B(B2p)=");
  for(int i=0;i<nu;i++) printf("%d(%d) ",B[i],B2p[i]);
  printf("\n");
}

//=================================================================================
//=================================================================================
//=================================================================================

//=================================================================================
InnerCodebook::InnerCodebook(int _beta, int *_B, int _nu){
  beta  = _beta;
  nu    = _nu;
  beta4p= (int)pow(4,beta);
  B   = new int [nu];
  B2p = new int [nu];
  memcpy(B,_B,sizeof(int)*nu);
  for(int i=0;i<nu;i++){
    B2p[i] = (int)pow(2,B[i]);
    assert(B[i]>=0 && B[i]<=beta*2);
  } // for i
  PrintParam();  
  assert(beta>0 && beta<(int)sizeof(int)*4);
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
  beta   = (int)log2(beta4p)/2;
  assert((int)pow(4,beta)==beta4p);
  //---
  nu     = CWM->getN();
  assert(nu>0);
  //---
  B   = new int [nu];
  B2p = new int [nu];
  for(int i=0;i<nu;i++){
    B2p[i] = GetNumCW(i);
    assert(B2p[i]>0);
    B[i] = (int)log2(B2p[i]);
    assert((int)pow(2,B[i])==B2p[i]);
    assert(B[i]<=beta*2);
  } // for i
  PrintParam();
  //---
  GenMap();
  assert(check());
  FuncTest();
}

//=================================================================================
InnerCodebook::InnerCodebook(class InnerCodebook *ICB){
  beta  = ICB->Get_beta();
  nu    = ICB->Get_nu();
  beta4p= ICB->Get_beta4p();
  B   = new int [nu];
  B2p = new int [nu];
  CWM = new class bmatrix(beta4p,nu);
  ICB->Get_B(B);
  ICB->Get_B2p(B2p);
  CWM->setV(ICB->CWM);
  //GenMap();
  //assert(check());
  FuncTest();
}

//=================================================================================
InnerCodebook::~InnerCodebook(){
  delete CWM;
  delete [] B;
  delete [] B2p;
  if(EncMap!=NULL){
    for(int i=0;i<nu;i++) delete [] EncMap[i];
    delete [] EncMap;
  } // if EncMap
  if(DecMap!=NULL){
    for(int i=0;i<nu;i++) delete [] DecMap[i];
    delete [] DecMap;
  } // if DecMap
  //printf("# InnerCodebook: deleted\n");
}

//=================================================================================
int  InnerCodebook::Get_beta(){  return beta;  }
int  InnerCodebook::Get_nu(){    return nu;    }
int  InnerCodebook::Get_beta4p(){return beta4p;}
//=================================================================================
void InnerCodebook::Get_B(int *_B){    memcpy(_B,B,sizeof(int)*nu);     }
void InnerCodebook::Get_B2p(int *_B2p){memcpy(_B2p,B2p,sizeof(int)*nu); }

//=================================================================================
int InnerCodebook::GetNumCW(int cb){
  assert(cb>=0 && cb<nu);
  int cnt=0;
  for(int i=0;i<beta4p;i++)
    cnt+=CWM->getV(i,cb);
  return cnt;
}

//=================================================================================
void InnerCodebook::GetNumCW(int *NumCW){
  for(int cb=0;cb<nu;cb++)
    NumCW[cb] = GetNumCW(cb);
}

//=================================================================================
void InnerCodebook::Encode(unsigned char *ICout, const int *ICin, int Nseg){
  int c;
  for(int i=0;i<Nseg;i++){
    assert(ICin[i]>=0 && ICin[i]<B2p[i%nu]);
    c = EncMap[i%nu][ICin[i]];
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
    c = DecMap[i%nu][d];
    assert(c>=0 && c<B2p[i%nu]);
    IDout[i] = c;
  } // for i
}

//=================================================================================
void InnerCodebook::GenRndInfo(int *ICin, int Nseg){
  for(int i=0;i<Nseg;i++) ICin[i] = (int)(random()%B2p[i%nu]);
}

//=================================================================================
void InnerCodebook::GetPrior(double **Px, int Nseg){
  for(int j=0;j<Nseg;j++){
    for(int i=0;i<beta4p;i++){
      if(CWM->getV(i,j%nu)==1) Px[j][i]=(double)1.0/B2p[j%nu];
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
void InnerCodebook::CWMclear(){
  CWM->clear();
}

//=================================================================================
void InnerCodebook::CWMset(){
  for(int i=0;i<beta4p;i++)
    for(int j=0;j<nu;j++) CWM->setV(i,j,1);
}

//=================================================================================
void InnerCodebook::CWMwrite(const char *fn){
  CWM->write(fn);
}

//=================================================================================
void InnerCodebook::GenMap(){
  int cnt;
  
  // EncMap
  EncMap = new int * [nu];
  for(int i=0;i<nu;i++) EncMap[i] = new int [B2p[i]];
  for(int i=0;i<nu;i++){
    cnt = 0;
    for(int j=0;j<beta4p;j++){
      if(CWM->getV(j,i)==1){
	EncMap[i][cnt] = j;
	cnt++;
	assert(cnt<=B2p[i]);
      } // if
    } // for j
  } // for i
  
  // DecMap
  DecMap = new int * [nu];
  for(int i=0;i<nu;i++) DecMap[i] = new int [beta4p];
  for(int i=0;i<nu;i++){
    cnt = 0;
    for(int j=0;j<beta4p;j++){
      if(CWM->getV(j,i)==1){
	DecMap[i][j] = cnt;
	cnt++;
	assert(cnt<=B2p[i]);
      } else {
	DecMap[i][j] = -1;
      } // if
    } // for j
  } // for i
}

//=================================================================================
bool InnerCodebook::check(){
  int c;
  // CWM -----
  for(int j=0;j<nu;j++){
    if(GetNumCW(j)!=B2p[j]) return false;
  } // for j
  // EncMap, DecMap -----
  for(int i=0;i<nu;i++){
    for(int j=0;j<B2p[i%nu];j++){
      c = EncMap[i][j];
      //printf("%d %d %d\n",i,j,c);
      if(DecMap[i][c]!=j) return false;
    } // for j
  } // for i
  //printf("# InnerCodebook: checked CWM, EncMap, DecMap\n");
  return true;
}

//=================================================================================
void InnerCodebook::PrintCWM(){
  int fb,fr;
  unsigned char *V = new unsigned char [beta];
  for(int i=0;i<beta4p;i++){
    ConvIntVect(i,V,beta);
    fb = GCbalance(V,beta);
    fr = MaxRunLength(V,beta);
    printf("%04X ",i);
    PrintVect4(V,beta);
    printf(" [%+d %d] ",fb,fr);
    for(int j=0;j<nu;j++) printf("%d",CWM->getV(i,j));
    printf("\n");
  } // for i
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
      for(int j=0;j<nu;j++) printf("%+04d ",DecMap[j][i]);
    } // if DecMap
    printf("\n");
  } // for i
  //---
  printf("EncMAP\n");
  for(int i=0;i<nu;i++){
    printf("%d: ",i);
    for(int j=0;j<B2p[i];j++) printf("%02X ",EncMap[i][j]);
    printf("\n");
  } // for i
  
}
