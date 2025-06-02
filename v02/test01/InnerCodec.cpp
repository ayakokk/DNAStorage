#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "InnerCodec.hpp"

//================================================================================
void InnerCodec::ConvLongStat(int *Qx, long val, int len){
  assert(len>=0);
  for(int i=0;i<len;i++){
    Qx[i] = val%4;
    val /= 4;
  } // for i
}

//================================================================================
long InnerCodec::ConvStatLong(const int *Qx, int len){
  assert(len>=0);
  long val=0;
  for(int i=len-1;i>=0;i--){
    assert(Qx[i]>=0 && Qx[i]<=3);
    val *= 4;
    val += Qx[i];
  } // for i
  return val;
}

//================================================================================
void InnerCodec::ConvCheck(int len){
  assert(len>=0 && (len*2<(int)sizeof(long)*8));
  long N = (long)pow(4,len);
  long val;
  int *Qx = new int [len];
  for(long i=0;i<N;i++){
    ConvLongStat(Qx,i,len);
    val = ConvStatLong(Qx,len);
    //printf("%06ld > ",i); PrintVect(S,len,""," > "); printf("%06ld\n",val);
    assert(val==i);
  } // for i
  delete [] Qx;
}

//================================================================================
int InnerCodec::MaxRunLen(const int *Qx, int len){
  int RL, maxRL = 1;
  assert(len>0);
  for(int i=0;i<len;i++){
    assert(Qx[i]>=0 && Qx[i]<=3);
    RL = 0;
    while(i+RL<len){
      if(Qx[i+RL]!=Qx[i]) break;
      RL++;
    } // while
    if(RL>maxRL) maxRL = RL;
  } // for i
  return maxRL;
}

//================================================================================
int InnerCodec::LocalBalance(const int *Qx, int len){
  int cnt01=0, cnt23=0;
  assert(len>0);
  for(int i=0;i<len;i++){
    assert(Qx[i]>=0 && Qx[i]<=3);
    if(Qx[i]==0 || Qx[i]==1) cnt01++;
    else                     cnt23++;
  } // for i
  return cnt23-cnt01;
}

//================================================================================
long InnerCodec::NextState(long Q, int a){
  assert(Q>=0 && Q<numQ);
  assert(a>=0 && a<=3);
  int *Qx  = new int [ell];
  int *Q2x = new int [ell];
  long Q2;
  ConvLongStat(Qx,Q,ell);
  for(int i=0;i<ell-1;i++) Q2x[i]=Qx[i+1];
  Q2x[ell-1] = a;
  Q2 = ConvStatLong(Q2x,ell);
  //(dbg)
  // PrintVect(Qx, ell,"Qx:", " ");
  // printf("(%d) ",a);
  // PrintVect(Q2x,ell,"Q2x:","\n");
  //---
  delete [] Qx;
  delete [] Q2x;
  return Q2;
}

//================================================================================
void InnerCodec::SetTable(){
  int  *Qx = new int [ell];
  long Q,Q2;
  int RL,LB;
  // ----- Vflg
  SetVect(Vflg,numQ,true);
  for(Q=0;Q<numQ;Q++){
    ConvLongStat(Qx,Q,ell);
    // Run-length
    RL = MaxRunLen(Qx,ell);
    if(RL>rlk) Vflg[Q] = false;
    // GC-balance
    LB = LocalBalance(Qx,ell);
    if(abs(LB)>eps) Vflg[Q] = false;
    // (dbg)
    printf("# %06ld F%d RL=%d LB=%+d ",Q,Vflg[Q],RL,LB);
    PrintVect(Qx,ell,"","\n");
  } // for Q
  // ----- STT
  for(Q=0;Q<numQ;Q++){
    ConvLongStat(Qx,Q,ell);
    if(!Vflg[Q]){
      SetVect(STT[Q],4,-1);
    } else {
      for(int a=0;a<=3;a++){
	Q2 = NextState(Q,a);
	assert(Q2>=0 && Q2<numQ);
	STT[Q][a] = (Vflg[Q2])? Q2 : -1;
      } // for a
    } // if !Vflg
    // (dbg)
    if(Vflg[Q]){
      printf("# %06ld ",Q);
      PrintVect(Qx,ell,""," ");
      PrintVect(STT[Q],4,"STT:","\n");
    } // if
  } // for Q  
  
  delete [] Qx;
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
void InnerCodec::PrintVect(const int *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d ",V[i]);
  printf("%s",post);
}

//================================================================================
void InnerCodec::SetVect(bool *V, long len, bool val){
  for(long i=0;i<len;i++) V[i]=val;
}

//================================================================================
void InnerCodec::SetVect(int *V, long len, int val){
  for(long i=0;i<len;i++) V[i]=val;
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
InnerCodec::InnerCodec(int _ell, int _rlk, int _eps){
  ell = _ell;
  rlk = _rlk;
  eps = _eps;
  numQ= (long)pow(4,ell);
  printf("# InnerCodec: ell(win size):%d  rlk(max RL):%d  eps(balance):%d  numQ:%ld\n",ell,rlk,eps,numQ);
  printf("# InnerCodec: sizeof(long):%lu\n",sizeof(long));
  assert(ell>0);
  assert(rlk>0 && rlk<ell);
  assert(eps>=0 && eps<ell);
  assert(ell*2<(int)sizeof(long)*8);
  ConvCheck(ell);

  Vflg = new bool [numQ];
  STT = new int * [numQ];
  for(long i=0;i<numQ;i++) STT[i] = new int [4];
  SetTable();
}

//================================================================================
InnerCodec::~InnerCodec(){
  for(long i=0;i<numQ;i++) delete [] STT[i];
  delete [] STT;
  delete [] Vflg;
  printf("# InnerCodec: deleted\n");
}
