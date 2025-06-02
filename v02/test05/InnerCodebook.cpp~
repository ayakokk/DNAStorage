#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "InnerCodebook.hpp"

#define BSIZE 4096

//================================================================================
int InnerCodebook::argmax(const int *V, int len){
  assert(len>0);
  int pos=0;
  for(int i=1;i<len;i++){
    if(V[i]>V[pos]) pos=i;
  } // for i
  return pos;
}

//================================================================================
int InnerCodebook::argmin(const int *V, int len){
  assert(len>0);
  int pos=0;
  for(int i=1;i<len;i++){
    if(V[i]<V[pos]) pos=i;
  } // for i
  return pos;
}

//================================================================================
int InnerCodebook::max(const int *V, int len){ return V[ argmax(V,len) ]; }
int InnerCodebook::min(const int *V, int len){ return V[ argmin(V,len) ]; }

//================================================================================
void InnerCodebook::VectInv(unsigned char *VI, const unsigned char *V, int len){
  for(int i=0;i<len;i++){
    assert(V[i]==0 || V[i]==1);
    VI[i] = ( V[i]==0 )? 1 : 0;
  } // for i
}

//================================================================================
long InnerCodebook::VectToLong(const unsigned char *V, int len){
  assert(len>0);
  long val = 0;
  for(int i=0;i<len;i++){
    assert(V[i]==0 || V[i]==1);
    val <<= 1;
    if(V[i]==1) val |= 0x1;
  } // for i
  return val;
}

//================================================================================
void InnerCodebook::LongToVect(unsigned char *V, long val, int len){
  assert(len>0 && val>=0);
  long mask = 0x1 << (len-1);
  for(int i=0;i<len;i++){
    V[i] = ( (val & mask)==0 )? 0 : 1;
    mask >>= 1;
  } // for i
}

//================================================================================
bool InnerCodebook::CheckVectConv(){
  bool ret = true;
  long val;
  unsigned char *V = new unsigned char [Nu];
  for(val=0;val<Nu2p;val++){
    LongToVect(V,val,Nu);
    // printf("%ld: ",val);
    // PrintVect(V,Nu,"","\n");
    if(VectToLong(V,Nu)!=val){
      ret = false;
      break;
    } // if
  } // for val
  delete [] V;
  return ret;
}

//================================================================================
void InnerCodebook::PrintVect(const unsigned char *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%u",V[i]);
  printf("%s",post);
}

//================================================================================
void InnerCodebook::PrintVect(const int *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d ",V[i]);
  printf("%s",post);
}

//================================================================================
bool InnerCodebook::IsEqual(const unsigned char *V0, const unsigned char *V1, int len){
  assert(len>0);
  for(int i=0;i<len;i++){
    assert(V0[i]==0 || V0[i]==1);
    assert(V1[i]==0 || V1[i]==1);
    if(V0[i]!=V1[i]) return false;
  } // for i
  return true;
}

//================================================================================
bool InnerCodebook::IsInvEqual(const unsigned char *V0, const unsigned char *V1, int len){
  assert(len>0);
  for(int i=0;i<len;i++){
    assert(V0[i]==0 || V0[i]==1);
    assert(V1[i]==0 || V1[i]==1);
    if(V0[i]==V1[i]) return false;
  } // for i
  return true;
}

//================================================================================
int InnerCodebook::Balance01(const unsigned char *V, int len){
  assert(len>0);
  int ret=0;
  for(int i=0;i<len;i++){
    assert(V[i]==0 || V[i]==1);
    if(V[i]==0) ret--;
    else        ret++;
  } // for i
  return ret;
}

//================================================================================
int InnerCodebook::MaxRunLen(const unsigned char *V, int len){
  assert(len>0);
  int L, Lmax=1;
  for(int i=0;i<len;i++){
    for(L=1;i+L<len;L++){
      if(V[i]!=V[i+L]) break;
    } // for L
    if(L>Lmax) Lmax = L;
    i += (L-1);
  } // for i
  return Lmax;
}

//================================================================================
int InnerCodebook::RunLenF(const unsigned char *V, int pos, int len){
  assert(pos>=0 && pos<len);
  int L;
  for(L=1;pos+L<len;L++){
    if(V[pos]!=V[pos+L]) break;
  } // for L
  return L;
}

//================================================================================
int InnerCodebook::RunLenB(const unsigned char *V, int pos, int len){
  assert(pos>=0 && pos<len);
  int L;
  for(L=1;pos-L>=0;L++){
    if(V[pos]!=V[pos-L]) break;
  } // for pos
  return L;
}

//================================================================================
void InnerCodebook::ReadFile(const char *fn){
  FILE *fp;
  char *buf = new char [BSIZE];
  if((fp=fopen(fn,"r"))==NULL){
    fprintf(stderr,"Cannot open %s\n",fn);
    exit(1);
  } // if
  // L1
  assert(fgets(buf,BSIZE,fp)!=NULL);
  Nu    = atoi(strtok(buf, " \t"));
  numCW = atoi(strtok(NULL," \t\n")); 
  assert(Nu>0 && numCW>0);
  RLmax   = new int [numCW]; 
  RLleft  = new int [numCW]; 
  RLright = new int [numCW]; 
  CW      = new unsigned char * [numCW];
  Wleft   = new int * [numCW];
  Wright  = new int * [numCW];
  for(long i=0;i<numCW;i++){
    CW[i]     = new unsigned char [Nu];
    Wleft[i]  = new int [Nu];
    Wright[i] = new int [Nu];
  } // for i
  // L2+
  for(long i=0;i<numCW;i++){
    assert(fgets(buf,BSIZE,fp)!=NULL);
    for(int j=0;j<Nu;j++){
      CW[i][j] = buf[j]-'0';
      assert(CW[i][j]==0 || CW[i][j]==1);
    } // for j
    RLmax[i]   = MaxRunLen(CW[i],   Nu);
    RLleft[i]  = RunLenF(CW[i],0,   Nu);
    RLright[i] = RunLenB(CW[i],Nu-1,Nu);
    SetWLR(CW[i],Wleft[i],Wright[i],Nu);
  } // for i
  fclose(fp);
  delete [] buf;
}

//================================================================================
void InnerCodebook::SetFlg(){
  bool flg;
  // FlgUnique
  FlgUnique = true;
  for(long i=0;i<numCW;i++){
    for(long j=i+1;j<numCW;j++){
      if(IsEqual(CW[i],CW[j],Nu)){
	FlgUnique = false;
	break;
      } // if
    } // for j
    if(!FlgUnique) break;
  } // for i
  // FlgInvertible
  FlgInvertible = true;
  for(long i=0;i<numCW;i++){
    flg = false;
    for(long j=0;j<numCW;j++){
      if(IsInvEqual(CW[i],CW[j],Nu)){
	flg = true;
	break;
      } // if
    } // for j
    if(!flg){
      FlgInvertible = false;
      break;
    } // if
  } // for i
  // FlgBalanced
  FlgBalanced = true;
  for(long i=0;i<numCW;i++){
    if(Balance01(CW[i],Nu)!=0){
      FlgBalanced = false;
      break;
    } // if
  } // for i
}

//================================================================================
void InnerCodebook::SetCWL(){
  long val;
  for(long i=0;i<Nu2p;i++) CWL[i]=0;
  for(long i=0;i<numCW;i++){
    val = VectToLong(CW[i],Nu);
    CWL[val]++;
  } // for i
}

//================================================================================
void InnerCodebook::SetWLR(const unsigned char *V, int *WL, int *WR, int len){
  assert(len>0);
  int s=0;
  // left
  for(int i=0;i<len;i++){
    assert(V[i]==0 || V[i]==1);
    s += V[i];
    WL[i] = s;
  } // for i
  // right
  s=0;
  for(int i=len-1;i>=0;i--){
    assert(V[i]==0 || V[i]==1);
    s += V[i];
    WR[i] = s;
  } // for i
}

//================================================================================
bool InnerCodebook::EncodeRunLen(const unsigned char *V, int idx, int len){
  assert(idx>=0 && idx<len);
  int pos  = idx*Nu;
  int posL, posR, rlen;
  for(posL=pos;posL>=0;posL--){
    if(V[posL]!=V[pos]) break;
  } // for posL
  for(posR=pos; posR<len*Nu; posR++){
    if(V[posR]!=V[pos]) break;
  } // for posR
  rlen = posR-posL-1;
  //printf("  idx=%d pos=%d posL=%d posR=%d rlen=%d\n",idx,pos,posL,posR,rlen);
  if(rlen>Rho) return false;
  return true;
}

//================================================================================
bool InnerCodebook::EncodeBalance(const unsigned char *V, int idx, int len){
  assert(idx>=0 && idx<len);
  int pos,posL,lb;
  for(pos=idx*Nu; pos<(idx+1)*Nu; pos++){
    posL = pos-ell+1;
    if(posL<0) continue;
    lb = Balance01(&V[posL],ell);
    printf("  idx=%d pos=%d posL=%d lb=%d\n",idx,pos,posL,lb);
    if(abs(lb)>Delta) return false;
  } // for pos
  return true;
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
InnerCodebook::InnerCodebook(const char *fn, int _Rho, int _ell, int _Delta){
  Rho   = _Rho;
  ell   = _ell;
  Delta = _Delta;
  ReadFile(fn);
  Nu2p  = (long)pow(2,Nu);
  L0    = (int)ceil((double)(ell-1)/Nu);
  SetFlg();
  printf("# InnerCodebook: input %s\n",fn);
  printf("# InnerCodebook: Nu=%d Nu2p=%ld numCW=%ld FlgUnique=%d FlgInvertible=%d FlgBalanced=%d\n",
	 Nu,Nu2p,numCW,FlgUnique,FlgInvertible,FlgBalanced);
  printf("# InnerCodebook: Rho=%d (ell,Delta)=(%d,%d) L0=%d\n",Rho,ell,Delta,L0);
  assert( Nu < (int)sizeof(long)*8 );
  assert( Nu <= ell );
  assert( CheckVectConv() );
  assert( FlgInvertible );
  assert( FlgBalanced );
  assert( max(RLmax,numCW) <= Rho );
  CWL = new int [Nu2p];
  SetCWL();
  
  //(dbg)
  PrintCodebook();
}

//================================================================================
InnerCodebook::~InnerCodebook(){
  for(long i=0;i<numCW;i++){
    delete [] CW[i];
    delete [] Wleft[i];
    delete [] Wright[i];
  } // for i
  delete [] CW;
  delete [] Wleft;
  delete [] Wright;
  delete [] CWL;
  delete [] RLmax;
  delete [] RLleft;
  delete [] RLright;
  printf("# InnerCodebook: deleted\n");
}

//================================================================================
bool InnerCodebook::Encode(unsigned char *CV, const int *IV, int Ns){
  assert(Ns>0);
  bool ret = true;
  bool f0, f1;
  unsigned char *U0 = new unsigned char [Nu];
  unsigned char *U1 = new unsigned char [Nu];
  for(int idx=0;idx<Ns;idx++){
    assert(IV[idx]>=0 && IV[idx]<numCW);
    memcpy(U0, CW[IV[idx]], sizeof(unsigned char)*Nu);
    VectInv(U1,U0,Nu);
    //--- no-invert
    memcpy(&CV[idx*Nu], U0, sizeof(unsigned char)*Nu);
    f0 =  EncodeRunLen( CV,idx,Ns);
    f0 &= EncodeBalance(CV,idx,Ns);
    //--- invert
    memcpy(&CV[idx*Nu], U1, sizeof(unsigned char)*Nu);
    f1 =  EncodeRunLen( CV,idx,Ns);
    f1 &= EncodeBalance(CV,idx,Ns);
    //--- select
    if(f0){
      memcpy(&CV[idx*Nu], U0, sizeof(unsigned char)*Nu);
    } else if(f1){
      memcpy(&CV[idx*Nu], U1, sizeof(unsigned char)*Nu);
    } else {
      memcpy(&CV[idx*Nu], U0, sizeof(unsigned char)*Nu); //TMP
      ret=false;
    } // if f0
    //(dbg)
    printf("IV[%03d]: %03d-> ",idx,IV[idx]);
    PrintVect(U0,Nu,""," ");
    PrintVect(U1,Nu,""," ");
    printf("f0=%d f1=%d\n",f0,f1);
  } // for idx
  delete [] U0;
  delete [] U1;
  return ret;
}

//================================================================================
void InnerCodebook::PrintCodebook(){
  int w=10, cnt;
  printf("# CW\n");
  for(long i=0;i<numCW;i++){
    printf("# %03ld:",i);
    PrintVect(CW[i],Nu,""," ");
    printf("(%d,L%d,R%d) ",RLmax[i],RLleft[i],RLright[i]);
    PrintVect(Wleft[i], Nu,"WL:"," ");
    PrintVect(Wright[i],Nu,"WR:"," ");
    printf("\n");
  } // for i
  printf("# CWL\n");
  cnt = 0;
  for(long i=0;i<Nu2p;i++){
    if(CWL[i]>0){
      if(cnt%w==0) printf("# ");
      printf("%04lX:",i);
      printf("%d ",CWL[i]);
      if(cnt%w==w-1) printf("\n");
      cnt++;
    } // if
  } // for i
  printf("\n");
}

//================================================================================
int  InnerCodebook::Get_Nu(){           return Nu;}
int  InnerCodebook::Get_Nu2p(){         return Nu2p;}
int  InnerCodebook::Get_numCW(){        return numCW;}
bool InnerCodebook::Get_FlgUnique(){    return FlgUnique;}
bool InnerCodebook::Get_FlgInvertible(){return FlgInvertible;}
