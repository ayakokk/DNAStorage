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
void InnerCodebook::VectInv4(unsigned char *VI, const unsigned char *V, int len){
  unsigned char *upper = new unsigned char[len];
  unsigned char *upper_inv = new unsigned char[len];

  // 1. 上位ビット抽出
  Extract4toUpper(upper, V, len);
  
  // 2. 上位ビットを2元反転
  VectInv(upper_inv, upper, len);
  
  // 3. 下位ビット保持 + 上位ビット反転で4元結合
  for(int i=0;i<len;i++){
    unsigned char lower = V[i] & 0x1;        // 下位ビット保持
    VI[i] = (upper_inv[i] << 1) | lower;     // 上位反転ビット + 下位ビット
  } // for i
  
  delete[] upper;
  delete[] upper_inv;
}

//================================================================================
void InnerCodebook::Extract4toUpper(unsigned char *upper, const unsigned char *V4, int len){
  for(int i=0;i<len;i++){
    if(V4[i] > 3) {
      printf("DEBUG: Extract4toUpper received invalid value: V4[%d]=%d, clamping to 3\n", i, V4[i]);
      upper[i] = 1; // 3のupper bit: (3 >> 1) & 0x1 = 1
    } else {
      upper[i] = (V4[i] >> 1) & 0x1; // 上位ビット抽出: 0,1→0, 2,3→1
    }
  } // for i
}

//================================================================================
long InnerCodebook::VectToLong(const unsigned char *V, int len){
  assert(len>0);
  long val = 0;
  for(int i=0;i<len;i++){
    if(V[i]!=0 && V[i]!=1){
      printf("ERROR: VectToLong called with 4-element data: V[%d]=%d\n", i, V[i]);
      for(int j=0;j<len;j++) printf("%d ", V[j]);
      printf("\n");
      // 4元データの場合は上位ビットを使用
      val <<= 1;
      if((V[i] >> 1) & 0x1) val |= 0x1;
    } else {
      val <<= 1;
      if(V[i]==1) val |= 0x1;
    }
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
long InnerCodebook::VectToLong4(const unsigned char *V4, int len){
  unsigned char *upper = new unsigned char[len];
  Extract4toUpper(upper, V4, len);  // 4元→上位ビット抽出
  long result = VectToLong(upper, len);  // 上位ビットでVectToLong
  delete[] upper;
  return result;
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
bool InnerCodebook::IsEqual4(const unsigned char *V0, const unsigned char *V1, int len){
  unsigned char *upper0 = new unsigned char[len];
  unsigned char *upper1 = new unsigned char[len];
  
  Extract4toUpper(upper0, V0, len);
  Extract4toUpper(upper1, V1, len);
  
  bool result = IsEqual(upper0, upper1, len);  // 既存関数再利用
  
  delete[] upper0;
  delete[] upper1;
  return result;
}

//================================================================================
bool InnerCodebook::IsInvEqual4(const unsigned char *V0, const unsigned char *V1, int len){
  unsigned char *upper0 = new unsigned char[len];
  unsigned char *upper1 = new unsigned char[len];
  
  Extract4toUpper(upper0, V0, len);
  Extract4toUpper(upper1, V1, len);
  
  bool result = IsInvEqual(upper0, upper1, len);  // 既存関数再利用
  
  delete[] upper0;
  delete[] upper1;
  return result;
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
int InnerCodebook::Balance014(const unsigned char *V4, int len){
  unsigned char *upper = new unsigned char[len];
  Extract4toUpper(upper, V4, len);  // 4元→上位ビット抽出
  int result = Balance01(upper, len);  // 既存関数再利用
  delete[] upper;
  return result;
}

//================================================================================
int InnerCodebook::BalanceSW(const unsigned char *V, int len, int ws){
  assert(ws>0 && len>=ws);
  int b, bmax=0;
  for(int i=0;i<len-ws+1;i++){
    b = abs( Balance01(&V[i],ws) );
    if(b>bmax) bmax=b;
  } // for i
  return bmax;
}

//================================================================================
int InnerCodebook::BalanceSW4(const unsigned char *V4, int len, int ws){
  assert(ws>0 && len>=ws);
  int b, bmax=0;
  for(int i=0;i<len-ws+1;i++){
    b = abs( Balance014(&V4[i],ws) );
    if(b>bmax) bmax=b;
  } // for i
  return bmax;
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
int InnerCodebook::MaxRunLen4(const unsigned char *V4, int len){
  unsigned char *upper = new unsigned char[len];
  Extract4toUpper(upper, V4, len);  // 4元→上位ビット抽出
  int result = MaxRunLen(upper, len);  // 既存関数再利用
  delete[] upper;
  return result;
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
int InnerCodebook::RunLenF4(const unsigned char *V4, int pos, int len){
  unsigned char *upper = new unsigned char[len];
  Extract4toUpper(upper, V4, len);  // 4元→上位ビット抽出
  int result = RunLenF(upper, pos, len);  // 既存関数再利用
  delete[] upper;
  return result;
}

//================================================================================
int InnerCodebook::RunLenB4(const unsigned char *V4, int pos, int len){
  unsigned char *upper = new unsigned char[len];
  Extract4toUpper(upper, V4, len);  // 4元→上位ビット抽出
  int result = RunLenB(upper, pos, len);  // 既存関数再利用
  delete[] upper;
  return result;
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
      assert(CW[i][j]>=0 && CW[i][j]<=3);  // 4元コードブック対応
    } // for j
    RLmax[i]   = MaxRunLen4(CW[i],   Nu);
    RLleft[i]  = RunLenF4(CW[i],0,   Nu);
    RLright[i] = RunLenB4(CW[i],Nu-1,Nu);
    SetWLR4(CW[i],Wleft[i],Wright[i],Nu);
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
      if(IsEqual4(CW[i],CW[j],Nu)){
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
      if(IsInvEqual4(CW[i],CW[j],Nu)){
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
    if(Balance014(CW[i],Nu)!=0){
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
void InnerCodebook::SetCWI(){
  long val;
  for(long i=0;i<Nu2p;i++) CWI[i]=-1;
  for(long i=0;i<numCW;i++){
    val = VectToLong(CW[i],Nu);
    assert(CWI[val]==-1);
    CWI[val]=i;
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
void InnerCodebook::SetWLR4(const unsigned char *V4, int *WL, int *WR, int len){
  unsigned char *upper = new unsigned char[len];
  Extract4toUpper(upper, V4, len);  // 4元→上位ビット抽出
  SetWLR(upper, WL, WR, len);       // 既存関数を再利用
  delete[] upper;
}

//================================================================================
void InnerCodebook::SetRstSymb(){
  bool f;
  Rst01=-1;
  Rst10=-1;
  for(int i=0;i<numCW;i++){
    f = true;
    for(int j=0;j<Nu-1;j++){
      if(CW[i][j]==CW[i][j+1]){
	f = false;
	break;
      } // if
    } // for j
    if(f){
      if(CW[i][0]==0) Rst01 = i;  // 0101...
      else            Rst10 = i;  // 1010...
    } // if
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
  //printf("## idx=%d pos=%d posL=%d posR=%d rlen=%d\n",idx,pos,posL,posR,rlen);
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
    //printf("## idx=%d pos=%d posL=%d lb=%d\n",idx,pos,posL,lb);
    if(abs(lb)>Delta) return false;
  } // for pos
  return true;
}

//================================================================================
bool InnerCodebook::EncodeRunLen4(const unsigned char *V4, int idx, int len){
  unsigned char *upper = new unsigned char[len*Nu];
  Extract4toUpper(upper, V4, len*Nu);  // 全体を上位ビット変換
  bool result = EncodeRunLen(upper, idx, len);  // 既存関数を再利用
  delete[] upper;
  return result;
}

//================================================================================
bool InnerCodebook::EncodeBalance4(const unsigned char *V4, int idx, int len){
  unsigned char *upper = new unsigned char[len*Nu];
  Extract4toUpper(upper, V4, len*Nu);  // 全体を上位ビット変換
  bool result = EncodeBalance(upper, idx, len);  // 既存関数を再利用
  delete[] upper;
  return result;
}

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
  SetRstSymb();
  printf("# InnerCodebook: input %s\n",fn);
  printf("# InnerCodebook: Nu=%d Nu2p=%ld numCW=%ld FlgUnique=%d FlgInvertible=%d FlgBalanced=%d\n",
	 Nu,Nu2p,numCW,FlgUnique,FlgInvertible,FlgBalanced);
  printf("# InnerCodebook: Rho=%d (ell,Delta)=(%d,%d) L0=%d\n",Rho,ell,Delta,L0);
  printf("# InnerCodebook: Rst01=%d Rst10=%d\n",Rst01,Rst10);
  assert( Nu<(int)sizeof(long)*8 );
  assert( Nu%2==0 && Nu<=ell );
  //assert( CheckVectConv() );  // 4元コードブック: VectToLong/LongToVectは未対応
  //assert( FlgUnique );  // 4元コードブック: 上位ビット重複は正常
  assert( FlgInvertible );
  assert( FlgBalanced );
  //assert( max(RLmax,numCW)<=Rho-1 );
  assert( max(RLmax,  numCW)<=Rho   );
  assert( max(RLright,numCW)<=Rho-1 );
  assert( Rst01>=0 && Rst10>=0 );
  CWL = new int [Nu2p];
  CWI = new int [Nu2p];
  //SetCWL();  // 4元コードブック: VectToLong未対応のため無効化
  //SetCWI();  // 4元コードブック: VectToLong未対応のため無効化
  
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
  delete [] CWI;
  delete [] RLmax;
  delete [] RLleft;
  delete [] RLright;
  printf("# InnerCodebook: deleted\n");
}

//================================================================================
void InnerCodebook::Encode(unsigned char *CV, const int *IV, int Ns){
  assert(Ns>0);
  unsigned char *UV = new unsigned char [Nu*Ns];
  unsigned char *U0 = new unsigned char [Nu];
  // ----- phi0
  for(int idx=0;idx<Ns;idx++){
    assert(IV[idx]>=0 && IV[idx]<numCW);
    memcpy(&UV[idx*Nu], CW[IV[idx]], sizeof(unsigned char)*Nu); 
  } // for idx
  // ----- phi1
  for(int idx=0;idx<Ns;idx++){
    // (1) copy
    memcpy(&CV[idx*Nu], &UV[idx*Nu], sizeof(unsigned char)*Nu);
    if( EncodeRunLen4(CV,idx,Ns) && EncodeBalance4(CV,idx,Ns)) continue;
    // (2) invert
    //printf("idx=%d invert\n",idx);
    VectInv4(&CV[idx*Nu], &UV[idx*Nu], Nu);
    if( EncodeRunLen4(CV,idx,Ns) && EncodeBalance4(CV,idx,Ns)) continue;
    // (3) Rst0
    //printf("idx=%d Rst0\n",idx);
    if(UV[idx*Nu]==0) memcpy(&CV[idx*Nu], CW[Rst10], sizeof(unsigned char)*Nu);
    else              memcpy(&CV[idx*Nu], CW[Rst01], sizeof(unsigned char)*Nu);
    if( EncodeRunLen4(CV,idx,Ns) && EncodeBalance4(CV,idx,Ns)) continue;
    // (4) Rst1
    //printf("idx=%d Rst1\n",idx);
    if(UV[idx*Nu]==0) memcpy(&CV[idx*Nu], CW[Rst01], sizeof(unsigned char)*Nu);
    else              memcpy(&CV[idx*Nu], CW[Rst10], sizeof(unsigned char)*Nu);
    if( EncodeRunLen4(CV,idx,Ns) && EncodeBalance4(CV,idx,Ns)) continue;
    // encoding failure
    printf("Encoding failure: idx=%d\n",idx);
    PrintVect(UV,Nu*Ns,"","\n");
    assert(false);
  } // for idx
  delete [] UV;
  delete [] U0;
  assert( MaxRunLen4(CV,Nu*Ns) <= Rho);
  assert( BalanceSW4(CV,Nu*Ns,ell) <= Delta);
  //printf("MaxRunLen=%d BalanceSW=%d\n",MaxRunLen4(CV,Nu*Ns), BalanceSW4(CV,Nu*Ns,ell));
}

//================================================================================
int InnerCodebook::CWindex(const unsigned char *V){
  long val = VectToLong4(V,Nu);  // 4元対応: 上位ビット抽出してVectToLong
  return CWI[val];
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
  printf("# CWI\n");
  cnt = 0;
  for(long i=0;i<Nu2p;i++){
    if(CWI[i]>=0){
      if(cnt%w==0) printf("# ");
      printf("%04lX:",i);
      printf("%03d ",CWI[i]);
      if(cnt%w==w-1) printf("\n");
      cnt++;
    } // if
  } // for i
  printf("\n");
}

//================================================================================
void InnerCodebook::Get_CW(unsigned char *V, int idx){
  assert(idx>=0 && idx<numCW);
  memcpy(V,CW[idx],sizeof(unsigned char)*Nu);
}

//================================================================================
int  InnerCodebook::Get_Nu(){           return Nu;}
int  InnerCodebook::Get_Nu2p(){         return Nu2p;}
int  InnerCodebook::Get_numCW(){        return numCW;}
bool InnerCodebook::Get_FlgUnique(){    return FlgUnique;}
bool InnerCodebook::Get_FlgInvertible(){return FlgInvertible;}
