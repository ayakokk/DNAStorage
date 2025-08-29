#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <string>

#include "InnerCodebook.hpp"
#include "ChannelMatrix.hpp"
#include "IDSchannel.hpp"
#include "SLFBAdec.hpp"

#define NuMax 10

//================================================================================
int SLFBAdec::max(int a, int b){return (a>b)? a : b;}
int SLFBAdec::min(int a, int b){return (a<b)? a : b;}

//================================================================================
double SLFBAdec::Psub(unsigned char a, unsigned char b){
  assert( a==0 || a==1 );
  assert( b==0 || b==1 );
  return (a==b)? 1.0-Ps : Ps;
}

//================================================================================
void SLFBAdec::ClearMat(double **X, int M0, int N0){
  for(int i=0;i<M0;i++){
    for(int j=0;j<N0;j++) X[i][j]=0.0;
  } // for i
}

//================================================================================
void SLFBAdec::CopyMat(double **Y, const double **X, int M0, int N0){
  assert(M0>0 && N0>0);
  for(int i=0;i<M0;i++) memcpy(Y[i],X[i],sizeof(double)*N0);
}

//================================================================================
void SLFBAdec::PrintMat(const double **X, int M0, int N0, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<M0;i++){
    printf("%03d: ",i);
    for(int j=0;j<N0;j++) printf("%.2e ",X[i][j]);
    printf("\n");
  } // for i
  printf("%s",post);
}

//================================================================================
void SLFBAdec::MultMat(double **Y, const double **X0, const double **X1, int M0, int N0, int M1, int N1){
  assert(M0>0 && N0>0 && M1>0 && N1>0);
  assert(N0==M1);
  double **Y0 = new double * [M0];
  for(int i=0;i<M0;i++){
    Y0[i] = new double [N1];
    for(int j=0;j<N1;j++) Y0[i][j]=0.0;
  } // for i
  //--- calc
  for(int i=0;i<M0;i++){
    for(int j=0;j<N1;j++){
      for(int k=0;k<N0;k++){
	Y0[i][j] += X0[i][k]*X1[k][j];
      } // for k
    } // for j
  } // for i
  //--- Y <- Y0
  for(int i=0;i<M0;i++) memcpy(Y[i],Y0[i],sizeof(double)*N1);
  //--- del
  for(int i=0;i<M0;i++) delete [] Y0[i];
  delete [] Y0;
}

//================================================================================
void SLFBAdec::PrintVect(const unsigned char *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%u",V[i]);
  printf("%s",post);
}

//================================================================================
void SLFBAdec::PrintVect(const int *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d ",V[i]);
  printf("%s",post);
}

//================================================================================
void SLFBAdec::PrintVect(const double *V, int len, const char *pre, const char *post){
  int w=20;
  printf("%s",pre);
  for(int i=0;i<len;i++){
    printf("%.2e ",V[i]);
    if(i%w==w-1) printf("\n");
  } // for i
  printf("%s",post);
}

//================================================================================
long SLFBAdec::VectToLong(const unsigned char *V, int len){
  assert(len>=0);
  long val = 0;
  for(int i=0;i<len;i++){
    assert(V[i]==0 || V[i]==1);
    val <<= 1;
    if(V[i]==1) val |= 0x1;
  } // for i
  return val;
}

//================================================================================
void SLFBAdec::LongToVect(unsigned char *V, long val, int len){
  assert(len>0 && val>=0);
  long mask = 0x1 << (len-1);
  for(int i=0;i<len;i++){
    V[i] = ( (val & mask)==0 )? 0 : 1;
    mask >>= 1;
  } // for i
}

//--------------------------------------------------------------------------------
void SLFBAdec::normalize(double *V, int len){
  assert(len>0);
  double s=0;
  for(int i=0;i<len;i++){
    assert( V[i]>=0.0 );
    s += V[i];
  } // for i
  assert( s>0 );
  for(int i=0;i<len;i++) V[i] /= s;
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

//================================================================================
void SLFBAdec::SetGD(){
  double **GDX = new double * [Drng];  // for calc
  GD = new double * [Drng];
  for(int i=0;i<Drng;i++){
    GD[i]  = new double [Drng];
    GDX[i] = new double [Drng];
  } // for i
  // --- GDX (1-bit)
  ClearMat(GDX,Drng,Drng);
  GDX[0][0] = 1.0-Pi;
  GDX[Drng-1][Drng-1] = 1.0-Pd;
  for(int i=0;i<Drng-1;i++) GDX[i][i+1] = Pi;
  for(int i=1;i<Drng;  i++) GDX[i][i-1] = Pd;
  for(int i=1;i<Drng-1;i++) GDX[i][i]   = 1.0-(Pi+Pd);
  // --- GD <- identity mat
  ClearMat(GD,Drng,Drng);
  for(int i=0;i<Drng;i++) GD[i][i] = 1.0;
  // --- calc
  for(int i=0;i<Nu;i++) MultMat(GD,(const double**)GD,(const double **)GDX,Drng,Drng,Drng,Drng);
  // (dbg)
  // PrintMat((const double **)GDX,Drng,Drng,"GDX\n","");
  // PrintMat((const double **)GD, Drng,Drng,"GD\n", "");
  // --- del
  for(int i=0;i<Drng;i++) delete [] GDX[i];
  delete [] GDX;
}

//================================================================================
void SLFBAdec::DelGD(){
  for(int i=0;i<Drng;i++) delete [] GD[i];
  delete [] GD;
}

//================================================================================
void SLFBAdec::SetGX(){
  long ly2p,x;
  unsigned char *X = new unsigned char [Nu];
  GX = new double ** [Nu*2+1];
  for(int ly=0;ly<=Nu*2;ly++){
    if(ly<Nu2min || ly>Nu2max){
      // approximate
      GX[ly]    = new double * [1];
      GX[ly][0] = new double [1];
      GX[ly][0][0] = (ly<Nu)? pow(Pd,Nu-ly) : pow(Pi,ly-Nu); //?
      //printf("GX[%d]=%e\n",ly,GX[ly][0][0]);
    } else {
      // exact
      ly2p = (long)pow(2,ly);
      GX[ly] = new double * [ly2p];
      for(long y=0;y<ly2p;y++){
	GX[ly][y] = new double [Q];
	for(int xi=0;xi<Q;xi++){
	  ICB->Get_CW(X,xi);
	  x = VectToLong(X,Nu);
	  GX[ly][y][xi] = CalcPyx(y,x,ly,Nu);
	  //(dbg)
	  //printf("%e ly=%d y=%ld xi=%d x=%ld ",GX[ly][y][xi],ly,y,xi,x);
	  //PrintVect(X,Nu,"X=","\n");
	} // for x
      } // for y
    } // if ly
  } // for ly
  delete [] X;
}

//================================================================================
void SLFBAdec::DelGX(){
  int i2p;
  for(int i=0;i<=Nu*2;i++){
    if(i<Nu2min || i>Nu2max){
      delete [] GX[i][0];
      delete [] GX[i];
    } else {
      i2p = (long)pow(2,i);
      for(long y=0;y<i2p;y++) delete [] GX[i][y];
      delete [] GX[i];
    } // if i
  } // for i
  delete [] GX;
}

//================================================================================
double SLFBAdec::CalcPyx(long y, long x, int ly, int lx){
  assert(lx>0 && ly>=0);
  //printf("y=%ld x=%ld ly=%d lx=%d\n",y,x,ly,lx);
  if(ly==0) return pow(Pd,lx);
  long   x1,y1;
  double ret,qt,qi,qd;
  unsigned char *X = new unsigned char [lx];
  unsigned char *Y = new unsigned char [ly];
  LongToVect(X,x,lx);
  LongToVect(Y,y,ly);
  if(lx==1){
    if(ly==1){
      ret = Psub(X[0],Y[0]) * Pt;
    } else if(ly==2) {
      ret = Psub(X[0],Y[0]) * Psub(X[0],Y[1]) * Pi;
    } else {
      assert(ly>=3); 
      ret = 0.0;
    } // if ly
  } else {
    x1  = VectToLong(&X[1],lx-1);
    // -----trans
    qt = Psub(X[0],Y[0]) * Pt;
    y1 = (ly==1)? 0 : VectToLong(&Y[1],ly-1);
    qt *= CalcPyx(y1,x1,ly-1,lx-1);
    // -----ins
    if(ly<2){
      qi = 0.0;
    } else {
      qi = Psub(X[0],Y[0]) * Psub(X[0],Y[1]) * Pi;
      y1 = (ly==2)? 0 : VectToLong(&Y[2],ly-2);
      qi *= CalcPyx(y1,x1,ly-2,lx-1);	
    } // if ly<2
    // -----del
    qd = Pd;
    y1 = y;
    qd *= CalcPyx(y1,x1,ly,lx-1);
    // ---
    ret = qt + qi + qd;
  } // if lx
  delete [] X;
  delete [] Y;
  //printf("ret=%e\n",ret);
  return ret;
}

//================================================================================
double SLFBAdec::GetGX(int Nu2, long y, long xi){
  assert(Nu2>=0 && Nu2<=2*Nu);
  assert(y >=0);
  assert(xi>=0 && xi<Q);
  if(Nu2<Nu2min || Nu2>Nu2max) return GX[Nu2][0][0 ]; // approx
  else                         return GX[Nu2][y][xi]; 
}

//================================================================================
void SLFBAdec::SetFG(){
  PF = new double * [Ns+1];
  PB = new double * [Ns+1];
  for(int i=0;i<Ns+1;i++){
    PF[i] = new double [Drng];
    PB[i] = new double [Drng];
  } // for i
  PU = new double * [Ns];
  PD = new double * [Ns];
  PI = new double * [Ns];
  PO = new double * [Ns];
  for(int i=0;i<Ns;i++){
    PU[i] = new double [Q];
    PD[i] = new double [Q];
    PI[i] = new double [Q];
    PO[i] = new double [Q];
  } // for i
  Yin = new unsigned char [Nb+Dmax];
}

//================================================================================
void SLFBAdec::DelFG(){
  for(int i=0;i<Ns+1;i++){
    delete [] PF[i];
    delete [] PB[i];
  } // for i
  delete [] PF;
  delete [] PB;
  for(int i=0;i<Ns;i++){
    delete [] PU[i];
    delete [] PD[i];
    delete [] PI[i];
    delete [] PO[i];
  } // for i
  delete [] PU;
  delete [] PD;
  delete [] PI;
  delete [] PO;
  delete [] Yin;
}

//================================================================================
void SLFBAdec::InitFG(const unsigned char *RW, const double **Pin, int Nb2){
  ClearMat(PF,Ns+1,Drng);
  ClearMat(PB,Ns+1,Drng);
  ClearMat(PU,Ns,Q);
  ClearMat(PD,Ns,Q);
  ClearMat(PO,Ns,Q);
  CopyMat(PI,Pin,Ns,Q);
  PF[0 ][0     -Dmin] = 1.0;
  PB[Ns][Nb2-Nb-Dmin] = 1.0;
  memcpy(Yin,RW,sizeof(unsigned char)*Nb2);
}

//================================================================================
void SLFBAdec::CalcPU(int idx){
  assert(idx>=0 && idx<Ns);
  int Vin,Vup;
  for(Vup=0;Vup<Q;Vup++){
    PU[idx][Vup] = 0.0;
    for(Vin=0;Vin<Q;Vin++){
      PU[idx][Vup] += (PI[idx][Vin] * ECM->GetPxy(Vin,Vup));
    } // for Vin
  } // for Vup
  normalize(PU[idx], Q);
}

//================================================================================
void SLFBAdec::CalcPO(int idx){
  assert(idx>=0 && idx<Ns);
  int Vout,Vdn;
  for(Vout=0;Vout<Q;Vout++){
    PO[idx][Vout] = 0.0;
    for(Vdn=0;Vdn<Q;Vdn++){
      PO[idx][Vout] += (PD[idx][Vdn] * ECM->GetPxy(Vout,Vdn));
    } // for Vdn
  } // for Vout
  normalize(PO[idx], Q);
}

//================================================================================
void SLFBAdec::CalcPF(int idx, int Nb2){
  assert(idx>=0 && idx<Ns);
  int  d0,d1;
  int  xi;
  int  Nu2, iL, iR;
  //long x;
  long y;
  double s,ss;
  //unsigned char *X = new unsigned char [Nu];
  for(d1=Dmin;d1<=Dmax;d1++){
    PF[idx+1][d1-Dmin] = 0.0;
    for(d0=Dmin;d0<=Dmax;d0++){
      Nu2 = Nu+d1-d0;
      iL  = idx    *Nu + d0;
      iR  = (idx+1)*Nu + d1 - 1;
      if( (Nu2<0) || (Nu2>2*Nu) || (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
      if( Nu2>=Nu2min && Nu2<=Nu2max ){
	y = VectToLong( &Yin[iL], Nu2 );
	s = 0.0;
	for(xi=0;xi<Q;xi++){
	  //ICB->Get_CW(X,xi);
	  //x  = VectToLong(X,Nu);
	  ss = PU[idx][xi] * GetGX(Nu2,y,xi);
	  s += ss;
	  //(dbg)
	  //printf("%e d1=%d d0=%d Nu2=%d iL=%d iR=%d y=%ld xi=%d x=%ld\n",ss,d1,d0,Nu2,iL,iR,y,xi,x);
	} // for xi
      } else {
	s = GetGX(Nu2,0,0);
      } // if Nu2
      PF[idx+1][d1-Dmin] += (s*PF[idx][d0-Dmin]);
    } // for d0
  } // for d1
  normalize(PF[idx+1],Drng);
  //delete [] X;
}

//================================================================================
void SLFBAdec::CalcPB(int idx, int Nb2){
  assert(idx>=0 && idx<Ns);
  int  d0,d1;
  int  xi;
  int  Nu2, iL, iR;
  //long x;
  long y;
  double s,ss;
  //unsigned char *X = new unsigned char [Nu];
  for(d0=Dmin;d0<=Dmax;d0++){
    PB[idx][d0-Dmin] = 0.0;
    for(d1=Dmin;d1<=Dmax;d1++){
      Nu2 = Nu+d1-d0;
      iL  = idx    *Nu + d0;
      iR  = (idx+1)*Nu + d1 - 1;
      if( (Nu2<0) || (Nu2>2*Nu) || (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
      if( Nu2>=Nu2min && Nu2<=Nu2max ){
	y = VectToLong( &Yin[iL], Nu2 );
	s = 0.0;
	for(xi=0;xi<Q;xi++){
	  //ICB->Get_CW(X,xi);
	  //x  = VectToLong(X,Nu);
	  ss = PU[idx][xi] * GetGX(Nu2,y,xi);
	  s += ss;
	  //(dbg)
	  //printf("%e d1=%d d0=%d Nu2=%d iL=%d iR=%d y=%ld xi=%d x=%ld\n",ss,d1,d0,Nu2,iL,iR,y,xi,x);
	} // for xi
      } else {
	s = GetGX(Nu2,0,0);
      } // if Nu2
      PB[idx][d0-Dmin] += (s*PB[idx+1][d1-Dmin]);
    } // for d1
  } // for d0
  normalize(PB[idx],Drng);
  //delete [] X;
}

//================================================================================
void SLFBAdec::CalcPD(int idx, int Nb2){
  assert(idx>=0 && idx<Ns);
  int  d0,d1;
  int  xi;
  int  Nu2, iL, iR;
  //long x;
  long y;
  double s,ss;
  //unsigned char *X = new unsigned char [Nu];
  for(xi=0;xi<Q;xi++){
    PD[idx][xi] = 0.0;
    //ICB->Get_CW(X,xi);
    //x  = VectToLong(X,Nu);
    for(d1=Dmin;d1<=Dmax;d1++){
      s = 0.0;
      for(d0=Dmin;d0<=Dmax;d0++){
	Nu2 = Nu+d1-d0;
	iL  = idx    *Nu + d0;
	iR  = (idx+1)*Nu + d1 - 1;
	if( (Nu2<0) || (Nu2>2*Nu) || (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
	y = VectToLong( &Yin[iL], Nu2 );
	ss = PF[idx][d0-Dmin] * GetGX(Nu2,y,xi);
	s += ss;
	//(dbg)
	//printf("%e d1=%d d0=%d Nu2=%d iL=%d iR=%d y=%ld xi=%d x=%ld\n",ss,d1,d0,Nu2,iL,iR,y,xi,x);
      } // for d0
      PD[idx][xi] += (s*PB[idx+1][d1-Dmin]);
    } // for d1
  } // for xi
  normalize(PD[idx],Q);
  //delete [] X;
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
SLFBAdec::SLFBAdec(class InnerCodebook *_ICB, class ChannelMatrix *_ECM, class IDSchannel *_CH){
  ICB  = _ICB;
  ECM  = _ECM;
  CH   = _CH;
  Nu   = ICB->Get_Nu();
  Q    = ICB->Get_numCW();
  Nb   = CH->GetN();
  Dmax = CH->GetDmax();
  Dmin = CH->GetDmin();
  Pi   = CH->GetPi();
  Pd   = CH->GetPd();
  Ps   = CH->GetPs();
  Pt   = 1.0-Pi-Pd;
  Nu2p = (long) pow(2,Nu);
  Ns   = Nb/Nu;
  Drng = Dmax-Dmin+1;
  Nu2max = Nu + (int)ceil( (double)Nu*Pi ) + 2;
  Nu2min = Nu - (int)ceil( (double)Nu*Pd ) - 2;
  Nu2max = min(Nu2max, Nu*2);
  Nu2min = max(Nu2min, 0   );
  k_mer_length = 4;  // 初期値としてk=4を設定
  printf("# SLFBAdec: Ns=%d (Nu;Nu2min,Nu2max)=(%d(%ld);%d,%d) Nb=%d Q=%d\n",Ns,Nu,Nu2p,Nu2min,Nu2max,Nb,Q);
  printf("# SLFBAdec: (Pi,Pd,Ps,Pt)=(%e,%e,%e,%e) (Dmin,Dmax:Drng)=(%d,%d:%d)\n",
	 Pi,Pd,Ps,Pt,Dmin,Dmax,Drng);
  assert(Nb%Nu==0);
  assert(Nu<=NuMax);
  assert(ECM->GetM()==Q && ECM->GetN()==Q);
  assert(Pt>0.0 && Pt<=1.0);
  //----- set tables
  printf("# SLFBAdec: generating GD & GX\n");
  SetGD();
  SetGX();
  printf("# SLFBAdec: generating factor graph\n");
  SetFG();
  printf("# SLFBAdec: generating error state extended factor graph\n");
  SetFGE();
}

//================================================================================
SLFBAdec::~SLFBAdec(){
  DelGD();
  DelGX();
  DelFG();
  DelFGE();
  printf("# SLFBAdec: deleted\n");
}

//================================================================================
void SLFBAdec::Decode(double **Pout, const unsigned char *RW, int Nb2, const int *dbgIW, const double **Pin){
  int idx;
  assert( Nb2>=Nb+Dmin && Nb2<=Nb+Dmax );
  
  // デコーダー選択: 環境変数DECODER_MODEで制御
  // DECODER_MODE=LEGACY: 従来の2Dデコーダー（ベースライン）
  // DECODER_MODE=DEC2:   新しい3Dラティスデコーダー（Dec2）
  const char* decoder_mode = getenv("DECODER_MODE");
  
  if (decoder_mode && strcmp(decoder_mode, "LEGACY") == 0) {
    printf("# Using LEGACY 2D decoder (baseline)\n");
    InitFG(RW,Pin,Nb2);
    for(idx=0;   idx<Ns;idx++) CalcPU(idx);
    for(idx=0;   idx<Ns;idx++) CalcPF(idx,Nb2);
    for(idx=Ns-1;idx>=0;idx--) CalcPB(idx,Nb2);
    for(idx=0;   idx<Ns;idx++) CalcPD(idx,Nb2);
    for(idx=0;   idx<Ns;idx++) CalcPO(idx);
  } else {
    printf("# Using Dec2 3D lattice decoder with error state memory\n");
    InitFGE(RW,Pin,Nb2);
    for(idx=0;   idx<Ns;idx++) CalcPU(idx);
    for(idx=0;   idx<Ns;idx++) CalcPFE(idx,Nb2);   // 3D前進確率
    for(idx=Ns-1;idx>=0;idx--) CalcPBE(idx,Nb2);   // 3D後進確率
    for(idx=0;   idx<Ns;idx++) CalcPDE(idx,Nb2);   // 3D事後確率
    for(idx=0;   idx<Ns;idx++) CalcPO(idx);
  }
  CopyMat(Pout,(const double**)PO,Ns,Q);
  
  // PrintNode(dbgIW);
  // PrintNode(0,   dbgIW[0]);
  // PrintNode(1,   dbgIW[1]);
  // PrintNode(2,   dbgIW[2]);
  // PrintNode(Ns-4,dbgIW[Ns-2]);
  // PrintNode(Ns-3,dbgIW[Ns-2]);
  // PrintNode(Ns-2,dbgIW[Ns-2]);
  // PrintNode(Ns-1,dbgIW[Ns-1]);
}

//================================================================================
void SLFBAdec::Decode(double **Pout, const unsigned char *RW, int Nb2, const int *dbgIW){
  double **Pin = new double * [Ns];
  for(int i=0;i<Ns;i++){
    Pin[i] = new double [Q];
    for(int x=0;x<Q;x++) Pin[i][x] = (double)1.0/Q;
  } // for i
  Decode(Pout,RW,Nb2,dbgIW,(const double **)Pin);
  for(int i=0;i<Ns;i++) delete [] Pin[i];
  delete [] Pin;
}

//================================================================================
void SLFBAdec::PrintNode(int idx, int iw){
  assert(idx>=0 && idx<Ns);
  printf("[%04d] %02d\n",idx,iw);
  PrintVect(PF[idx  ],Drng,"PF: ","\n");
  PrintVect(PB[idx+1],Drng,"PB: ","\n");
  PrintVect(PI[idx  ],Q,   "PI: ","\n");
  PrintVect(PU[idx  ],Q,   "PU: ","\n");
  PrintVect(PD[idx  ],Q,   "PD: ","\n");
  PrintVect(PO[idx  ],Q,   "PO: ","\n");
}

//================================================================================
void SLFBAdec::PrintNode(const int *dbgIW){
  for(int i=0;i<Ns;i++) PrintNode(i,dbgIW[i]);
}

//================================================================================
// テーブル出力機能の実装
//================================================================================
void SLFBAdec::exportGDTable(const char* filename) {
  printf("# Exporting GD table to %s\n", filename);
  
  FILE* file = fopen(filename, "w");
  if (!file) {
    printf("# Error: Cannot open file %s\n", filename);
    return;
  }
  
  // ヘッダー情報
  fprintf(file, "# GD (Drift Probability) Table\n");
  fprintf(file, "# Parameters: Pi=%.6f Pd=%.6f Ps=%.6f Pt=%.6f\n", Pi, Pd, Ps, Pt);
  fprintf(file, "# Drng=%d Nu=%d\n", Drng, Nu);
  fprintf(file, "# Format: d0 d1 probability\n");
  fprintf(file, "d0\td1\tprobability\n");
  
  // GDテーブル出力
  for (int d0 = 0; d0 < Drng; d0++) {
    for (int d1 = 0; d1 < Drng; d1++) {
      fprintf(file, "%d\t%d\t%.12e\n", d0 - Dmin, d1 - Dmin, GD[d0][d1]);
    }
  }
  
  fclose(file);
  printf("# GD table exported successfully\n");
}

//================================================================================
void SLFBAdec::exportGXTable(const char* filename) {
  printf("# Exporting GX table to %s\n", filename);
  
  FILE* file = fopen(filename, "w");
  if (!file) {
    printf("# Error: Cannot open file %s\n", filename);
    return;
  }
  
  // ヘッダー情報
  fprintf(file, "# GX (Channel Output Probability) Table\n");
  fprintf(file, "# Parameters: Pi=%.6f Pd=%.6f Ps=%.6f Pt=%.6f\n", Pi, Pd, Ps, Pt);
  fprintf(file, "# Nu=%d Q=%d Nu2min=%d Nu2max=%d\n", Nu, Q, Nu2min, Nu2max);
  fprintf(file, "# Format: ly binary(y) xi probability\n");
  fprintf(file, "ly\tbinary(y)\txi\tprobability\n");
  
  unsigned char *Y = new unsigned char[Nu*2];
  
  // GXテーブル出力
  for (int ly = 0; ly <= Nu*2; ly++) {
    if (ly < Nu2min || ly > Nu2max) {
      // 近似値の場合
      fprintf(file, "%d\t(approx)\t-1\t%.12e\n", ly, GX[ly][0][0]);
    } else {
      // 正確値の場合
      long ly2p = (long)pow(2, ly);
      for (long y = 0; y < ly2p; y++) {
        // yをバイナリ文字列に変換
        LongToVect(Y, y, ly);
        std::string y_binary = "";
        for (int i = 0; i < ly; i++) {
          y_binary += (char)('0' + Y[i]);
        }
        
        for (int xi = 0; xi < Q; xi++) {
          fprintf(file, "%d\t%s\t%d\t%.12e\n", ly, y_binary.c_str(), xi, GX[ly][y][xi]);
        }
      }
    }
  }
  
  delete[] Y;
  fclose(file);
  printf("# GX table exported successfully\n");
}

//================================================================================
void SLFBAdec::exportAllTables(const char* output_dir) {
  printf("# Exporting all probability tables to directory: %s\n", output_dir);
  
  // ディレクトリ作成
  char mkdir_cmd[1024];
  snprintf(mkdir_cmd, sizeof(mkdir_cmd), "mkdir -p %s", output_dir);
  system(mkdir_cmd);
  
  // ファイル名生成
  char gd_file[1024], gx_file[1024], summary_file[1024];
  snprintf(gd_file, sizeof(gd_file), "%s/GD_table.txt", output_dir);
  snprintf(gx_file, sizeof(gx_file), "%s/GX_table.txt", output_dir);
  snprintf(summary_file, sizeof(summary_file), "%s/sim01_tables_summary.txt", output_dir);
  
  // テーブル出力
  exportGDTable(gd_file);
  exportGXTable(gx_file);
  
  // サマリーファイル作成
  FILE* summary = fopen(summary_file, "w");
  if (summary) {
    fprintf(summary, "# SLFBAdec sim01 Probability Tables Summary\n");
    time_t rawtime;
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    fprintf(summary, "# Generated on: %s", asctime(timeinfo));
    fprintf(summary, "\n## Parameters\n");
    fprintf(summary, "Nu (symbol length): %d\n", Nu);
    fprintf(summary, "Q (codewords): %d\n", Q);
    fprintf(summary, "Channel parameters: Pi=%.6f Pd=%.6f Ps=%.6f Pt=%.6f\n", Pi, Pd, Ps, Pt);
    fprintf(summary, "Drift range: %d to %d (Drng=%d)\n", Dmin, Dmax, Drng);
    fprintf(summary, "Nu2 range: %d to %d\n", Nu2min, Nu2max);
    fprintf(summary, "\n## Generated Files\n");
    fprintf(summary, "1. GD_table.txt - Drift probability table P(d1|d0)\n");
    fprintf(summary, "2. GX_table.txt - Channel output probability table P(y|x)\n");
    fprintf(summary, "\n## Table Characteristics\n");
    fprintf(summary, "- GD table: %dx%d matrix\n", Drng, Drng);
    fprintf(summary, "- GX table: Variable size based on ly and xi\n");
    fprintf(summary, "- Theoretical IDS channel model\n");
    fclose(summary);
  }
  
  printf("# All probability tables exported to: %s\n", output_dir);
}

//================================================================================
// エラー状態拡張FGデータ構造の設定
//================================================================================
void SLFBAdec::SetFGE(){
  PFE = new double ** [Ns+1];
  PBE = new double ** [Ns+1];
  for(int i=0;i<Ns+1;i++){
    PFE[i] = new double * [Drng];
    PBE[i] = new double * [Drng];
    for(int d=0;d<Drng;d++){
      PFE[i][d] = new double [NUM_ERROR_STATES];
      PBE[i][d] = new double [NUM_ERROR_STATES];
    }
  }
  
  PE = new double ** [Ns];
  for(int i=0;i<Ns;i++){
    PE[i] = new double * [NUM_ERROR_STATES];
    for(int e=0;e<NUM_ERROR_STATES;e++){
      PE[i][e] = new double [NUM_ERROR_STATES];
    }
  }
  
  printf("# SLFBAdec: SetFGE completed (3D arrays allocated)\n");
}

//================================================================================
// エラー状態拡張FGデータ構造の削除
//================================================================================
void SLFBAdec::DelFGE(){
  for(int i=0;i<Ns+1;i++){
    for(int d=0;d<Drng;d++){
      delete [] PFE[i][d];
      delete [] PBE[i][d];
    }
    delete [] PFE[i];
    delete [] PBE[i];
  }
  delete [] PFE;
  delete [] PBE;
  
  for(int i=0;i<Ns;i++){
    for(int e=0;e<NUM_ERROR_STATES;e++){
      delete [] PE[i][e];
    }
    delete [] PE[i];
  }
  delete [] PE;
  
  printf("# SLFBAdec: DelFGE completed (3D arrays deallocated)\n");
}

//================================================================================
// エラー状態拡張FGの初期化
//================================================================================
void SLFBAdec::InitFGE(const unsigned char *RW, const double **Pin, int Nb2){
  // 既存の初期化を実行
  InitFG(RW, Pin, Nb2);
  
  // エラー状態確率を初期化
  for(int i=0;i<Ns+1;i++){
    for(int d=0;d<Drng;d++){
      for(int e=0;e<NUM_ERROR_STATES;e++){
        PFE[i][d][e] = 0.0;
        PBE[i][d][e] = 0.0;
      }
    }
  }
  
  // 初期エラー状態: 論文に従いMatch状態で開始 (σ_0 = (0, 0, M))
  PFE[0][0-Dmin][ERROR_MATCH] = 1.0;
  PFE[0][0-Dmin][ERROR_INSERTION] = 0.0;
  PFE[0][0-Dmin][ERROR_DELETION] = 0.0;
  PFE[0][0-Dmin][ERROR_SUBSTITUTION] = 0.0;
  
  // 最終エラー状態も均等分布で初期化
  PBE[Ns][Nb2-Nb-Dmin][ERROR_MATCH] = 0.25;
  PBE[Ns][Nb2-Nb-Dmin][ERROR_INSERTION] = 0.25;
  PBE[Ns][Nb2-Nb-Dmin][ERROR_DELETION] = 0.25;
  PBE[Ns][Nb2-Nb-Dmin][ERROR_SUBSTITUTION] = 0.25;
  
  printf("# SLFBAdec: InitFGE completed\n");
}

//================================================================================
// エラー状態遷移確率の計算（P(e_{i+1}|e_i)）
//================================================================================
void SLFBAdec::CalcPE(int idx){
  assert(idx>=0 && idx<Ns);
  
  // 段階1a: 単純なエラー状態メモリモデル
  // 実際の実装では、ここでDNArSim-mainの確率テーブルを使用予定
  
  for(int e0=0; e0<NUM_ERROR_STATES; e0++){
    for(int e1=0; e1<NUM_ERROR_STATES; e1++){
      if(e0 == e1){
        // 同一エラー状態に留まる確率（高い）
        PE[idx][e0][e1] = 0.7;
      } else {
        // 他のエラー状態に遷移する確率（低く均等）
        PE[idx][e0][e1] = 0.1;
      }
    }
  }
  
  // 正規化
  for(int e0=0; e0<NUM_ERROR_STATES; e0++){
    double sum = 0.0;
    for(int e1=0; e1<NUM_ERROR_STATES; e1++){
      sum += PE[idx][e0][e1];
    }
    for(int e1=0; e1<NUM_ERROR_STATES; e1++){
      PE[idx][e0][e1] /= sum;
    }
  }
}

//================================================================================
// エラー状態を含む前進確率計算（3Dラティス）
//================================================================================
void SLFBAdec::CalcPFE(int idx, int Nb2){
  assert(idx >= 0 && idx < Ns);

  // エラー状態遷移確率を事前計算
  CalcPE(idx);

  // 次の状態 t+1 の確率を初期化
  for(int d1 = Dmin; d1 <= Dmax; d1++){
    for(int e1 = 0; e1 < NUM_ERROR_STATES; e1++){
      PFE[idx+1][d1-Dmin][e1] = 0.0;
    }
  }
  
  // 現在の状態 t から 次の状態 t+1 への全ての遷移を計算
  for (int d0 = Dmin; d0 <= Dmax; d0++) {
    for (int e0 = 0; e0 < NUM_ERROR_STATES; e0++) {
      
      // もし現在の状態(d0, e0)に至る確率が0なら、このパスからの寄与はない
      if (PFE[idx][d0-Dmin][e0] == 0.0) continue;

      for (int d1 = Dmin; d1 <= Dmax; d1++) {
        
        // --- ここからが分岐メトリック γ_t の計算 ---
        int Nu2 = Nu + d1 - d0;
        int iL = idx * Nu + d0;
        
        if ((Nu2 < 0) || (Nu2 > 2 * Nu) || (iL < 0) || (iL + Nu2 > Nb2)) continue;

        double *s_per_codeword = new double[Q];
        // 全ての符号語xiについて、観測確率 P(y|xi,d0,d1) を計算
        for (int xi = 0; xi < Q; xi++) {
          if (Nu2 >= Nu2min && Nu2 <= Nu2max) {
            long y = VectToLong(&Yin[iL], Nu2);
            s_per_codeword[xi] = GetGX(Nu2, y, xi);
          } else {
            s_per_codeword[xi] = GetGX(Nu2, 0, 0); // 近似値
          }
        }

        // 全ての次のエラー状態e1への遷移を計算
        for (int e1 = 0; e1 < NUM_ERROR_STATES; e1++) {
          double total_gamma = 0.0;
          // 全ての符号語xiについて、分岐メトリックを計算し、周辺化(合計)する
          for (int xi = 0; xi < Q; xi++) {
            // 論文の3Dラティス計算: P(y|xi) * P(xi) * P(e1|e0)
            double gamma_xi = s_per_codeword[xi] * PU[idx][xi] * PE[idx][e0][e1];
            total_gamma += gamma_xi;
          }
          
          // 前進確率の更新
          PFE[idx+1][d1-Dmin][e1] += PFE[idx][d0-Dmin][e0] * total_gamma;
        }
        
        delete[] s_per_codeword;
      }
    }
  }

  // 正規化
  double total_sum = 0.0;
  for(int d1 = Dmin; d1 <= Dmax; d1++){
    for(int e1 = 0; e1 < NUM_ERROR_STATES; e1++){
      total_sum += PFE[idx+1][d1-Dmin][e1];
    }
  }
  if(total_sum > 0.0){
    for(int d1 = Dmin; d1 <= Dmax; d1++){
      for(int e1 = 0; e1 < NUM_ERROR_STATES; e1++){
        PFE[idx+1][d1-Dmin][e1] /= total_sum;
      }
    }
  }
}

//================================================================================
// エラー状態を含む後進確率計算（3Dラティス）
//================================================================================
void SLFBAdec::CalcPBE(int idx, int Nb2){
  assert(idx >= 0 && idx < Ns);

  // エラー状態遷移確率を事前計算
  CalcPE(idx);
  
  // 現在の状態 t の確率を初期化
  for (int d0 = Dmin; d0 <= Dmax; d0++) {
    for (int e0 = 0; e0 < NUM_ERROR_STATES; e0++) {
      PBE[idx][d0-Dmin][e0] = 0.0;
    }
  }
  
  // 次の状態 t+1 から 現在の状態 t への全ての遷移を計算
  for (int d1 = Dmin; d1 <= Dmax; d1++) {
    for (int e1 = 0; e1 < NUM_ERROR_STATES; e1++) {
      
      if (PBE[idx+1][d1-Dmin][e1] == 0.0) continue;

      for (int d0 = Dmin; d0 <= Dmax; d0++) {
          
        // --- 分岐メトリック γ_t の計算 (CalcPFEと同一) ---
        int Nu2 = Nu + d1 - d0;
        int iL = idx * Nu + d0;
        
        if ((Nu2 < 0) || (Nu2 > 2 * Nu) || (iL < 0) || (iL + Nu2 > Nb2)) continue;

        double *s_per_codeword = new double[Q];
        for (int xi = 0; xi < Q; xi++) {
          if (Nu2 >= Nu2min && Nu2 <= Nu2max) {
            long y = VectToLong(&Yin[iL], Nu2);
            s_per_codeword[xi] = GetGX(Nu2, y, xi);
          } else {
            s_per_codeword[xi] = GetGX(Nu2, 0, 0);
          }
        }
        
        // 全ての前のエラー状態e0からの遷移を計算
        for (int e0 = 0; e0 < NUM_ERROR_STATES; e0++) {
          double total_gamma = 0.0;
          for (int xi = 0; xi < Q; xi++) {
            double gamma_xi = s_per_codeword[xi] * PU[idx][xi] * PE[idx][e0][e1];
            total_gamma += gamma_xi;
          }
          
          // 後進確率の更新
          PBE[idx][d0-Dmin][e0] += PBE[idx+1][d1-Dmin][e1] * total_gamma;
        }
        
        delete[] s_per_codeword;
      }
    }
  }

  // 正規化
  double total_sum = 0.0;
  for(int d0 = Dmin; d0 <= Dmax; d0++){
    for(int e0 = 0; e0 < NUM_ERROR_STATES; e0++){
      total_sum += PBE[idx][d0-Dmin][e0];
    }
  }
  if(total_sum > 0.0){
    for(int d0 = Dmin; d0 <= Dmax; d0++){
      for(int e0 = 0; e0 < NUM_ERROR_STATES; e0++){
        PBE[idx][d0-Dmin][e0] /= total_sum;
      }
    }
  }
}

//================================================================================
// エラー状態を含む事後確率計算（3Dラティス）
//================================================================================
void SLFBAdec::CalcPDE(int idx, int Nb2){
  assert(idx >= 0 && idx < Ns);
  
  // エラー状態遷移確率を事前計算
  CalcPE(idx);
  
  // 各符号語xiの事後確率を計算
  for (int xi = 0; xi < Q; xi++) {
    PD[idx][xi] = 0.0;
    
    // 全てのドリフト状態とエラー状態の組み合わせについて周辺化
    for (int d0 = Dmin; d0 <= Dmax; d0++) {
      for (int e0 = 0; e0 < NUM_ERROR_STATES; e0++) {
        for (int d1 = Dmin; d1 <= Dmax; d1++) {
          for (int e1 = 0; e1 < NUM_ERROR_STATES; e1++) {
            
            // 基本チャネルパラメータ計算
            int Nu2 = Nu + d1 - d0;
            int iL = idx * Nu + d0;
            
            if ((Nu2 < 0) || (Nu2 > 2 * Nu) || (iL < 0) || (iL + Nu2 > Nb2)) continue;
            
            // 観測確率 P(y|xi, d0, d1)
            double obs_prob;
            if (Nu2 >= Nu2min && Nu2 <= Nu2max) {
              long y = VectToLong(&Yin[iL], Nu2);
              obs_prob = GetGX(Nu2, y, xi);
            } else {
              obs_prob = GetGX(Nu2, 0, 0);
            }
            
            // 3D前進・後進確率を使った事後確率計算
            // P(xi|y,d,e) ∝ α(d0,e0) × [ P(y|xi,d) × P(xi) × P(e1|e0) ] × β(d1,e1)
            double posterior_contrib = PFE[idx][d0-Dmin][e0] * 
                                     obs_prob * 
                                     PU[idx][xi] * 
                                     PE[idx][e0][e1] * 
                                     PBE[idx+1][d1-Dmin][e1];
            
            PD[idx][xi] += posterior_contrib;
          }
        }
      }
    }
  }
  
  // 正規化（符号語確率の合計を1にする）
  normalize(PD[idx], Q);
}

