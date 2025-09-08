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
// 4値シンボルの置換確率を計算 (論文Table 2対応)
//================================================================================
double SLFBAdec::Psub_quaternary(unsigned char a, unsigned char b) {
  assert(a <= 3 && b <= 3);
  
  if (a == b) {
    // 置換が起こらない確率 (Match)
    return 1.0 - Ps; 
  } else {
    // aがbに置換される確率
    // Ps(置換イベントが起こる確率) * P(その結果bになる条件付き確率)
    return Ps * SubMatrix[a][b];
  }
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
    assert(V[i]>=0 && V[i]<=3);  // 4元シンボル(0-3)に対応
    val <<= 2;  // 4進数なので2ビットシフト
    val |= V[i];  // シンボル値をそのまま使用
  } // for i
  return val;
}

//================================================================================
void SLFBAdec::LongToVect(unsigned char *V, long val, int len){
  assert(len>0 && val>=0);
  long mask = 0x3 << ((len-1)*2);  // 4進数なので2ビット×(len-1)
  for(int i=0;i<len;i++){
    V[i] = (val & mask) >> ((len-1-i)*2);  // 2ビットシフトで4元シンボル抽出
    mask >>= 2;  // 次の2ビット位置へ
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

  printf("# SetGX: Starting pre-calculation of the GX table. This may take a very long time...\n");
  
  for(int ly=0;ly<=Nu*2;ly++){
    if(ly<Nu2min || ly>Nu2max){
      // approximate
      GX[ly]    = new double * [1];
      GX[ly][0] = new double [1];
      GX[ly][0][0] = (ly<Nu)? pow(Pd,Nu-ly) : pow(Pi,ly-Nu); //?
      //printf("GX[%d]=%e\n",ly,GX[ly][0][0]);
    } else {
      // exact
      // ly2p = (long)pow(2,ly);
      ly2p = 1;          // ✅ 修正：4^lyを計算
      for(int i=0; i<ly; i++) ly2p *= 4;

      // ▼▼▼ 進捗表示：ステップA (どの `ly` を処理中か表示) ▼▼▼
      printf("# SetGX: Processing received length ly = %d (y has %ld patterns)...\n", ly, ly2p);
      fflush(stdout); // バッファを強制的に出力して、メッセージをすぐに表示させる

      GX[ly] = new double * [ly2p];
      for(long y=0;y<ly2p;y++){

        // ▼▼▼ 進捗表示：ステップB (パーセンテージ表示) ▼▼▼
        // 計算回数が多い場合のみ、約1%ごとに進捗を表示する
        if (ly2p > 1000 && y % (ly2p / 100) == 0) {
            double percent = (double)y / ly2p * 100.0;
            printf("\r  -> Progress: [%.1f %%]", percent); // \r でカーソルを行頭に戻す
            fflush(stdout); // バッファを強制的に出力
        }

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
      // ▼▼▼ 進捗表示：ステップC (完了表示) ▼▼▼
      printf("\r  -> Progress: [100.0 %%] ... Done.                     \n");
    } // if ly
  } // for ly
  
  printf("# SetGX: Pre-calculation of the GX table is complete.\n");
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
      ret = Psub_quaternary(X[0],Y[0]) * Pt;
    } else if(ly==2) {
      ret = Psub_quaternary(X[0],Y[0]) * Psub_quaternary(X[0],Y[1]) * Pi;
    } else {
      assert(ly>=3); 
      ret = 0.0;
    } // if ly
  } else {
    x1  = VectToLong(&X[1],lx-1);
    // -----trans
    qt = Psub_quaternary(X[0],Y[0]) * Pt;
    y1 = (ly==1)? 0 : VectToLong(&Y[1],ly-1);
    qt *= CalcPyx(y1,x1,ly-1,lx-1);
    // -----ins
    if(ly<2){
      qi = 0.0;
    } else {
      qi = Psub_quaternary(X[0],Y[0]) * Psub_quaternary(X[0],Y[1]) * Pi;
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
  // Pi, Pd, Ps は Legacy/Dec2 用のメンバ変数として保持
  this->Pi = CH->GetPi();
  this->Pd = CH->GetPd();
  this->Ps = CH->GetPs();
  this->Pt = 1.0 - this->Pi - this->Pd;

  // Nu2p = (long) pow(2,Nu);
  Ns   = Nb/Nu;
  Drng = Dmax-Dmin+1;
  Nu2max = Nu + (int)ceil( (double)Nu*Pi ) + 2;
  Nu2min = Nu - (int)ceil( (double)Nu*Pd ) - 2;
  Nu2max = min(Nu2max, Nu*2);
  Nu2min = max(Nu2min, 0   );
  // Dec3: k-mer関連の初期化
  num_kmers = 1;
  for(int i = 0; i < KMER_LENGTH; i++) {
    num_kmers *= 4;                  // 4^KMER_LENGTH = 256 for k=4
  }
  
  // Dec3: 4次元ラティス用ポインタの初期化
  PFE4D = nullptr;
  PBE4D = nullptr;
  kmer_error_probs = nullptr;
  
  // 論文 Table 2 に基づく4値置換確率行列の初期化
  // P(b|a) = aがbに置換される条件付き確率 (A=0, C=1, G=2, T=3)
  //           A      C      G      T
  SubMatrix[0][0]=0.0;   SubMatrix[0][1]=0.149; SubMatrix[0][2]=0.675; SubMatrix[0][3]=0.176; // Aから
  SubMatrix[1][0]=0.351; SubMatrix[1][1]=0.0;   SubMatrix[1][2]=0.173; SubMatrix[1][3]=0.476; // Cから
  SubMatrix[2][0]=0.756; SubMatrix[2][1]=0.076; SubMatrix[2][2]=0.0;   SubMatrix[2][3]=0.168; // Gから
  SubMatrix[3][0]=0.328; SubMatrix[3][1]=0.424; SubMatrix[3][2]=0.248; SubMatrix[3][3]=0.0;   // Tから
  
  assert(Nb%Nu==0);
  assert(Nu<=NuMax);
  assert(ECM->GetM()==Q && ECM->GetN()==Q);
  assert(Pt>0.0 && Pt<=1.0);
  //----- set tables
  const char* decoder_mode = getenv("DECODER_MODE");
  if (decoder_mode == nullptr) decoder_mode = "DEC3"; // デフォルトをDEC3に
  
   if (strcmp(decoder_mode, "DEC3") == 0) {
      // --- Dec3用の初期化 ---
      printf("# SLFBAdec: Initializing for DEC3 mode.\n");
      printf("# SLFBAdec: Skipping SetGD() and SetGX() pre-calculation.\n");
      SetFG();
      SetFGE();
      SetFGE4D(); 
      LoadKmerErrorProbabilities("DNArSim-main"); 
  } else {
      // --- Legacy/Dec2用の初期化 ---
      printf("# SLFBAdec: Initializing for LEGACY/DEC2 mode.\n");
      printf("# SLFBAdec: Performing SetGD() and SetGX() pre-calculation. This may take a very long time...\n");
      Nu2max = Nu + (int)ceil((double)Nu*Pi) + 2;
      Nu2min = Nu - (int)ceil((double)Nu*Pd) - 2;
      Nu2max = min(Nu2max, Nu*2);
      Nu2min = max(Nu2min, 0);

      SetGD();
      SetGX(); // 巨大なメモリ確保と計算爆発が起こる可能性がある
      SetFG();
      SetFGE();
  }
}

//================================================================================
SLFBAdec::~SLFBAdec(){
  DelGD();
  DelGX();
  DelFG();
  DelFGE();
  DelFGE4D();  // Dec3: 4次元ラティスメモリ解放
  printf("# SLFBAdec: deleted (including Dec3 4D lattice)\n");
}

//================================================================================
void SLFBAdec::Decode(double **Pout, const unsigned char *RW, int Nb2, const int *dbgIW, const double **Pin){
  int idx;
  assert( Nb2>=Nb+Dmin && Nb2<=Nb+Dmax );
  
  // デコーダー選択: 環境変数DECODER_MODEで制御
  // DECODER_MODE=LEGACY: 従来の2Dデコーダー（ベースライン）
  // DECODER_MODE=DEC2:   新しい3Dラティスデコーダー（Dec2）
  // DECODER_MODE=DEC3:   究極の4Dラティスデコーダー（Dec3、k-mer依存）
  const char* decoder_mode = getenv("DECODER_MODE");
  if (decoder_mode == nullptr) {
    decoder_mode = "DEC2"; // デフォルトをDec2に設定
  }
  
  if (strcmp(decoder_mode, "LEGACY") == 0) {
    printf("# Using LEGACY 2D decoder (baseline)\n");
    InitFG(RW,Pin,Nb2);
    for(idx=0;   idx<Ns;idx++) CalcPU(idx);
    for(idx=0;   idx<Ns;idx++) CalcPF(idx,Nb2);
    for(idx=Ns-1;idx>=0;idx--) CalcPB(idx,Nb2);
    for(idx=0;   idx<Ns;idx++) CalcPD(idx,Nb2);
    for(idx=0;   idx<Ns;idx++) CalcPO(idx);
  } else if (strcmp(decoder_mode, "DEC2") == 0) {
    printf("# Using Dec2 3D lattice decoder with error state memory\n");
    InitFGE(RW,Pin,Nb2);
    for(idx=0;   idx<Ns;idx++) CalcPU(idx);
    for(idx=0;   idx<Ns;idx++) CalcPFE(idx,Nb2);   // 3D前進確率
    for(idx=Ns-1;idx>=0;idx--) CalcPBE(idx,Nb2);   // 3D後進確率
    for(idx=0;   idx<Ns;idx++) CalcPDE(idx,Nb2);   // 3D事後確率
    for(idx=0;   idx<Ns;idx++) CalcPO(idx);
  } else if (strcmp(decoder_mode, "DEC3") == 0) {
    printf("# Using Dec3 4D lattice decoder with k-mer dependency (ultimate decoder)\n");
    // ✅ メモ化：新しいデコードの開始前にキャッシュをクリア
    pyx_cache.clear();
    InitFGE4D(RW,Pin,Nb2);  // 4D lattice initialization
    for(idx=0;   idx<Ns;idx++) CalcPU(idx);
    for(idx=0;   idx<Ns;idx++) CalcPFE4D(idx,Nb2);   // 4D前進確率
    for(idx=Ns-1;idx>=0;idx--) CalcPBE4D(idx,Nb2);   // 4D後進確率
    for(idx=0;   idx<Ns;idx++) CalcPDE4D(idx,Nb2);   // 4D事後確率
    for(idx=0;   idx<Ns;idx++) CalcPO(idx);
    DelFGE4D();  // 4D lattice memory deallocation
  } else {
    printf("# Error: Unknown DECODER_MODE '%s', defaulting to Dec2\n", decoder_mode);
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

//================================================================================
//================================================================================
//=== Dec3: k-mer依存4次元ラティスデコーダ実装 ===================================
//================================================================================
//================================================================================

#include <map>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <iostream>

//================================================================================
int SLFBAdec::GetKmerIndex(const unsigned char *kmer_seq) {
  int idx = 0;
  for(int i = 0; i < KMER_LENGTH; i++) {
    assert(kmer_seq[i] <= 3);  // 0=A, 1=C, 2=G, 3=T
    idx = idx * 4 + kmer_seq[i];
  }
  return idx;
}

//================================================================================
void SLFBAdec::GetKmerFromIndex(int kmer_idx, unsigned char *kmer_seq) {
  for(int i = KMER_LENGTH - 1; i >= 0; i--) {
    kmer_seq[i] = kmer_idx % 4;
    kmer_idx /= 4;
  }
}

//================================================================================
// 決定論的k-mer遷移の計算
// 入力: current_kmer (k_t), codeword_xi (符号語インデックス)
// 出力: next_kmer (k_{t+1}) - 一意に決まる
//================================================================================
int SLFBAdec::ComputeNextKmer(int current_kmer, int codeword_xi) {
  assert(current_kmer >= 0 && current_kmer < num_kmers);
  assert(codeword_xi >= 0 && codeword_xi < Q);
  
  // 現在のk-merをbit配列に変換
  unsigned char current_kmer_seq[KMER_LENGTH];
  GetKmerFromIndex(current_kmer, current_kmer_seq);
  
  // 符号語xiをbit配列に取得
  unsigned char codeword[Nu];
  ICB->Get_CW(codeword, codeword_xi);
  
  // 次のk-mer計算: 左シフト + 新しい塩基追加
  unsigned char next_kmer_seq[KMER_LENGTH];
  
  // KMER_LENGTH-1個の古い塩基を左シフト
  for(int i = 0; i < KMER_LENGTH - 1; i++) {
    next_kmer_seq[i] = current_kmer_seq[i + 1];
  }
  
  // ✅ 修正：新しい塩基を末尾に追加
  // 符号語はすでに4値なので、先頭の1シンボルをそのまま使用する
  next_kmer_seq[KMER_LENGTH - 1] = codeword[0];
  
  // 次のk-merインデックスに変換
  return GetKmerIndex(next_kmer_seq);
}

//================================================================================
void SLFBAdec::LoadKmerErrorProbabilities(const char* dir_path) {
  // std::mapを動的に作成
  std::map<std::pair<int,int>, std::vector<double>>* prob_map = 
    new std::map<std::pair<int,int>, std::vector<double>>();
  kmer_error_probs = (void*)prob_map;
  
  printf("# Dec3: Loading REAL P(e_{t+1}|e_t, η_{t+1}) from DNArSim-main\n");
  
  // DNArSim-mainの確率ファイルパスを構築  
  char prob_dir[512];
  snprintf(prob_dir, sizeof(prob_dir), "%s/simulator/probEdit/k%d", dir_path, KMER_LENGTH);
  printf("# Dec3: Reading from: %s/\n", prob_dir);
  
  // 4つのエラー状態に対応するファイル名
  const char* error_files[NUM_ERROR_STATES] = {
    "KmerYi_prevYiM_RatesAvg.txt",  // ERROR_MATCH = 0
    "KmerYi_prevYiI_RatesAvg.txt",  // ERROR_INSERTION = 1  
    "KmerYi_prevYiD_RatesAvg.txt",  // ERROR_DELETION = 2
    "KmerYi_prevYiS_RatesAvg.txt"   // ERROR_SUBSTITUTION = 3
  };
  
  int loaded_entries = 0;
  int file_errors = 0;
  
  // 各エラー状態に対してファイルを読み込み
  for(int prev_error = 0; prev_error < NUM_ERROR_STATES; prev_error++) {
    char filepath[1024];
    snprintf(filepath, sizeof(filepath), "%s/%s", prob_dir, error_files[prev_error]);
    
    FILE* file = fopen(filepath, "r");
    if(!file) {
      printf("# Dec3: Warning - Could not open %s\n", filepath);
      file_errors++;
      continue;
    }
    
    printf("# Dec3: Reading %s...\n", error_files[prev_error]);
    
    char line[256];
    char kmer_str[16];
    double ins_prob, del_prob, subst_prob, err_prob, match_prob;
    
    while(fgets(line, sizeof(line), file)) {
      // ファイル形式: "KMER Ins Del Subst Err Match"
      if(sscanf(line, "%s %lf %lf %lf %lf %lf", 
                kmer_str, &ins_prob, &del_prob, &subst_prob, &err_prob, &match_prob) == 6) {
        
        // k-mer文字列をインデックスに変換
        if(strlen(kmer_str) != KMER_LENGTH) continue;
        
        unsigned char kmer_seq[KMER_LENGTH];
        bool valid_kmer = true;
        for(int i = 0; i < KMER_LENGTH; i++) {
          switch(kmer_str[i]) {
            case 'A': kmer_seq[i] = 0; break;
            case 'C': kmer_seq[i] = 1; break;
            case 'G': kmer_seq[i] = 2; break;
            case 'T': kmer_seq[i] = 3; break;
            default: valid_kmer = false; break;
          }
        }
        
        if(!valid_kmer) continue;
        
        int kmer_idx = GetKmerIndex(kmer_seq);
        
        // 確率ベクトルを構築
        std::pair<int,int> key = std::make_pair(prev_error, kmer_idx);
        std::vector<double> probs(NUM_ERROR_STATES);
        
        probs[ERROR_MATCH] = match_prob;           // P(e_{t+1}=Match | e_t, η_{t+1})
        probs[ERROR_INSERTION] = ins_prob;         // P(e_{t+1}=Insertion | e_t, η_{t+1})
        probs[ERROR_DELETION] = del_prob;          // P(e_{t+1}=Deletion | e_t, η_{t+1})
        probs[ERROR_SUBSTITUTION] = subst_prob;    // P(e_{t+1}=Substitution | e_t, η_{t+1})
        
        (*prob_map)[key] = probs;
        loaded_entries++;
      }
    }
    
    fclose(file);
  }
  
  // 読み込めなかったエントリに対してフォールバック値を設定
  for(int prev_error = 0; prev_error < NUM_ERROR_STATES; prev_error++) {
    for(int kmer_idx = 0; kmer_idx < num_kmers; kmer_idx++) {
      std::pair<int,int> key = std::make_pair(prev_error, kmer_idx);
      
      if(prob_map->find(key) == prob_map->end()) {
        // フォールバック: 理論的確率
        std::vector<double> probs(NUM_ERROR_STATES);
        probs[ERROR_MATCH] = 0.85;
        probs[ERROR_INSERTION] = 0.05;
        probs[ERROR_DELETION] = 0.05; 
        probs[ERROR_SUBSTITUTION] = 0.05;
        (*prob_map)[key] = probs;
        file_errors++;
      }
    }
  }
  
  printf("# Dec3: Successfully loaded %d REAL probability entries from DNArSim-main\n", loaded_entries);
  printf("# Dec3: Used fallback for %d entries\n", file_errors);
  printf("# Dec3: TRUE P(e_{t+1}|e_t, η_{t+1}) table ready for ultimate decoding!\n");
}

//================================================================================
double SLFBAdec::GetKmerErrorProb(int prev_error, int next_kmer, int next_error) {
  assert(prev_error >= 0 && prev_error < NUM_ERROR_STATES);
  assert(next_kmer >= 0 && next_kmer < num_kmers);
  assert(next_error >= 0 && next_error < NUM_ERROR_STATES);
  
  std::map<std::pair<int,int>, std::vector<double>>* prob_map = 
    (std::map<std::pair<int,int>, std::vector<double>>*)kmer_error_probs;
  
  std::pair<int,int> key = std::make_pair(prev_error, next_kmer);
  auto it = prob_map->find(key);
  
  if(it != prob_map->end()) {
    return it->second[next_error];
  } else {
    // フォールバック: 対応するエントリがない場合のデフォルト確率
    if(next_error == ERROR_MATCH) return 0.70;
    else return 0.10;
  }
}

//================================================================================
// Dec3究極統合: 動的確率を使った観測確率計算 (CalcPyxの進化版)
// GetGXの完全な代替として、DNArSim現実確率に基づくCalcPyx動的計算
//================================================================================
double SLFBAdec::CalcPyx_dynamic(long y, long x, int ly, int lx, int prev_error, int kmer, int codeword_xi) {
  assert(lx > 0 && ly >= 0);
  assert(prev_error >= 0 && prev_error < NUM_ERROR_STATES);
  assert(kmer >= 0 && kmer < num_kmers);
  assert(codeword_xi >= 0 && codeword_xi < Q);
  // ▼▼▼ メモ化：ステップA（キャッシュ確認） ▼▼▼
  // 1. 現在の引数の組み合わせをキーとして作成
  std::tuple<long, long, int, int, int, int, int> key = 
      std::make_tuple(y, x, ly, lx, prev_error, kmer, codeword_xi);

  // 2. キャッシュにキーが存在するか確認
  if (pyx_cache.count(key)) {
      return pyx_cache[key]; // 存在すれば、計算せずにすぐに結果を返す
  }
  // ▲▲▲ メモ化：ステップAここまで ▲▲▲

  // ★ 状態に応じた動的な確率をDNArSim-mainから取得
  double pi = GetKmerErrorProb(prev_error, kmer, ERROR_INSERTION);
  double pd = GetKmerErrorProb(prev_error, kmer, ERROR_DELETION);  
  double ps = GetKmerErrorProb(prev_error, kmer, ERROR_SUBSTITUTION);
  double pt = 1.0 - pi - pd; // 伝送確率（Match + 実際のSubstitutionを含む）

  // 削除のみの場合
  if (ly == 0) return pow(pd, lx);

  long x1, y1;
  double ret, qt, qi, qd;
  unsigned char *X = new unsigned char[lx];
  unsigned char *Y = new unsigned char[ly];
  LongToVect(X, x, lx);
  LongToVect(Y, y, ly);

  if (lx == 1) {
    if (ly == 1) {
      // 1対1の場合：動的psを使ったPsub計算
      ret = Psub_quaternary(X[0], Y[0], ps) * pt;
    } else if (ly == 2) {
      // 1対2の場合：挿入を含む
      ret = Psub_quaternary(X[0], Y[0], ps) * Psub_quaternary(X[0], Y[1], ps) * pi;
    } else {
      // 1対3以上は確率0
      ret = 0.0;
    }
  } else {
    // 再帰的計算
    x1 = VectToLong(&X[1], lx - 1);
    
    // ✅ 修正：次のk-merは一度だけ計算
    int next_kmer = ComputeNextKmer(kmer, codeword_xi);
    // -----伝送 (transmission)
    qt = Psub_quaternary(X[0], Y[0], ps) * pt;
    y1 = (ly == 1) ? 0 : VectToLong(&Y[1], ly - 1);
    qt *= CalcPyx_dynamic(y1, x1, ly - 1, lx - 1, ERROR_MATCH, next_kmer, codeword_xi);

    // -----挿入 (insertion)
    if (ly < 2) {
      qi = 0.0;
    } else {
      qi = Psub_quaternary(X[0], Y[0], ps) * Psub_quaternary(X[0], Y[1], ps) * pi;
      y1 = (ly == 2) ? 0 : VectToLong(&Y[2], ly - 2);
      qi *= CalcPyx_dynamic(y1, x1, ly - 2, lx - 1, ERROR_INSERTION, next_kmer, codeword_xi);
    }

    // -----削除 (deletion)
    qd = pd;
    y1 = y;
    qd *= CalcPyx_dynamic(y1, x1, ly, lx - 1, ERROR_DELETION, next_kmer, codeword_xi);

    ret = qt + qi + qd;
  }

  delete[] X;
  delete[] Y;
  
  // ▼▼▼ メモ化：ステップB（結果を保存） ▼▼▼
  // 3. 計算結果をキャッシュに保存してから返す
  pyx_cache[key] = ret;
  return ret;
  // ▲▲▲ メモ化：ステップBここまで ▲▲▲
}

//================================================================================
// 動的ps対応のPsub_quaternaryオーバーロード
//================================================================================
double SLFBAdec::Psub_quaternary(unsigned char a, unsigned char b, double dynamic_ps) {
  assert(a <= 3 && b <= 3);
  
  if (a == b) {
    // 置換が起こらない確率 (Match)
    return 1.0 - dynamic_ps; 
  } else {
    // aがbに置換される確率
    return dynamic_ps * SubMatrix[a][b];
  }
}

//================================================================================
// Dec3統合: DNAチャネル統合観測確率計算（CalcPyx_dynamicのラッパー）
//================================================================================
double SLFBAdec::ComputeObservationProbabilityFromDNA(int idx, int Nu2, int iL, int xi, int k0, int e0, int Nb2) {
  // 範囲チェック
  if ((Nu2 < 0) || (Nu2 > 2 * Nu) || (iL < 0) || (iL + Nu2 > Nb2)) {
    return 0.001;
  }
  
  // ★ 真の統合：CalcPyx_dynamicを呼び出してDNArSim現実確率で計算
  if (Nu2 >= Nu2min && Nu2 <= Nu2max) {
    long y = VectToLong(&Yin[iL], Nu2);
    
    // 符号語を取得してlongに変換
    unsigned char *codeword = new unsigned char[Nu];
    ICB->Get_CW(codeword, xi);
    long x = VectToLong(codeword, Nu);
    
    // k0, e0情報がある場合はそれを使用、ない場合はデフォルト値
    int actual_prev_error = (e0 >= 0) ? e0 : ERROR_MATCH;
    int actual_kmer = (k0 >= 0) ? k0 : 0;
    
    double result = CalcPyx_dynamic(y, x, Nu2, Nu, actual_prev_error, actual_kmer, xi);
    delete[] codeword;
    return result;
  } else {
    // Nu2が範囲外の場合は均等確率
    return 1.0 / Q;
  }
}

//================================================================================
void SLFBAdec::SetFGE4D() {
  printf("# Dec3: Setting up 4D lattice memory structures\n");
  printf("# Dec3: Dimensions: [Ns+1=%d][Drng=%d][NUM_ERROR_STATES=%d][num_kmers=%d]\n", 
         Ns+1, Drng, NUM_ERROR_STATES, num_kmers);
  
  // PFE4D: [Ns+1][Drng][NUM_ERROR_STATES][num_kmers]
  PFE4D = new double***[Ns+1];
  for(int i = 0; i < Ns+1; i++) {
    PFE4D[i] = new double**[Drng];
    for(int d = 0; d < Drng; d++) {
      PFE4D[i][d] = new double*[NUM_ERROR_STATES];
      for(int e = 0; e < NUM_ERROR_STATES; e++) {
        PFE4D[i][d][e] = new double[num_kmers];
        // 初期化
        for(int k = 0; k < num_kmers; k++) {
          PFE4D[i][d][e][k] = 0.0;
        }
      }
    }
  }
  
  // PBE4D: [Ns+1][Drng][NUM_ERROR_STATES][num_kmers]
  PBE4D = new double***[Ns+1];
  for(int i = 0; i < Ns+1; i++) {
    PBE4D[i] = new double**[Drng];
    for(int d = 0; d < Drng; d++) {
      PBE4D[i][d] = new double*[NUM_ERROR_STATES];
      for(int e = 0; e < NUM_ERROR_STATES; e++) {
        PBE4D[i][d][e] = new double[num_kmers];
        // 初期化
        for(int k = 0; k < num_kmers; k++) {
          PBE4D[i][d][e][k] = 0.0;
        }
      }
    }
  }
  
  size_t memory_size = (size_t)(Ns+1) * Drng * NUM_ERROR_STATES * num_kmers * sizeof(double) * 2;
  printf("# Dec3: 4D lattice memory allocated: %.2f MB\n", memory_size / (1024.0 * 1024.0));
}

//================================================================================
void SLFBAdec::DelFGE4D() {
  printf("# Dec3: Deleting 4D lattice memory structures\n");
  
  // PFE4D削除
  if(PFE4D != nullptr) {
    for(int i = 0; i < Ns+1; i++) {
      if(PFE4D[i] != nullptr) {
        for(int d = 0; d < Drng; d++) {
          if(PFE4D[i][d] != nullptr) {
            for(int e = 0; e < NUM_ERROR_STATES; e++) {
              delete[] PFE4D[i][d][e];
            }
            delete[] PFE4D[i][d];
          }
        }
        delete[] PFE4D[i];
      }
    }
    delete[] PFE4D;
    PFE4D = nullptr;
  }
  
  // PBE4D削除
  if(PBE4D != nullptr) {
    for(int i = 0; i < Ns+1; i++) {
      if(PBE4D[i] != nullptr) {
        for(int d = 0; d < Drng; d++) {
          if(PBE4D[i][d] != nullptr) {
            for(int e = 0; e < NUM_ERROR_STATES; e++) {
              delete[] PBE4D[i][d][e];
            }
            delete[] PBE4D[i][d];
          }
        }
        delete[] PBE4D[i];
      }
    }
    delete[] PBE4D;
    PBE4D = nullptr;
  }
  
  // k-mer確率テーブル削除
  if(kmer_error_probs != nullptr) {
    std::map<std::pair<int,int>, std::vector<double>>* prob_map = 
      (std::map<std::pair<int,int>, std::vector<double>>*)kmer_error_probs;
    delete prob_map;
    kmer_error_probs = nullptr;
  }
}

//================================================================================
void SLFBAdec::InitFGE4D(const unsigned char *RW, const double **Pin, int Nb2) {
  // 既存の初期化を実行
  InitFGE(RW, Pin, Nb2);
  
  // 4次元確率配列を初期化
  for(int i = 0; i < Ns+1; i++) {
    for(int d = 0; d < Drng; d++) {
      for(int e = 0; e < NUM_ERROR_STATES; e++) {
        for(int k = 0; k < num_kmers; k++) {
          PFE4D[i][d][e][k] = 0.0;
          PBE4D[i][d][e][k] = 0.0;
        }
      }
    }
  }
  
  // 4次元初期状態: 論文に従い(d=0, e=Match, η=0(AAAA))で開始
  PFE4D[0][0-Dmin][ERROR_MATCH][0] = 1.0;  // σ_0 = (0, M, AAAA)
  
  // 4次元最終状態: 全k-mer状態に均等分布
  for(int k = 0; k < num_kmers; k++) {
    for(int e = 0; e < NUM_ERROR_STATES; e++) {
      PBE4D[Ns][Nb2-Nb-Dmin][e][k] = 1.0 / (num_kmers * NUM_ERROR_STATES);
    }
  }
  
  printf("# Dec3: InitFGE4D completed - 4D lattice initialized\n");
}

//================================================================================
// Dec3核心：4次元前進確率計算 CalcPFE4D
// 状態: (d_t, e_t, η_t) → (d_{t+1}, e_{t+1}, η_{t+1})
// 遷移確率: P(e_{t+1}|e_t, η_{t+1}) を使用
//================================================================================
void SLFBAdec::CalcPFE4D(int idx, int Nb2) {
  assert(idx >= 0 && idx < Ns);

  // 次の状態 t+1 の確率を初期化
  for(int d1 = Dmin; d1 <= Dmax; d1++) {
    for(int e1 = 0; e1 < NUM_ERROR_STATES; e1++) {
      for(int k1 = 0; k1 < num_kmers; k1++) {
        PFE4D[idx+1][d1-Dmin][e1][k1] = 0.0;
      }
    }
  }
  
  // 4次元状態遷移：全ての (d0,e0,k0) → (d1,e1,k1) を計算
  for (int d0 = Dmin; d0 <= Dmax; d0++) {
    for (int e0 = 0; e0 < NUM_ERROR_STATES; e0++) {
      for (int k0 = 0; k0 < num_kmers; k0++) {
        
        // 現在の状態への確率が0なら、このパスからの寄与はない
        if (PFE4D[idx][d0-Dmin][e0][k0] == 0.0) continue;

        for (int d1 = Dmin; d1 <= Dmax; d1++) {
          
          // チャネルパラメータ計算
          int Nu2 = Nu + d1 - d0;
          int iL = idx * Nu + d0;
          
          if ((Nu2 < 0) || (Nu2 > 2 * Nu) || (iL < 0) || (iL + Nu2 > Nb2)) continue;

          // ✅ Dec3新方式: DNAチャネル統合動的確率計算
          // GetGX依存を排除し、実際のDNA配列とk-mer確率のみを使用
          double *s_per_codeword = new double[Q];
          for (int xi = 0; xi < Q; xi++) {
            // 実DNA配列から直接観測確率を計算
            s_per_codeword[xi] = ComputeObservationProbabilityFromDNA(idx, Nu2, iL, xi, k0, e0, Nb2);
          }

          // ✅ Dec3効率化修正：全符号語xiで決定論的k-mer遷移を直接計算
          for (int xi = 0; xi < Q; xi++) {
            int k1 = ComputeNextKmer(k0, xi); // 決定論的遷移: k0 + xi → k1

            for (int e1 = 0; e1 < NUM_ERROR_STATES; e1++) {
              double kmer_error_prob = GetKmerErrorProb(e0, k1, e1);
              double gamma_xi = s_per_codeword[xi] * PU[idx][xi] * kmer_error_prob;
              
              // 4D前進確率の更新（符号語xiによって決まるk1状態に蓄積）
              PFE4D[idx+1][d1-Dmin][e1][k1] += PFE4D[idx][d0-Dmin][e0][k0] * gamma_xi;
            }
          }
          
          delete[] s_per_codeword;
        }
      }
    }
  }

  // 4次元正規化
  double total_sum = 0.0;
  for(int d1 = Dmin; d1 <= Dmax; d1++) {
    for(int e1 = 0; e1 < NUM_ERROR_STATES; e1++) {
      for(int k1 = 0; k1 < num_kmers; k1++) {
        total_sum += PFE4D[idx+1][d1-Dmin][e1][k1];
      }
    }
  }
  if(total_sum > 0.0) {
    for(int d1 = Dmin; d1 <= Dmax; d1++) {
      for(int e1 = 0; e1 < NUM_ERROR_STATES; e1++) {
        for(int k1 = 0; k1 < num_kmers; k1++) {
          PFE4D[idx+1][d1-Dmin][e1][k1] /= total_sum;
        }
      }
    }
  }
  
  printf("# Dec3: CalcPFE4D[%d] completed with P(e_{t+1}|e_t,η_{t+1})\n", idx);
}

//================================================================================
// Dec3核心：4次元後進確率計算 CalcPBE4D
//================================================================================
void SLFBAdec::CalcPBE4D(int idx, int Nb2) {
  assert(idx >= 0 && idx < Ns);

  // 現在の状態 t の確率を初期化
  for (int d0 = Dmin; d0 <= Dmax; d0++)
    for (int e0 = 0; e0 < NUM_ERROR_STATES; e0++)
      for (int k0 = 0; k0 < num_kmers; k0++)
        PBE4D[idx][d0-Dmin][e0][k0] = 0.0;

  // ✅ 修正：CalcPFE4Dと同じ効率的なループ構造
  for (int d0 = Dmin; d0 <= Dmax; d0++) {
    for (int e0 = 0; e0 < NUM_ERROR_STATES; e0++) {
      for (int k0 = 0; k0 < num_kmers; k0++) {
        for (int d1 = Dmin; d1 <= Dmax; d1++) {

          int Nu2 = Nu + d1 - d0;
          int iL = idx * Nu + d0;
          if ((Nu2 < 0) || (Nu2 > 2 * Nu) || (iL < 0) || (iL + Nu2 > Nb2)) continue;

          double s_per_codeword[Q];
          for(int xi=0; xi < Q; xi++) {
            s_per_codeword[xi] = ComputeObservationProbabilityFromDNA(idx, Nu2, iL, xi, k0, e0, Nb2);
          }

          for (int xi = 0; xi < Q; xi++) {
            int k1 = ComputeNextKmer(k0, xi);
            for (int e1 = 0; e1 < NUM_ERROR_STATES; e1++) {
              if (PBE4D[idx+1][d1-Dmin][e1][k1] == 0.0) continue;

              double kmer_error_prob = GetKmerErrorProb(e0, k1, e1);
              double gamma_xi = s_per_codeword[xi] * PU[idx][xi] * kmer_error_prob;

              PBE4D[idx][d0-Dmin][e0][k0] += PBE4D[idx+1][d1-Dmin][e1][k1] * gamma_xi;
            }
          }
        }
      }
    }
  }

  // 4次元正規化
  double total_sum = 0.0;
  for(int d0 = Dmin; d0 <= Dmax; d0++) {
    for(int e0 = 0; e0 < NUM_ERROR_STATES; e0++) {
      for(int k0 = 0; k0 < num_kmers; k0++) {
        total_sum += PBE4D[idx][d0-Dmin][e0][k0];
      }
    }
  }
  if(total_sum > 0.0) {
    for(int d0 = Dmin; d0 <= Dmax; d0++) {
      for(int e0 = 0; e0 < NUM_ERROR_STATES; e0++) {
        for(int k0 = 0; k0 < num_kmers; k0++) {
          PBE4D[idx][d0-Dmin][e0][k0] /= total_sum;
        }
      }
    }
  }
  
  printf("# Dec3: CalcPBE4D[%d] completed with P(e_{t+1}|e_t,η_{t+1})\n", idx);
}

//================================================================================
// Dec3核心：4次元事後確率計算 CalcPDE4D
// 4次元前進・後進確率を使って符号語の事後確率を計算
//================================================================================
void SLFBAdec::CalcPDE4D(int idx, int Nb2) {
  assert(idx >= 0 && idx < Ns);
  
  // ✅ Dec3効率化修正：8重ネストループを劇的に簡略化
  for (int xi = 0; xi < Q; xi++) {
    PD[idx][xi] = 0.0;
    
    for (int d0 = Dmin; d0 <= Dmax; d0++) {
      for (int e0 = 0; e0 < NUM_ERROR_STATES; e0++) {
        for (int k0 = 0; k0 < num_kmers; k0++) {

          if (PFE4D[idx][d0-Dmin][e0][k0] == 0.0) continue;

          for (int d1 = Dmin; d1 <= Dmax; d1++) {
            int Nu2 = Nu + d1 - d0;
            int iL = idx * Nu + d0;
            if ((Nu2 < 0) || (Nu2 > 2 * Nu) || (iL < 0) || (iL + Nu2 > Nb2)) continue;
            
            // ✅ Dec3新方式: DNAチャネル統合動的確率計算（CalcPDE4D用）
            double obs_prob = ComputeObservationProbabilityFromDNA(idx, Nu2, iL, xi, k0, e0, Nb2);

            int k1 = ComputeNextKmer(k0, xi); // 決定論的遷移: k0 + xi → k1

            for (int e1 = 0; e1 < NUM_ERROR_STATES; e1++) {
              double kmer_error_prob = GetKmerErrorProb(e0, k1, e1);
              double posterior_contrib = PFE4D[idx][d0-Dmin][e0][k0] * obs_prob * PU[idx][xi] * kmer_error_prob * PBE4D[idx+1][d1-Dmin][e1][k1];
              
              PD[idx][xi] += posterior_contrib;
            }
          }
        }
      }
    }
  }
  
  // 正規化（符号語確率の合計を1にする）
  normalize(PD[idx], Q);
  
  printf("# Dec3: CalcPDE4D[%d] completed with 4D marginalization\n", idx);
}

