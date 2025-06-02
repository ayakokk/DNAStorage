#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

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
long SLFBAdec::VectToLong(const unsigned char *V, int len){
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
void SLFBAdec::LongToVect(unsigned char *V, long val, int len){
  assert(len>0 && val>=0);
  long mask = 0x1 << (len-1);
  for(int i=0;i<len;i++){
    V[i] = ( (val & mask)==0 )? 0 : 1;
    mask >>= 1;
  } // for i
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
	  printf("%e ly=%d y=%ld xi=%d x=%ld ",GX[ly][y][xi],ly,y,xi,x);
	  PrintVect(X,Nu,"X=","\n");
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
  double ret,qt;
  unsigned char *X = new unsigned char [lx];
  unsigned char *Y = new unsigned char [ly];
  LongToVect(X,x,lx);
  LongToVect(Y,y,ly);
  if(lx==1){
    if(ly==1){
      ret = Psub(X[0],Y[0]);
    } else if(ly==2) {
      ret = Psub(X[0],Y[0]) * Psub(X[0],Y[1]);
    } else {
      ret = 0.0;
    } // if ly
  } else {
    x1  = VectToLong(&X[1],lx-1);
    // -----trans
    qt = Psub(X[0],Y[0]);
    y1 = (ly==1)? 0 : VectToLong(&Y[1],ly-1);
    qt *= CalcPyx(y1,x1,ly-1,lx-1);
    // -----ins
    /*
    if(ly<2){
      qi = 0.0;
    } else {
      qi = Psub(X[0],Y[0]) * Psub(X[0],Y[1]);
      if(ly>2){
	y1  = VectToLong(&Y[2],ly-2);
	qi *= CalcPyx(y1,x1,ly-2,lx-1);	
      } // if ly>2
    } // if ly<2
    */
    // -----del

    // ---
    ret = qt; //TMP
  } // if lx
  delete [] X;
  delete [] Y;
  //printf("ret=%e\n",ret);
  return ret;
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
  Nu2p = (long) pow(2,Nu);
  Ns   = Nb/Nu;
  Drng = Dmax-Dmin+1;
  Nu2max = Nu + (int)ceil( (double)Nu*Pi ) + 2;
  Nu2min = Nu - (int)ceil( (double)Nu*Pd ) - 2;
  Nu2max = min(Nu2max, Nu*2);
  Nu2min = max(Nu2min, 0   );
  printf("# SLFBAdec: Ns=%d (Nu;Nu2min,Nu2max)=(%d(%ld);%d,%d) Nb=%d Q=%d\n",Ns,Nu,Nu2p,Nu2min,Nu2max,Nb,Q);
  printf("# SLFBAdec: (Pi,Pd,Ps)=(%e,%e,%e) (Dmin,Dmax:Drng)=(%d,%d:%d)\n",
	 Pi,Pd,Ps,Dmin,Dmax,Drng);
  assert(Nb%Nu==0);
  assert(Nu<=NuMax);
  assert(ECM->GetM()==Q && ECM->GetN()==Q);
  //----- set tables
  SetGD();
  SetGX();
}

//================================================================================
SLFBAdec::~SLFBAdec(){
  DelGD();
  DelGX();
  printf("# SLFBAdec: deleted\n");
}

