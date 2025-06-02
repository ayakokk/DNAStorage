#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "InnerCodebook.hpp"
#include "ChannelMatrix.hpp"
#include "IDSchannel.hpp"
#include "CSFBAdec.hpp"

#define NuMax 10

//================================================================================
int CSFBAdec::max(int a, int b){return (a>b)? a : b;}
int CSFBAdec::min(int a, int b){return (a<b)? a : b;}

//================================================================================
int CSFBAdec::max(const int *a, int len){
  int p=0;
  for(int i=1;i<len;i++){
    if(a[i]>a[p]) p=i;
  } // for i
  return a[p];
}

//================================================================================
int CSFBAdec::min(const int *a, int len){
  int p=0;
  for(int i=1;i<len;i++){
    if(a[i]<a[p]) p=i;
  } // for i
  return a[p];
}

//================================================================================
double CSFBAdec::Psub(unsigned char a, unsigned char b){
  assert( a==0 || a==1 );
  assert( b==0 || b==1 );
  return (a==b)? 1.0-Ps : Ps;
}

//================================================================================
void CSFBAdec::ClearMat(double **X, int M0, int N0){
  for(int i=0;i<M0;i++){
    for(int j=0;j<N0;j++) X[i][j]=0.0;
  } // for i
}

//================================================================================
void CSFBAdec::ClearMat(double ***X, int M0, int N0, int K0){
  for(int i=0;i<M0;i++){
    for(int j=0;j<N0;j++)
      for(int k=0;k<K0;k++) X[i][j][k]=0.0;
  } // for i
}

//================================================================================
void CSFBAdec::SetMat(double **X, int M0, int N0, double val){
  for(int i=0;i<M0;i++){
    for(int j=0;j<N0;j++) X[i][j]=val;
  } // for i
}

//================================================================================
void CSFBAdec::CopyMat(double **Y, const double **X, int M0, int N0){
  assert(M0>0 && N0>0);
  for(int i=0;i<M0;i++) memcpy(Y[i],X[i],sizeof(double)*N0);
}

//================================================================================
void CSFBAdec::PrintMat(const double **X, int M0, int N0, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<M0;i++){
    printf("%03d: ",i);
    for(int j=0;j<N0;j++) printf("%.2e ",X[i][j]);
    printf("\n");
  } // for i
  printf("%s",post);
}

//================================================================================
void CSFBAdec::MultMat(double **Y, const double **X0, const double **X1, int M0, int N0, int M1, int N1){
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
void CSFBAdec::PrintVect(const unsigned char *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%u",V[i]);
  printf("%s",post);
}

//================================================================================
void CSFBAdec::PrintVect(const int *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d ",V[i]);
  printf("%s",post);
}

//================================================================================
void CSFBAdec::PrintVect(const double *V, int len, const char *pre, const char *post){
  int w=20;
  printf("%s",pre);
  for(int i=0;i<len;i++){
    printf("%.2e ",V[i]);
    if(i%w==w-1) printf("\n");
  } // for i
  printf("%s",post);
}

//================================================================================
long CSFBAdec::VectToLong(const unsigned char *V, int len){
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
void CSFBAdec::LongToVect(unsigned char *V, long val, int len){
  assert(len>0 && val>=0);
  long mask = 0x1 << (len-1);
  for(int i=0;i<len;i++){
    V[i] = ( (val & mask)==0 )? 0 : 1;
    mask >>= 1;
  } // for i
}

//--------------------------------------------------------------------------------
void CSFBAdec::normalize(double *V, int len){
  assert(len>0);
  double s=0;
  for(int i=0;i<len;i++){
    assert( V[i]>=0.0 );
    s += V[i];
  } // for i
  if( s>0 ){
    for(int i=0;i<len;i++) V[i] /= s;
  } else {
    for(int i=0;i<len;i++) V[i] = 1.0/(double)len; //??
  } // if
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

//================================================================================
void CSFBAdec::SetGD(){
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
void CSFBAdec::DelGD(){
  for(int i=0;i<Drng;i++) delete [] GD[i];
  delete [] GD;
}

//================================================================================
void CSFBAdec::SetGX(){
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
void CSFBAdec::DelGX(){
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
double CSFBAdec::CalcPyx(long y, long x, int ly, int lx){
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
double CSFBAdec::GetGX(int Nu2, long y, long xi){
  assert(Nu2>=0 && Nu2<=2*Nu);
  assert(y >=0);
  assert(xi>=0 && xi<Q);
  if(Nu2<Nu2min || Nu2>Nu2max) return GX[Nu2][0][0 ]; // approx
  else                         return GX[Nu2][y][xi]; 
}

//================================================================================
void CSFBAdec::SetFG(){
  PF = new double ** [Nseq];
  PB = new double ** [Nseq];
  PU = new double ** [Nseq];
  PD = new double ** [Nseq];
  Yin= new unsigned char * [Nseq];
  for(int j=0;j<Nseq;j++){
    PF[j] = new double * [Ns+1];
    PB[j] = new double * [Ns+1];
    for(int i=0;i<Ns+1;i++){
      PF[j][i] = new double [Drng];
      PB[j][i] = new double [Drng];
    } // for i
    PU[j] = new double * [Ns];
    PD[j] = new double * [Ns];
    for(int i=0;i<Ns;i++){
      PU[j][i] = new double [Q];
      PD[j][i] = new double [Q];
    } // for i
    Yin[j] = new unsigned char [Nb+Dmax];
  } // for j
  PI = new double * [Ns];
  PO = new double * [Ns];
  PU0= new double * [Ns];
  PM = new double * [Ns];
  for(int i=0;i<Ns;i++){
    PI[i]  = new double [Q];
    PO[i]  = new double [Q];
    PU0[i] = new double [Q];
    PM[i]  = new double [Q];
  } // for i
}

//================================================================================
void CSFBAdec::DelFG(){
  for(int j=0;j<Nseq;j++){
    for(int i=0;i<Ns+1;i++){
      delete [] PF[j][i];
      delete [] PB[j][i];
    } // for i
    delete [] PF[j];
    delete [] PB[j];
    for(int i=0;i<Ns;i++){
      delete [] PU[j][i];
      delete [] PD[j][i];
    } // for i
    delete [] PU[j];
    delete [] PD[j];
    delete [] Yin[j];
  } // for j
  for(int i=0;i<Ns;i++){
    delete [] PI[i];
    delete [] PO[i];
    delete [] PU0[i];
    delete [] PM[i];
  } // for i
  delete [] PU;
  delete [] PD;
  delete [] PI;
  delete [] PO;
  delete [] PU0;
  delete [] PM;
  delete [] Yin;
}

//================================================================================
void CSFBAdec::InitFG(const unsigned char **RW, const double **Pin, int *Nb2){
  ClearMat(PF, Nseq, Ns+1, Drng);
  ClearMat(PB, Nseq, Ns+1, Drng);
  ClearMat(PU, Nseq, Ns,   Q);
  ClearMat(PD, Nseq, Ns,   Q);
  ClearMat(PO,       Ns,   Q);
  SetMat(  PM,       Ns,   Q, 1.0);
  CopyMat(PI, Pin, Ns, Q);
  for(int j=0;j<Nseq;j++){
    PF[j][0 ][0        -Dmin] = 1.0;
    PB[j][Ns][Nb2[j]-Nb-Dmin] = 1.0;
    memcpy(Yin[j],RW[j],sizeof(unsigned char)*Nb2[j]);
  } // for j
}

//================================================================================
void CSFBAdec::CalcPU0(int idx){
  assert(idx>=0 && idx<Ns);
  int Vin,Vup;
  for(Vup=0;Vup<Q;Vup++){
    PU0[idx][Vup] = 0.0;
    for(Vin=0;Vin<Q;Vin++){
      PU0[idx][Vup] += (PI[idx][Vin] * ECM->GetPxy(Vin,Vup));
    } // for Vin
  } // for Vup
  normalize(PU0[idx], Q);
}

//================================================================================
void CSFBAdec::CalcPU(int idx){
  assert(idx>=0 && idx<Ns);
  double *Ptmp = new double [Q];
  for(int xi=0;xi<Q;xi++) Ptmp[xi] = PU0[idx][xi]*PM[idx][xi];
  for(int jdx=0;jdx<Nseq;jdx++){
    for(int xi=0;xi<Q;xi++) Ptmp[xi] *= PD[jdx][idx][xi];
  } // for jdx
  for(int jdx=0;jdx<Nseq;jdx++){
    for(int xi=0;xi<Q;xi++){
      if(PD[jdx][idx][xi]==0.0) PU[jdx][idx][xi] = 0.0; 
      else                      PU[jdx][idx][xi] = Ptmp[xi]/PD[jdx][idx][xi];
    } // for xi
    normalize(PU[jdx][idx],Q);
  } // for jdx
  delete [] Ptmp;
}

//================================================================================
void CSFBAdec::CalcPM(int idx){
  assert(idx>=0 && idx<Ns);
  for(int xi=0;xi<Q;xi++) PM[idx][xi] = 1.0;
  for(int jdx=0;jdx<Nseq;jdx++){
    for(int xi=0;xi<Q;xi++) PM[idx][xi] *= PD[jdx][idx][xi];
  } // for jdx
  normalize(PM[idx],Q);
}

//================================================================================
void CSFBAdec::CalcPO(int idx){
  assert(idx>=0 && idx<Ns);
  int Vout,Vdn;
  double *Ptmp = new double [Q];
  for(int xi=0;xi<Q;xi++) Ptmp[xi] = 1.0;
  for(int jdx=0;jdx<Nseq;jdx++){
    for(int xi=0;xi<Q;xi++) Ptmp[xi] *= PD[jdx][idx][xi];
  } // for jdx
  for(Vout=0;Vout<Q;Vout++){
    PO[idx][Vout] = 0.0;
    for(Vdn=0;Vdn<Q;Vdn++){
      PO[idx][Vout] += (Ptmp[Vdn] * ECM->GetPxy(Vout,Vdn));
    } // for Vdn
  } // for Vout
  normalize(PO[idx], Q);
  delete [] Ptmp;
}

//================================================================================
void CSFBAdec::CalcPF(int jdx, int idx, int Nb2){
  assert(idx>=0 && idx<Ns);
  assert(jdx>=0 && jdx<Nseq);
  int  d0,d1;
  int  xi;
  int  Nu2, iL, iR;
  long y;
  double s,ss;
  for(d1=Dmin;d1<=Dmax;d1++){
    PF[jdx][idx+1][d1-Dmin] = 0.0;
    for(d0=Dmin;d0<=Dmax;d0++){
      Nu2 = Nu+d1-d0;
      iL  = idx    *Nu + d0;
      iR  = (idx+1)*Nu + d1 - 1;
      if( (Nu2<0) || (Nu2>2*Nu) || (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
      if( Nu2>=Nu2min && Nu2<=Nu2max ){
	y = VectToLong( &Yin[jdx][iL], Nu2 );
	s = 0.0;
	for(xi=0;xi<Q;xi++){
	  ss = PU[jdx][idx][xi] * GetGX(Nu2,y,xi);
	  s += ss;
	  //(dbg)
	  //printf("%e d1=%d d0=%d Nu2=%d iL=%d iR=%d y=%ld xi=%d x=%ld\n",ss,d1,d0,Nu2,iL,iR,y,xi,x);
	} // for xi
      } else {
	s = GetGX(Nu2,0,0);
      } // if Nu2
      PF[jdx][idx+1][d1-Dmin] += (s*PF[jdx][idx][d0-Dmin]);
    } // for d0
  } // for d1
  normalize(PF[jdx][idx+1],Drng);
}

//================================================================================
void CSFBAdec::CalcPB(int jdx, int idx, int Nb2){
  assert(idx>=0 && idx<Ns);
  assert(jdx>=0 && jdx<Nseq);
  int  d0,d1;
  int  xi;
  int  Nu2, iL, iR;
  long y;
  double s,ss;
  for(d0=Dmin;d0<=Dmax;d0++){
    PB[jdx][idx][d0-Dmin] = 0.0;
    for(d1=Dmin;d1<=Dmax;d1++){
      Nu2 = Nu+d1-d0;
      iL  = idx    *Nu + d0;
      iR  = (idx+1)*Nu + d1 - 1;
      if( (Nu2<0) || (Nu2>2*Nu) || (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
      if( Nu2>=Nu2min && Nu2<=Nu2max ){
	y = VectToLong( &Yin[jdx][iL], Nu2 );
	s = 0.0;
	for(xi=0;xi<Q;xi++){
	  ss = PU[jdx][idx][xi] * GetGX(Nu2,y,xi);
	  s += ss;
	  //(dbg)
	  //printf("%e d1=%d d0=%d Nu2=%d iL=%d iR=%d y=%ld xi=%d x=%ld\n",ss,d1,d0,Nu2,iL,iR,y,xi,x);
	} // for xi
      } else {
	s = GetGX(Nu2,0,0);
      } // if Nu2
      PB[jdx][idx][d0-Dmin] += (s*PB[jdx][idx+1][d1-Dmin]);
    } // for d1
  } // for d0
  normalize(PB[jdx][idx],Drng);
}

//================================================================================
void CSFBAdec::CalcPD(int jdx, int idx, int Nb2){
  assert(jdx>=0 && jdx<Nseq);
  assert(idx>=0 && idx<Ns);
  int  d0,d1;
  int  xi;
  int  Nu2, iL, iR;
  long y;
  double s,ss;
  for(xi=0;xi<Q;xi++){
    PD[jdx][idx][xi] = 0.0;
    for(d1=Dmin;d1<=Dmax;d1++){
      s = 0.0;
      for(d0=Dmin;d0<=Dmax;d0++){
	Nu2 = Nu+d1-d0;
	iL  = idx    *Nu + d0;
	iR  = (idx+1)*Nu + d1 - 1;
	if( (Nu2<0) || (Nu2>2*Nu) || (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
	y = VectToLong( &Yin[jdx][iL], Nu2 );
	ss = PF[jdx][idx][d0-Dmin] * GetGX(Nu2,y,xi);
	s += ss;
	//(dbg)
	//printf("%e d1=%d d0=%d Nu2=%d iL=%d iR=%d y=%ld xi=%d x=%ld\n",ss,d1,d0,Nu2,iL,iR,y,xi,x);
      } // for d0
      PD[jdx][idx][xi] += (s*PB[jdx][idx+1][d1-Dmin]);
    } // for d1
  } // for xi
  normalize(PD[jdx][idx],Q);
}

//================================================================================
void CSFBAdec::CalcPDf(int jdx, int idx, int Nb2){
  assert(jdx>=0 && jdx<Nseq);
  assert(idx>=0 && idx<Ns);
  int  d0,d1;
  int  xi;
  int  Nu2, iL, iR;
  long y;
  double s,ss;
  for(xi=0;xi<Q;xi++){
    PD[jdx][idx][xi] = 0.0;
    for(d0=Dmin;d0<=Dmax;d0++){
      s = 0.0;
      for(d1=Dmin;d1<=Dmax;d1++){
	Nu2 = Nu+d1-d0;
	iL  = idx    *Nu + d0;
	iR  = (idx+1)*Nu + d1 - 1;
	if( (Nu2<0) || (Nu2>2*Nu) || (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
	y = VectToLong( &Yin[jdx][iL], Nu2 );
	ss = GetGX(Nu2,y,xi);
	s += ss;
	//(dbg)
	//printf("%e d1=%d d0=%d Nu2=%d iL=%d iR=%d y=%ld xi=%d x=%ld\n",ss,d1,d0,Nu2,iL,iR,y,xi,x);
      } // for d1
      PD[jdx][idx][xi] += (s*PF[jdx][idx][d0-Dmin]); 
    } // for d0
  } // for xi
  normalize(PD[jdx][idx],Q);
}

//================================================================================
void CSFBAdec::CalcPDb(int jdx, int idx, int Nb2){
  assert(jdx>=0 && jdx<Nseq);
  assert(idx>=0 && idx<Ns);
  int  d0,d1;
  int  xi;
  int  Nu2, iL, iR;
  long y;
  double s,ss;
  for(xi=0;xi<Q;xi++){
    PD[jdx][idx][xi] = 0.0;
    for(d1=Dmin;d1<=Dmax;d1++){
      s = 0.0;
      for(d0=Dmin;d0<=Dmax;d0++){
	Nu2 = Nu+d1-d0;
	iL  = idx    *Nu + d0;
	iR  = (idx+1)*Nu + d1 - 1;
	if( (Nu2<0) || (Nu2>2*Nu) || (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
	y = VectToLong( &Yin[jdx][iL], Nu2 );
	ss = GetGX(Nu2,y,xi);
	s += ss;
	//(dbg)
	//printf("%e d1=%d d0=%d Nu2=%d iL=%d iR=%d y=%ld xi=%d x=%ld\n",ss,d1,d0,Nu2,iL,iR,y,xi,x);
      } // for d0
      PD[jdx][idx][xi] += (s*PB[jdx][idx+1][d1-Dmin]);
    } // for d1
  } // for xi
  normalize(PD[jdx][idx],Q);
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
CSFBAdec::CSFBAdec(class InnerCodebook *_ICB, class ChannelMatrix *_ECM, class IDSchannel *_CH, int _Nseq, int _NumItr){
  ICB  = _ICB;
  ECM  = _ECM;
  CH   = _CH;
  Nseq = _Nseq;
  NumItr = _NumItr;
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
  printf("# CSFBAdec: Ns=%d (Nu;Nu2min,Nu2max)=(%d(%ld);%d,%d) Nb=%d Q=%d Nseq=%d\n",Ns,Nu,Nu2p,Nu2min,Nu2max,Nb,Q,Nseq);
  printf("# CSFBAdec: (Pi,Pd,Ps,Pt)=(%e,%e,%e,%e) (Dmin,Dmax:Drng)=(%d,%d:%d)\n",
	 Pi,Pd,Ps,Pt,Dmin,Dmax,Drng);
  assert(Nb%Nu==0);
  assert(Nu<=NuMax);
  assert(ECM->GetM()==Q && ECM->GetN()==Q);
  assert(Pt>0.0 && Pt<=1.0);
  //----- set tables
  printf("# CSFBAdec: generating GD & GX\n");
  SetGD();
  SetGX();
  printf("# CSFBAdec: generating factor graph\n");
  SetFG();
}

//================================================================================
CSFBAdec::~CSFBAdec(){
  DelGD();
  DelGX();
  DelFG();
  printf("# CSFBAdec: deleted\n");
}

//================================================================================
void CSFBAdec::Decode(double **Pout, const unsigned char **RW, int *Nb2, const int *dbgIW, const double **Pin){
  int idx, jdx, itr;
  int Nb2max = max(Nb2,Nseq);
  int Nb2min = min(Nb2,Nseq);
  assert( Nb2min>=Nb+Dmin && Nb2max<=Nb+Dmax );
  InitFG(RW,Pin,Nb2);
  // ----- init
  //printf("  Init\n");
  for(idx=0;idx<Ns;idx++) CalcPU0(idx);

  for(itr=0; itr<NumItr; itr++){
    // ----- forward
    for(idx=0; idx<Ns; idx++){
      //printf("  F: itr=%d idx=%d\n",itr,idx);
      for(jdx=0; jdx<Nseq; jdx++) CalcPDf(jdx,idx,Nb2[jdx]);
      CalcPU(idx);
      for(jdx=0; jdx<Nseq; jdx++) CalcPF( jdx,idx,Nb2[jdx]);
    } // for idx
    // ----- backward
    for(idx=Ns-1;idx>=0;idx--){
      //printf("  B: itr=%d idx=%d\n",itr,idx);
      for(jdx=0; jdx<Nseq; jdx++) CalcPDb(jdx,idx,Nb2[jdx]);
      CalcPU(idx);
      for(jdx=0; jdx<Nseq; jdx++) CalcPB( jdx,idx,Nb2[jdx]);
    } // for idx
    // ----- store
    for(idx=0; idx<Ns; idx++){
      //printf("  S: itr=%d idx=%d\n",itr,idx);
      for(jdx=0; jdx<Nseq; jdx++) CalcPD(jdx,idx,Nb2[jdx]);
      CalcPM(idx);
    } // for idx
  } // for itr
  
  // ----- output
  //printf("  Output\n");
  for(idx=0; idx<Ns; idx++) CalcPO(idx);
  CopyMat(Pout,(const double**)PO,Ns,Q);
  
  /*
  for(idx=0;   idx<Ns;idx++) CalcPF(idx,Nb2);
  for(idx=Ns-1;idx>=0;idx--) CalcPB(idx,Nb2);
  for(idx=0;   idx<Ns;idx++) CalcPD(idx,Nb2);
  for(idx=0;   idx<Ns;idx++) CalcPO(idx);
  CopyMat(Pout,(const double**)PO,Ns,Q);
  */
  
  // PrintNode(dbgIW);
  // PrintNode(0,0,   dbgIW[0]);
  // PrintNode(0,1,   dbgIW[1]);
  // PrintNode(0,2,   dbgIW[2]);
  // PrintNode(0,Ns-4,dbgIW[Ns-2]);
  // PrintNode(0,Ns-3,dbgIW[Ns-2]);
  // PrintNode(0,Ns-2,dbgIW[Ns-2]);
  // PrintNode(0,Ns-1,dbgIW[Ns-1]);
}

//================================================================================
void CSFBAdec::Decode(double **Pout, const unsigned char **RW, int *Nb2, const int *dbgIW){
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
void CSFBAdec::PrintNode(int jdx, int idx, int iw){
  assert(jdx>=0 && jdx<Nseq);
  assert(idx>=0 && idx<Ns);
  printf("[%04d:%04d] %02d\n",jdx,idx,iw);
  PrintVect(PF[jdx][idx  ],Drng,"PF:  ","\n");
  PrintVect(PB[jdx][idx+1],Drng,"PB:  ","\n");
  PrintVect(PU[jdx][idx  ],Q,   "PU:  ","\n");
  PrintVect(PD[jdx][idx  ],Q,   "PD:  ","\n");
  PrintVect(PI[idx  ],     Q,   "PI:  ","\n");
  PrintVect(PU0[idx ],     Q,   "PU0: ","\n");
  PrintVect(PM[idx  ],     Q,   "PM:  ","\n");
  PrintVect(PO[idx  ],     Q,   "PO:  ","\n");
}

//================================================================================
void CSFBAdec::PrintNode(const int *dbgIW){
  for(int i=0;i<Ns;i++)
    for(int j=0;j<Nseq;j++) PrintNode(j,i,dbgIW[i]);
}

