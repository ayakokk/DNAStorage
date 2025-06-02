#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "CPScodec.hpp"

#define BSIZE 4096

//================================================================================
int CPScodec::max(int a, int b){return (a>b)? a : b;}
int CPScodec::min(int a, int b){return (a<b)? a : b;}

//================================================================================
int CPScodec::max(const int *a, int len){
  int p=0;
  for(int i=1;i<len;i++){
    if(a[i]>a[p]) p=i;
  } // for i
  return a[p];
}

//================================================================================
int CPScodec::min(const int *a, int len){
  int p=0;
  for(int i=1;i<len;i++){
    if(a[i]<a[p]) p=i;
  } // for i
  return a[p];
}

//================================================================================
double CPScodec::Psub(unsigned char a, unsigned char b){
  assert( a==0 || a==1 );
  assert( b==0 || b==1 );
  return (a==b)? 1.0-Ps : Ps;
}

//================================================================================
void CPScodec::ClearVect(double *X, int M0){
  for(int i=0; i<M0; i++) X[i] = 0.0;
}

//================================================================================
void CPScodec::ClearMat(double **X, int M0, int N0){
  for(int i=0;i<M0;i++){
    for(int j=0;j<N0;j++) X[i][j]=0.0;
  } // for i
}

//================================================================================
void CPScodec::ClearMat(double ***X, int M0, int N0, int K0){
  for(int i=0;i<M0;i++){
    for(int j=0;j<N0;j++)
      for(int k=0;k<K0;k++) X[i][j][k]=0.0;
  } // for i
}

//================================================================================
void CPScodec::SetMat(double **X, int M0, int N0, double val){
  for(int i=0;i<M0;i++){
    for(int j=0;j<N0;j++) X[i][j]=val;
  } // for i
}

//================================================================================
void CPScodec::CopyMat(double **Y, const double **X, int M0, int N0){
  assert(M0>0 && N0>0);
  for(int i=0;i<M0;i++) memcpy(Y[i],X[i],sizeof(double)*N0);
}

//================================================================================
void CPScodec::PrintMat(const double **X, int M0, int N0, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<M0;i++){
    printf("%03d: ",i);
    for(int j=0;j<N0;j++) printf("%.2e ",X[i][j]);
    printf("\n");
  } // for i
  printf("%s",post);
}

//================================================================================
void CPScodec::MultMat(double **Y, const double **X0, const double **X1, int M0, int N0, int M1, int N1){
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
void CPScodec::PrintVect(const unsigned char *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%u",V[i]);
  printf("%s",post);
}

//================================================================================
void CPScodec::PrintVect(const int *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d ",V[i]);
  printf("%s",post);
}

//================================================================================
void CPScodec::PrintVect(const double *V, int len, const char *pre, const char *post){
  int w=20;
  printf("%s",pre);
  for(int i=0;i<len;i++){
    printf("%.2e ",V[i]);
    if(i%w==w-1) printf("\n");
  } // for i
  printf("%s",post);
}

//================================================================================
long CPScodec::VectToLong(const unsigned char *V, int len){
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
void CPScodec::LongToVect(unsigned char *V, long val, int len){
  assert(len>0 && val>=0);
  long mask = 0x1 << (len-1);
  for(int i=0;i<len;i++){
    V[i] = ( (val & mask)==0 )? 0 : 1;
    mask >>= 1;
  } // for i
}

//================================================================================
void CPScodec::normalize(double *V, int len){
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
void CPScodec::SetCSL(const char *fn){
  FILE *fp;
  char *buf = new char [BSIZE];
  if((fp=fopen(fn,"r"))==NULL){
    fprintf(stderr,"Cannot open %s\n",fn);
    exit(1);
  } // if
  
  // L1: Kc QB
  assert( fgets(buf,BSIZE,fp)!=NULL );
  Kc  = atoi(strtok(buf,  " \t"  ));
  QB  = atoi(strtok(NULL, " \t\n"));
  QB2 = QB/2; 
  printf("# CPScodec::SetCSL Kc=%d QB=%d QB2=%d\n",Kc,QB,QB2);
  assert(Kc>0);
  assert(QB>0 && QB<=2*(Kc+1) && QB%2==0);
  CSL = new int * [QB];
  for(int i=0;i<QB;i++) CSL[i] = new int [4];
  // 1st half: b=0
  for(int i=0;i<QB2;i++){
    assert( fgets(buf,BSIZE,fp)!=NULL );
    CSL[i][0] = atoi(strtok(buf,  " \t"));
    CSL[i][1] = atoi(strtok(NULL, " \t"));
    CSL[i][2] = atoi(strtok(NULL, " \t"));
    CSL[i][3] = atoi(strtok(NULL, " \t\n"));
    assert(CSL[i][0] + CSL[i][1] == Kc);
    assert(CSL[i][2]==0 && CSL[i][3]==0);
    for(int j=0;j<4;j++) assert( CSL[i][j]>=0 && CSL[i][j]<=Kc );
  } // for i
  // 2nd half: b=1
  for(int i=QB2;i<QB;i++){
    assert( fgets(buf,BSIZE,fp)!=NULL );
    CSL[i][0] = atoi(strtok(buf,  " \t"));
    CSL[i][1] = atoi(strtok(NULL, " \t"));
    CSL[i][2] = atoi(strtok(NULL, " \t"));
    CSL[i][3] = atoi(strtok(NULL, " \t\n"));
    assert(CSL[i][2] + CSL[i][3] == Kc);
    assert(CSL[i][0]==0 && CSL[i][1]==0);
    for(int j=0;j<4;j++) assert( CSL[i][j]>=0 && CSL[i][j]<=Kc );
  } // for i
  
  fclose(fp);
  delete [] buf;
}

//================================================================================
void CPScodec::DelCSL(){
  for(int i=0;i<QB;i++) delete [] CSL[i];
  delete [] CSL;
}

//================================================================================
void CPScodec::SetGS(){
  int    y0, y1, s, x;
  double p;
  // ----- GS1
  GS1 = new double * [4];
  for(int i=0; i<4; i++) GS1[i] = new double [QB];
  for(y0=0; y0<4; y0++){
    for(s=0; s<QB; s++){
      GS1[y0][s] = 0.0;
      for(x=0; x<4; x++){
	p = (double)CSL[s][x]/Kc * Psub(x,y0);
	GS1[y0][s] += p;
      } // for x
      GS1[y0][s] *= Pt;
    } // for s
  } // for y0
  // ----- GS2
  GS2 = new double ** [4];
  for(int i=0;i<4;i++){
    GS2[i] = new double * [4];
    for(int j=0;j<4;j++) GS2[i][j] = new double [QB];
  } // for i
  for(y0=0; y0<4; y0++){
    for(y1=0; y1<4; y1++){
      for(s=0; s<QB; s++){
	GS2[y0][y1][s] = 0.0;
	for(x=0; x<4; x++){
	  p = (double)CSL[s][x]/Kc * Psub(x,y0) * Psub(x,y1);
	  GS2[y0][y1][s] += p;
	} // for x
	GS2[y0][y1][s] *= Pi;
      } // for s
    } // for y1
  } // for y0
}

//================================================================================
void CPScodec::DelGS(){
  for(int i=0;i<4;i++) delete [] GS1[i];
  delete [] GS1;
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++) delete [] GS2[i][j];
    delete [] GS2[i];
  } // for i
  delete [] GS2;
}

//================================================================================
void CPScodec::SetFG(){
  // PF,PB
  PF = new double ** [Nseq];
  PB = new double ** [Nseq];
  for(int i=0; i<Nseq; i++){
    PF[i] = new double * [Nb+1];
    PB[i] = new double * [Nb+1];
    for(int j=0; j<Nb+1; j++){
      PF[i][j] = new double [Drng];
      PB[i][j] = new double [Drng];
    } // for j
  } // for i
  // PU,PD
  PU = new double ** [Nseq];
  PD = new double ** [Nseq];
  for(int i=0; i<Nseq; i++){
    PU[i] = new double * [Nb];
    PD[i] = new double * [Nb];
    for(int j=0; j<Nb; j++){
      PU[i][j] = new double [QB];
      PD[i][j] = new double [QB];      
    } // for j
  } // for i
  // PM, PXU, PXD, PIB, PIC, POC
  PM  = new double * [Nb];
  PXU = new double * [Nb];
  PXD = new double * [Nb];
  PIB = new double * [Nb];
  PIC = new double * [Nb];
  POC = new double * [Nb];
  for(int i=0;i<Nb;i++){
    PM[i]  = new double [QB];
    PXU[i] = new double [QB];
    PXD[i] = new double [QB];
    PIB[i] = new double [2];
    PIC[i] = new double [QB2];
    POC[i] = new double [QB2];
  } // for i
}

//================================================================================
void CPScodec::DelFG(){
  // PF,PB
  for(int i=0; i<Nseq; i++){
    for(int j=0; j<Nb+1; j++){
      delete [] PF[i][j];
      delete [] PB[i][j];
    } // for j
    delete [] PF[i];
    delete [] PB[i];
  } // for i
  delete [] PF;
  delete [] PB;
  // PU,PD
  for(int i=0; i<Nseq; i++){
    for(int j=0; j<Nb; j++){
      delete [] PU[i][j];
      delete [] PD[i][j];
    } // for j
    delete [] PU[i];
    delete [] PD[i];
  } // for i
  delete [] PU;
  delete [] PD;
  // PM, PXU, PXD, PIB, PIC, POC
  for(int i=0; i<Nb; i++){
    delete [] PM[i];
    delete [] PXU[i];
    delete [] PXD[i];
    delete [] PIB[i];
    delete [] PIC[i];
    delete [] POC[i];
  } // for i
  delete [] PM;
  delete [] PXU;
  delete [] PXD;
  delete [] PIB;
  delete [] PIC;
  delete [] POC;
}

//================================================================================
double CPScodec::Psub(int x, int y){
  assert(x>=0 && x<=3);
  assert(y>=0 && y<=3);
  return (x==y)? 1.0-Ps : Ps/3.0;
}

//================================================================================
double CPScodec::Psbc(int s, int b, int c){
  assert( s>=0 && s<QB  );
  assert( b==0 || b==1  );
  assert( c>=0 && c<QB2 );
  if(b==0) return (s == c)?     1.0 : 0.0;
  else     return (s == c+QB2)? 1.0 : 0.0; 
}

//================================================================================
void CPScodec::CalcPXU(int idx){
  assert(idx>=0 && idx<Nb);
  int s, b, c;
  for(s=0; s<QB; s++){
    PXU[idx][s] = 0.0;
    for(b=0; b<2; b++){
      for(c=0; c<QB2; c++){
	PXU[idx][s] += (PIB[idx][b] * PIC[idx][c] * Psbc(s,b,c));
      } // for c
    } // for b
  } // for s
  normalize( PXU[idx], QB );
}

//================================================================================
void CPScodec::CalcPOC(int idx){
  assert(idx>=0 && idx<Nb);
  int s, b, c;
  for(c=0; c<QB2; c++){
    POC[idx][c] = 0.0;
    for(b=0; b<2; b++){
      for(s=0; s<QB; s++){
	POC[idx][c] += (PIB[idx][b] * PXD[idx][s] * Psbc(s,b,c));
      } // for s
    } // for b
  } // for c
  normalize( POC[idx], QB2 );  
}

//================================================================================
void CPScodec::CalcPDf(int jdx, int idx, int Nb2){
  assert( jdx>=0 && jdx<Nseq );
  assert( idx>=0 && idx<Nb   );
  int s, d0, d1;
  int iL, iR;
  double p;
  for(s=0; s<QB; s++){
    PD[jdx][idx][s] = 0.0;
    for(d0=Dmin; d0<=Dmax; d0++){
      for(d1=max(d0-1,Dmin); d1<=min(d0+1,Dmax); d1++){
	iL = idx + d0;
	iR = idx + d1;
	if( (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
	if(      d1 == d0+1 ) p = GS2[ Yin[jdx][idx+d0] ][ Yin[jdx][idx+d1] ][s];
	else if( d1 == d0-1 ) p = Pd;
	else if( d1 == d0   ) p = GS1[ Yin[jdx][idx+d0] ][s];
	else                  assert( false );
	PD[jdx][idx][s] += (p * PF[jdx][idx][d0-Dmin]);
	//(dbg)
	//printf("PDf: jdx=%d idx=%d s=%d d0=%d d1=%d iL=%d iR=%d p=%e\n",jdx,idx,s,d0,d1,iL,iR,p);
      } // for d1
    } // for d0
  } // for s
  normalize( PD[jdx][idx], QB );
}

//================================================================================
void CPScodec::CalcPDb(int jdx, int idx, int Nb2){
  assert( jdx>=0 && jdx<Nseq );
  assert( idx>=0 && idx<Nb   );
  int s, d0, d1;
  int iL, iR;
  double p;
  for(s=0; s<QB; s++){
    PD[jdx][idx][s] = 0.0;
    for(d0=Dmin; d0<=Dmax; d0++){
      for(d1=max(d0-1,Dmin); d1<=min(d0+1,Dmax); d1++){
	iL = idx + d0;
	iR = idx + d1;
	if( (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
	if(      d1 == d0+1 ) p = GS2[ Yin[jdx][idx+d0] ][ Yin[jdx][idx+d1] ][s];
	else if( d1 == d0-1 ) p = Pd;
	else if( d1 == d0   ) p = GS1[ Yin[jdx][idx+d0] ][s];
	else                  assert( false );
	PD[jdx][idx][s] += (p * PB[jdx][idx+1][d1-Dmin]);
	//(dbg)
	//printf("PDb: jdx=%d idx=%d s=%d d0=%d d1=%d iL=%d iR=%d p=%e\n",jdx,idx,s,d0,d1,iL,iR,p);
      } // for d1
    } // for d0
  } // for s
  normalize( PD[jdx][idx], QB );
}

//================================================================================
void CPScodec::CalcPD(int jdx, int idx, int Nb2){
  assert( jdx>=0 && jdx<Nseq );
  assert( idx>=0 && idx<Nb   );
  int s, d0, d1;
  int iL, iR;
  double p;
  for(s=0; s<QB; s++){
    PD[jdx][idx][s] = 0.0;
    for(d0=Dmin; d0<=Dmax; d0++){
      for(d1=max(d0-1,Dmin); d1<=min(d0+1,Dmax); d1++){
	iL = idx + d0;
	iR = idx + d1;
	if( (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
	if(      d1 == d0+1 ) p = GS2[ Yin[jdx][idx+d0] ][ Yin[jdx][idx+d1] ][s];
	else if( d1 == d0-1 ) p = Pd;
	else if( d1 == d0   ) p = GS1[ Yin[jdx][idx+d0] ][s];
	else                  assert( false );
	PD[jdx][idx][s] += (p * PF[jdx][idx][d0-Dmin] * PB[jdx][idx+1][d1-Dmin]);
	//(dbg)
	//printf("PD: jdx=%d idx=%d s=%d d0=%d d1=%d iL=%d iR=%d p=%e\n",jdx,idx,s,d0,d1,iL,iR,p);
      } // for d1
    } // for d0
  } // for s
  normalize( PD[jdx][idx], QB );
}

//================================================================================
void CPScodec::CalcPF(int jdx, int idx, int Nb2){
  assert( jdx>=0 && jdx<Nseq );
  assert( idx>=0 && idx<Nb   );
  int s, d0, d1;
  int iL, iR;
  double p;
  for(d1=Dmin; d1<=Dmax; d1++){
    PF[jdx][idx+1][d1-Dmin] = 0.0;
    for(d0=max(d1-1,Dmin); d0<=min(d1+1,Dmax); d0++){
      iL = idx + d0;
      iR = idx + d1;
      if( (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
      for(s=0; s<QB; s++){
	if(      d1 == d0+1 ) p = GS2[ Yin[jdx][idx+d0] ][ Yin[jdx][idx+d1] ][s];
	else if( d1 == d0-1 ) p = Pd;
	else if( d1 == d0   ) p = GS1[ Yin[jdx][idx+d0] ][s];
	else                  assert( false );
	PF[jdx][idx+1][d1-Dmin] += (p * PF[jdx][idx][d0-Dmin] * PU[jdx][idx][s]);
	//(dbg)
	//printf("PF: jdx=%d idx=%d s=%d d0=%d d1=%d iL=%d iR=%d p=%e\n",jdx,idx,s,d0,d1,iL,iR,p);
      } // for s 
    } // for d0
  } // for d1
  normalize( PF[jdx][idx+1], Drng );
}

//================================================================================
void CPScodec::CalcPB(int jdx, int idx, int Nb2){
  assert( jdx>=0 && jdx<Nseq );
  assert( idx>=0 && idx<Nb   );
  int s, d0, d1;
  int iL, iR;
  double p;
  for(d0=Dmin; d0<=Dmax; d0++){
    PB[jdx][idx][d0-Dmin] = 0.0;
    for(d1=max(d0-1,Dmin); d1<=min(d0+1,Dmax); d1++){
      iL = idx + d0;
      iR = idx + d1;
      if( (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
      for(s=0; s<QB; s++){
	if(      d1 == d0+1 ) p = GS2[ Yin[jdx][idx+d0] ][ Yin[jdx][idx+d1] ][s];
	else if( d1 == d0-1 ) p = Pd;
	else if( d1 == d0   ) p = GS1[ Yin[jdx][idx+d0] ][s];
	else                  assert( false );
	PB[jdx][idx][d0-Dmin] += (p * PB[jdx][idx+1][d1-Dmin] * PU[jdx][idx][s]);
	//(dbg)
	//printf("PB: jdx=%d idx=%d s=%d d0=%d d1=%d iL=%d iR=%d p=%e\n",jdx,idx,s,d0,d1,iL,iR,p);
      } // for s 
    } // for d1
  } // for d0
  normalize( PB[jdx][idx], Drng );
}

//================================================================================
void CPScodec::CalcPU(int idx){
  assert(idx>=0 && idx<Nb);
  double *Ptmp = new double [QB];
  memcpy(Ptmp, PXU[idx], sizeof(double)*QB);
  for(int jdx=0; jdx<Nseq; jdx++){
    for(int s=0; s<QB; s++) Ptmp[s] *= PD[jdx][idx][s];
  } // for jdx
  for(int jdx=0;jdx<Nseq;jdx++){
    for(int s=0; s<QB; s++){
      if(PD[jdx][idx][s]==0.0) PU[jdx][idx][s] = 0.0; 
      else                     PU[jdx][idx][s] = Ptmp[s]/PD[jdx][idx][s];
    } // for xi
    normalize(PU[jdx][idx],QB);
  } // for jdx
  delete [] Ptmp;
}

//================================================================================
void CPScodec::CalcPXD(int idx){
  assert(idx>=0 && idx<Nb);
  memcpy(PXD[idx], PD[0][idx], sizeof(double)*QB);
  for(int jdx=1; jdx<Nseq; jdx++){
    for(int s=0; s<QB; s++) PXD[idx][s] *= PD[jdx][idx][s];
  } // for jdx
  normalize(PXD[idx], QB);  
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
CPScodec::CPScodec(int _Nb, double _Pi, double _Pd, double _Ps, int _Dmin, int _Dmax,
		   int _Nu, int _Nseq, const char *fn){
  Nu   = _Nu;
  Nseq = _Nseq;
  Nb   = _Nb;
  Pi   = _Pi;
  Pd   = _Pd;
  Ps   = _Ps;  
  Dmin = _Dmin;
  Dmax = _Dmax;
  Pt   = 1.0-Pi-Pd;
  Ns   = Nb/Nu;
  Drng = Dmax-Dmin+1;
  printf("# CPScodec: Ns=%d Nu=%d Nb=%d Nseq=%d\n",Ns,Nu,Nb,Nseq);
  printf("# CPScodec: (Pi,Pd,Ps,Pt)=(%e,%e,%e,%e) (Dmin,Dmax:Drng)=(%d,%d:%d)\n",
	 Pi,Pd,Ps,Pt,Dmin,Dmax,Drng);
  assert(Nseq > 0);
  assert(Nu > 0 && Nb%Nu==0);
  assert(Pt>0.0 && Pt<=1.0);
  Yin = new unsigned char * [Nseq];
  for(int i=0; i<Nseq; i++) Yin[i] = new unsigned char [Nb+Dmax];
  //---
  SetCSL(fn);
  printf("# CPScodec: Kc=%d QB=%d QB2=%d\n", Kc,QB,QB2);
  assert(Kc>0 && QB>0 && QB%2==0);
  PrintCSL();
  SetGS();
  SetFG();
  //---
  //PrintGS();
}

//================================================================================
CPScodec::~CPScodec(){
  for(int i=0; i<Nseq; i++) delete [] Yin[i];
  delete [] Yin;
  DelCSL();
  DelGS();
  DelFG();
  printf("# CPScodec; deleted\n");
}

//================================================================================
void CPScodec::Encode(int **CCW, const int **IWB, const unsigned char *CWA){
  int pos=0, v;
  for(int i=0;i<Ns;i++){
    for(int j=0;j<Nu;j++){
      assert(CWA[pos]==0 || CWA[pos]==1);
      assert(IWB[i][j]>=0 && IWB[i][j]<QB2);
      v = CWA[pos]*(QB2) + IWB[i][j];
      assert(v>=0 && v<QB);
      memcpy( CCW[pos], CSL[v], sizeof(int)*4 );
      pos++;
    } // for j
  } // for i
  assert(pos==Nb);
}

//================================================================================
void CPScodec::Decode(double **Pout, const unsigned char **RW, int *Nb2, const int **dbgIW, const double **Pin){
  int jdx, idx;
  int Nb2max = max(Nb2,Nseq);
  int Nb2min = min(Nb2,Nseq);
  assert( Nb2min>=Nb+Dmin && Nb2max<=Nb+Dmax );
  InitFG(RW, Pin, Nb2);
  // ----- PXU
  for(idx=0; idx<Nb; idx++) CalcPXU(idx);
  // ----- forward
  for(idx=0; idx<Nb; idx++){
    for(jdx=0; jdx<Nseq; jdx++) CalcPDf(jdx,idx,Nb2[jdx]);
    CalcPU(idx);
    if(idx%Nu != Nu-1){
      for(jdx=0; jdx<Nseq; jdx++) CalcPF(jdx,idx,Nb2[jdx]);
    } // if idx
  } // for idx
  // ----- backward
  for(idx=Nb-1; idx>=0; idx--){
    for(jdx=0; jdx<Nseq; jdx++) CalcPDb(jdx,idx,Nb2[jdx]);
    CalcPU(idx);
    if(idx%Nu != 0){
      for(jdx=0; jdx<Nseq; jdx++) CalcPB(jdx,idx,Nb2[jdx]);
    } // if idx
  } // for idx
  // ----- downward
  for(idx=0; idx<Nb; idx++){
    for(jdx=0; jdx<Nseq; jdx++) CalcPD(jdx,idx,Nb2[jdx]);
    CalcPXD(idx);
    CalcPOC(idx);
  } // for idx
  // ----- set Pout
  CopyMat(Pout, (const double **)POC, Nb, QB2);
  //(dbg)
  //PrintFG((const int **)dbgIW);
}

//================================================================================
void CPScodec::Decode(double **Pout, const unsigned char **RW, int *Nb2, const int **dbgIW){
  double **Pin = new double * [Nb];
  for(int i=0;i<Nb;i++){
    Pin[i] = new double [QB2];
    for(int j=0;j<QB2;j++) Pin[i][j] = (double)1.0 / QB2; // unif
  } // for i
  //-----
  Decode(Pout,RW,Nb2,dbgIW,(const double **)Pin);
  //-----
  for(int i=0;i<Nb;i++) delete [] Pin[i];
  delete [] Pin;
}

//================================================================================
void CPScodec::InitFG(){
  ClearMat(PF,  Nseq, Nb+1, Drng);
  ClearMat(PB,  Nseq, Nb+1, Drng);
  ClearMat(PU,  Nseq, Nb,  QB);
  ClearMat(PD,  Nseq, Nb,  QB);
  ClearMat(PM,  Nb, QB);
  ClearMat(PXU, Nb, QB);
  ClearMat(PXD, Nb, QB);
  ClearMat(PIB, Nb, 2);
  ClearMat(PIC, Nb, QB2);
  ClearMat(POC, Nb, QB2);
}

//================================================================================
void CPScodec::InitFG(const unsigned char **RW, const double **Pin, int *Nb2){
  CopyMat(PIC, Pin, Nb, QB2);
  for(int j=0;j<Nseq;j++){
    memcpy(Yin[j],RW[j],sizeof(unsigned char)*Nb2[j]);
  } // for j
}

//================================================================================
void CPScodec::SetFGD(const double ***PrD){
  for(int jdx=0; jdx<Nseq; jdx++){
    for(int idx=0; idx<Ns+1; idx++){
      memcpy( PF[jdx][idx*Nu], PrD[jdx][idx], sizeof(double)*Drng );
      memcpy( PB[jdx][idx*Nu], PrD[jdx][idx], sizeof(double)*Drng );
    } // idx
  } // for jdx
}

//================================================================================
void CPScodec::SetFGB(const double **PrB2){
  for(int i=0; i<Nb; i++){
    memcpy( PIB[i], PrB2[i], sizeof(double)*2 );
  } // for i
}

//================================================================================
void CPScodec::PrintCSL(){
  printf("# CPScodec::CSL\n");
  printf("#  [0] ");
  for(int i=0; i<QB2;i++) printf("%02d,%02d,%02d,%02d  ",CSL[i][0],CSL[i][1],CSL[i][2],CSL[i][3]);
  printf("\n");
  printf("#  [1] ");
  for(int i=QB2;i<QB;i++) printf("%02d,%02d,%02d,%02d  ",CSL[i][0],CSL[i][1],CSL[i][2],CSL[i][3]);
  printf("\n");
}

//================================================================================
void CPScodec::PrintGS(){
  int y0,y1;
  printf("# CPScodec::GS1 p(y|s)\n");
  for(y0=0; y0<4; y0++){
    printf("#  [%d] ",y0);
    PrintVect( GS1[y0], QB, "", "\n");
  } // for y0
  printf("# CPScodec::GS2 p(y0,y1|s)\n");
  for(y0=0; y0<4; y0++){
    for(y1=0; y1<4; y1++){
      printf("#  [%d,%d] ",y0,y1);
      PrintVect( GS2[y0][y1], QB, "", "\n");
    } // for y1
  } // for y0
}

//================================================================================
void CPScodec::PrintFG(int idx, int dbgIWi){
  assert(idx>=0 && idx<Nb);
  printf("[%04d] IW=%d\n", idx, dbgIWi);
  PrintVect( PIC[idx], QB2,  " PIC: ","\n");
  PrintVect( POC[idx], QB2,  " POC: ","\n");
  PrintVect( PIB[idx], 2,    " PIB: ","\n");
  PrintVect( PXU[idx], QB,   " PXU: ","\n");
  PrintVect( PXD[idx], QB,   " PXD: ","\n");
  PrintVect( PM[ idx], QB,   " PM:  ","\n");
  for(int jdx=0; jdx<Nseq; jdx++){
    printf(" (%d) Y=%u\n",jdx,Yin[jdx][idx]);
    PrintVect( PU[jdx][idx  ], QB,   " PU:  ", "\n" );
    PrintVect( PD[jdx][idx  ], QB,   " PD:  ", "\n" );
    PrintVect( PF[jdx][idx  ], Drng, " PF:  ", "\n" );
    PrintVect( PB[jdx][idx+1], Drng, " PB:  ", "\n" );
  } // for jdx
}

//================================================================================
void CPScodec::PrintFG(int idx){
  PrintFG(idx,0);
}

//================================================================================
void CPScodec::PrintFG(const int **dbgIW){
  for(int idx=0; idx<Nb; idx++) PrintFG(idx, dbgIW[idx/Nu][idx%Nu]);
}

//================================================================================
void CPScodec::PrintFG(){
  for(int idx=0; idx<Nb; idx++) PrintFG(idx);
}

//================================================================================
int    CPScodec::Get_Ns(){  return Ns;}
int    CPScodec::Get_Nu(){  return Nu;}
int    CPScodec::Get_Nb(){  return Nb;}
int    CPScodec::Get_Nseq(){return Nseq;}
int    CPScodec::Get_Kc(){  return Kc;}
int    CPScodec::Get_Dmin(){return Dmin;}
int    CPScodec::Get_Dmax(){return Dmax;}
int    CPScodec::Get_Drng(){return Drng;}
int    CPScodec::Get_QB(){  return QB;}
double CPScodec::Get_Pi(){  return Pi;}
double CPScodec::Get_Pd(){  return Pd;}
double CPScodec::Get_Ps(){  return Ps;}
double CPScodec::Get_Pt(){  return Pt;}
