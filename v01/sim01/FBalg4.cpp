#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "bmatrix.hpp"
#include "func.hpp"
#include "Gtable.hpp"
#include "FBalg4.hpp"

//=================================================================================
void FBalg4::ArrayCopy(double **Qt, const double **Qs, int len0, int len1){
  for(int i=0;i<len0;i++){
    memcpy(Qt[i],Qs[i],sizeof(double)*len1);
  }
}

//=================================================================================
void FBalg4::init(const double **Pin, const unsigned char *Yin, int N2){
  ArrayCopy(Pxin,Pin,Nseg,beta4p);
  memcpy(Y,Yin,sizeof(unsigned char)*(N+Dmax));
  for(int i=0;i<=Nseg;i++){
    for(int j=0;j<DmaxW;j++){
      PF[i][j]=0.0;
      PB[i][j]=0.0;
    } // for i
  } // for j
  PF[0   ][0   +Dmax]=1.0;
  PB[Nseg][N2-N+Dmax]=1.0;
  for(int i=0;i<Nseg;i++)
    for(int j=0;j<beta4p;j++) Pxout[i][j]=0.0;
}

//=================================================================================
void FBalg4::forward(int idx, int N2){
  int    d0,  d1;
  int    pos0,pos1;
  int    len,y;
  double psum;
  assert(idx>=0 && idx<Nseg);
  for(d1=-Dmax;d1<=Dmax;d1++){
    PF[idx+1][d1+Dmax]=0.0; // init
    pos1 = (idx*beta)+beta-1+d1;
    if(pos1>=-1 && pos1<N2){
      for(d0=-Dmax;d0<=Dmax;d0++){
	pos0 = (idx*beta)+d0;
	len  = pos1-pos0+1;
	if(pos0>=0 && pos0<N2 && len>=0 && len<=beta*2){
	  psum=0;
	  y = ConvVectInt(&Y[pos0],len);
	  for(auto xp=GTB->LIST[len][y].begin(); xp!=GTB->LIST[len][y].end();++xp){
	    psum += (Pxin[idx][(*xp).x] * (*xp).p);
	    //printf("%e idx=%d d1=%d pos1=%d d0=%d pos0=%d len=%d y=%d x=%d\n",(*xp).p,idx,d1,pos1,d0,pos0,len,y,(*xp).x);
	  } // for xp
	  psum *= PF[idx][d0+Dmax];
	  PF[idx+1][d1+Dmax] += psum;
	} // if pos0 && len
      } // for d0
    } // if pos1
  } // for d1
  Normalize(PF[idx+1],DmaxW);
}

//=================================================================================
void FBalg4::backward(int idx, int N2){
  int    d0,  d1;
  int    pos0,pos1;
  int    len,y;
  double psum;
  assert(idx>=0 && idx<Nseg);
  for(d0=-Dmax;d0<=Dmax;d0++){
    PB[idx][d0+Dmax]=0.0; // init
    pos0 = (idx*beta)+d0;
    if(pos0>=0 && pos0<N2){
      for(d1=-Dmax;d1<=Dmax;d1++){
	pos1 = (idx*beta)+beta-1+d1;
	len  = pos1-pos0+1;
	if(pos1>=-1 && pos0<N2 && len>=0 && len<=beta*2){
	  psum=0;
	  y = ConvVectInt(&Y[pos0],len);
	  for(auto xp=GTB->LIST[len][y].begin(); xp!=GTB->LIST[len][y].end();++xp){
	    psum += (Pxin[idx][(*xp).x] * (*xp).p);
	    //printf("%e idx=%d d0=%d pos0=%d d1=%d pos1=%d len=%d y=%d x=%d\n",(*xp).p,idx,d0,pos0,d1,pos1,len,y,(*xp).x);
	  } // for xp
	  psum *= PB[idx+1][d1+Dmax];
	  PB[idx][d0+Dmax] += psum;
	} // if pos1 && len
      } // for d1
    } // if pos0
  } // for d0
  Normalize(PB[idx],DmaxW);
}

//=================================================================================
void FBalg4::downward(int idx, int N2){
  int    d0,  d1;
  int    pos0,pos1;
  int    len,y;
  double pp;
  assert(idx>=0 && idx<Nseg);
  for(int x=0;x<beta4p;x++) Pxout[idx][x]=0.0;
  
  for(d1=-Dmax;d1<=Dmax;d1++){
    pos1 = (idx*beta)+beta-1+d1;
    if(pos1>=-1 && pos1<N2){
      for(d0=-Dmax;d0<=Dmax;d0++){
	pos0 = (idx*beta)+d0;
	len  = pos1-pos0+1;
	if(pos0>=0 && pos0<N2 && len>=0 && len<=beta*2){
	  pp= PF[idx][d0+Dmax]*PB[idx+1][d1+Dmax];
	  y = ConvVectInt(&Y[pos0],len);
	  for(auto xp=GTB->LIST[len][y].begin(); xp!=GTB->LIST[len][y].end();++xp){
	    Pxout[idx][(*xp).x] += (pp * (*xp).p);
	    //printf("%e idx=%d d1=%d pos1=%d d0=%d pos0=%d len=%d y=%d x=%d\n",(*xp).p,idx,d1,pos1,d0,pos0,len,y,(*xp).x);
	  } // for xp
	} // if pos0 && len
      } // for d0
    } // if pos1
  } // for d1
  Normalize(Pxout[idx],beta4p);
}

//=================================================================================
double FBalg4::CalcLmdY(int idx, int N2){
  int    d0,  d1;
  int    pos0,pos1;
  int    len,y;
  double psum,lmd;
  assert(idx>=0 && idx<Nseg);
  for(d1=-Dmax;d1<=Dmax;d1++){
    PF[idx+1][d1+Dmax]=0.0; // init
    pos1 = (idx*beta)+beta-1+d1;
    if(pos1>=-1 && pos1<N2){
      for(d0=-Dmax;d0<=Dmax;d0++){
	pos0 = (idx*beta)+d0;
	len  = pos1-pos0+1;
	if(pos0>=0 && pos0<N2 && len>=0 && len<=beta*2){
	  psum=0;
	  y = ConvVectInt(&Y[pos0],len);
	  for(auto xp=GTB->LIST[len][y].begin(); xp!=GTB->LIST[len][y].end();++xp){
	    psum += (Pxin[idx][(*xp).x] * (*xp).p);
	    //printf("%e idx=%d d1=%d pos1=%d d0=%d pos0=%d len=%d y=%d x=%d\n",(*xp).p,idx,d1,pos1,d0,pos0,len,y,(*xp).x);
	  } // for xp
	  psum *= PF[idx][d0+Dmax];
	  PF[idx+1][d1+Dmax] += psum;
	} // if pos0 && len
      } // for d0
    } // if pos1
  } // for d1
  lmd = Normalize(PF[idx+1],DmaxW);
  return 1.0/lmd;
}

//=================================================================================
double FBalg4::CalcLmdXY(int idx, int N2, const int *X){
  int    d0,  d1;
  int    pos0,pos1;
  int    len,y;
  double psum,lmd;
  assert(idx>=0 && idx<Nseg);
  for(d1=-Dmax;d1<=Dmax;d1++){
    PF[idx+1][d1+Dmax]=0.0; // init
    pos1 = (idx*beta)+beta-1+d1;
    if(pos1>=-1 && pos1<N2){
      for(d0=-Dmax;d0<=Dmax;d0++){
	pos0 = (idx*beta)+d0;
	len  = pos1-pos0+1;
	if(pos0>=0 && pos0<N2 && len>=0 && len<=beta*2){
	  psum=0;
	  y = ConvVectInt(&Y[pos0],len);
	  for(auto xp=GTB->LIST[len][y].begin(); xp!=GTB->LIST[len][y].end();++xp){
	    if((*xp).x==X[idx]){
	      psum = Pxin[idx][(*xp).x] * (*xp).p;
	      break;
	    } // if 
	    //psum += (Pxin[idx][(*xp).x] * (*xp).p);
	    //printf("%e idx=%d d1=%d pos1=%d d0=%d pos0=%d len=%d y=%d x=%d\n",(*xp).p,idx,d1,pos1,d0,pos0,len,y,(*xp).x);
	  } // for xp
	  if(psum==0.0) psum = Pth;
	  psum *= PF[idx][d0+Dmax];
	  PF[idx+1][d1+Dmax] += psum;
	} // if pos0 && len
      } // for d0
    } // if pos1
  } // for d1
  lmd = Normalize(PF[idx+1],DmaxW);
  return 1.0/lmd;
}


//=================================================================================
//=================================================================================
//=================================================================================

//=================================================================================
FBalg4::FBalg4(int _Nseg, int _Dmax, const char *fnGtable){
  Nseg = _Nseg;
  Dmax = _Dmax;
  DmaxW= Dmax*2+1;
  printf("# FBalg4: Nseg=%d Dmax=%d(w=%d)\n",Nseg,Dmax,DmaxW);
  printf("# FBalg4: Gtable: %s\n",fnGtable);
  GTB =  new class Gtable(fnGtable);
  beta  = GTB->getBeta();
  Pid   = GTB->getPid();
  Ps    = GTB->getPs();
  Pth   = GTB->getPth();
  beta4p= (int)pow(4,beta);
  N     = Nseg*beta;
  printf("# FBalg4: beta=%d(%d) (Pid,Ps)=(%e,%e) Pth=%e N=%d\n",beta,beta4p,Pid,Ps,Pth,N);

  Pxin = new double * [Nseg];
  Pxout= new double * [Nseg];
  for(int i=0;i<Nseg;i++){
    Pxin[i] = new double [beta4p];
    Pxout[i]= new double [beta4p];
  } // for i
  PF = new double * [Nseg+1];
  PB = new double * [Nseg+1];
  for(int i=0;i<Nseg+1;i++){
    PF[i] = new double [DmaxW];
    PB[i] = new double [DmaxW];
  } // for i
  Y   = new unsigned char [N+Dmax];
  dbgX= new unsigned char [N];
}

//=================================================================================
FBalg4::~FBalg4(){
  delete GTB;
  for(int i=0;i<Nseg;i++){
    delete [] Pxin[i];
    delete [] Pxout[i];
  } // for i
  delete [] Pxin;
  delete [] Pxout;
  for(int i=0;i<Nseg+1;i++){
    delete [] PF[i];
    delete [] PB[i];
  } // for i
  delete [] PF;
  delete [] PB;
  delete [] Y;
  delete [] dbgX;
  printf("# FBalg4: deleted\n");
}

//=================================================================================
int    FBalg4::getBeta(){return beta;}
double FBalg4::getPid(){ return Pid;}
double FBalg4::getPs(){  return Ps;}
double FBalg4::getPth(){ return Pth;}

//=================================================================================
void FBalg4::calc(double **Pout, const double **Pin, const unsigned char *Yin, int N2){
  assert(N2>=N-Dmax && N2<=N+Dmax);
  init(Pin,Yin,N2);
  for(int idx=0;     idx<Nseg;idx++) forward( idx,N2);
  for(int idx=Nseg-1;idx>=0;  idx--) backward(idx,N2);
  for(int idx=0;     idx<Nseg;idx++) downward(idx,N2);
  ArrayCopy(Pout,(const double **)Pxout,Nseg,beta4p);
  // dumpFB();
  // dump(50);
}

//=================================================================================
double FBalg4::calcIxy(const double **Pin, const unsigned char *Yin, int N2, const unsigned char *Xin, int nu, const int *B){
  assert(N2>=N-Dmax && N2<=N+Dmax);
  double lmdY, sumY =0.0;
  double lmdXY,sumXY=0.0;
  double lpx=0;
  double Ixy;
  int *X = new int [Nseg];
  for(int i=0;i<Nseg;i++) X[i] = ConvVectInt(&Xin[i*beta],beta);

  init(Pin,Yin,N2);
  for(int idx=0;idx<Nseg;idx++){
    lmdY = CalcLmdY(idx,N2);
    assert(lmdY>0);
    sumY += log2(lmdY);
    //printf("Y: %d %e %e\n",idx,lmdY,sumY);
  } // for idx

  init(Pin,Yin,N2);
  for(int idx=0;idx<Nseg;idx++){
    lmdXY = CalcLmdXY(idx,N2,X);
    assert(lmdXY>0);
    sumXY += log2(lmdXY);
    //printf("XY: %d %e %e\n",idx,lmdXY,sumXY);
  } // for idx

  for(int idx=0;idx<Nseg;idx++) lpx += log2(1.0/pow(2,B[idx%nu]));
  //printf("X: %e %e\n",pxx,log2(pxx));

  Ixy = sumY - lpx - sumXY;
  printf("(Y:%e X:%e XY:%e I:%e)\n",sumY/Nseg,-lpx/Nseg,sumXY/Nseg,Ixy/Nseg);
  delete [] X;
  return Ixy/Nseg;
}

//=================================================================================
void FBalg4::dump(int idx){
  int w=16;
  assert(idx>=0 && idx<Nseg);
  printf("%04d: ",idx);
  PrintVect(&dbgX[idx*beta],beta);
  printf("(%02X) > ",ConvVectInt(&dbgX[idx*beta],beta));
  PrintVect(&Y[idx*beta],beta);
  printf("(%02X)\n",ConvVectInt(&Y[idx*beta],beta));
  //---
  printf(" PF: ");
  for(int i=-Dmax;i<=Dmax;i++){
    if(i==0) printf("*");
    printf("%.2e ",PF[idx][i+Dmax]);
  } // for i
  printf("\n");
  //---
  printf(" PB: ");
  for(int i=-Dmax;i<=Dmax;i++){
    if(i==0) printf("*");
    printf("%.2e ",PB[idx][i+Dmax]);
  } // for i
  printf("\n");
  //----
  printf(" Pxin:  ");
  for(int i=0;i<beta4p;i++){
    if(i>0 && i%w==0) printf("        ");
    printf("%02X:%.2e ",i,Pxin[idx][i]);
    if(i%w==w-1) printf("\n");
  } // for i
  //----
  printf(" Pxout: ");
  for(int i=0;i<beta4p;i++){
    if(i>0 && i%w==0) printf("        ");
    printf("%02X:%.2e ",i,Pxout[idx][i]);
    if(i%w==w-1) printf("\n");
  } // for i
}

//=================================================================================
void FBalg4::dump(){
  for(int idx=0;idx<Nseg;idx++) dump(idx);
}

//=================================================================================
void FBalg4::dumpFB(){
  for(int idx=0;idx<=Nseg;idx++){
    printf("%03d: ",idx);
    PrintVect(&dbgX[idx*beta],beta);
    printf(" > ");
    PrintVect(&Y[idx*beta],beta);
    printf("\n");
    // PF
    for(int d=-Dmax;d<=Dmax;d++){
      if(d==0) printf("|%.2e| ",PF[idx][d+Dmax]);
      else printf("%.2e ",      PF[idx][d+Dmax]);
    }
    printf("\n");
    // PB
    for(int d=-Dmax;d<=Dmax;d++){
      if(d==0) printf("|%.2e| ",PB[idx][d+Dmax]);
      else printf("%.2e ",      PB[idx][d+Dmax]);
    }
    printf("\n");
  } // for idx
}

//=================================================================================
void FBalg4::dbgSetX(const unsigned char *X){
  memcpy(dbgX,X,sizeof(unsigned char)*N);
}
