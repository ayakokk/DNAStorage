#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "NBIDSchannel.hpp"
#include "bmatrix.hpp"
#include "func.hpp"
#include "InnerCodebook.hpp"
#include "Gtable.hpp"
#include "FBalg4.hpp"

#define WCmax 100
#define RLCmax   15  // run length conter
#define WinSize  10  // GC-balance window
#define Tout     10  // output interval 
#define HardDecThNum 10

double HardDecThList[HardDecThNum] = {0.95,0.90,0.85,0.80,0.75,0.70,0.65,0.60,0.55,0.50}; 

void HardDecErs(int *HDout, const double **Px, int Nseg, int beta4p, double Pers);
void dbg_print(const int *ICin, const unsigned char *ICout,
	       const unsigned char *IDin, const int *IDout, int Nseg, int beta);

//=================================================================================
int main(int argc, char *argv[]){
  double Pid, Ps;
  char   *fnI, *fnG;  // inner CB & Gtable
  int    Nseg, seed;
  int    beta, beta4p, N, Dmax, N2;
  
  if(argc!=7){
    fprintf(stderr,"Usage: %s <Pid> <Ps> <Nseg> <InnerCB.bmat> <Gtable.bin> <seed|-1>\n",argv[0]);
    return 1;
  }
  Pid = atof(argv[1]);
  Ps  = atof(argv[2]);
  Nseg= atoi(argv[3]);
  fnI =      argv[4];
  fnG =      argv[5];
  seed= atoi(argv[6]);
  if(seed==-1) seed=(int)time(NULL);
  srandom(seed);
  printf("# (Pid,Ps)=(%e,%e) Nseg=%d WinSize=%d [%d]\n",Pid,Ps,Nseg,WinSize,seed);
  printf("# InnerCB: %s\n",fnI);
  printf("# Gtable:  %s\n",fnG);
  assert(Pid>=0.0 && Pid<1.0);
  assert(Ps >=0.0 && Ps <1.0);
  PrintVect(HardDecThList,HardDecThNum,"# HardDecThList: ","\n");

  //-----
  class InnerCodebook *ICB = new class InnerCodebook(fnI);
  beta  = ICB->Get_beta();
  beta4p= ICB->Get_beta4p();
  N     = Nseg*beta;
  Dmax  = (int)ceil((double)N*Pid+5); //???
  //Dmax  = 20;
  printf("# beta=%d(%d) N=%d Dmax=%d\n",beta,beta4p,N,Dmax);
  class NBIDSchannel *CH = new class NBIDSchannel(N,4,Dmax,Pid,Ps);
  class FBalg4      *FBA = new class FBalg4(Nseg,Dmax,fnG);
  assert(FBA->getBeta()==beta && FBA->getPid()==Pid && FBA->getPs()==Ps);
  
  //-----
  int           *ICin  = new int [Nseg];
  unsigned char *ICout = new unsigned char [Nseg*beta];
  int           *ICoutI= new int [Nseg];  // for eval
  unsigned char *IDin  = new unsigned char [Nseg*beta+Dmax];
  //int           *IDout = new int [Nseg];
  int           *HDout = new int [Nseg];  // hard decision of Pout with erasure
  double **Px  = new double * [Nseg];
  double **Pout= new double * [Nseg];
  for(int i=0;i<Nseg;i++){
    Px[i]  = new double [beta4p];
    Pout[i]= new double [beta4p];
  } // for i
  ICB->GetPrior(Px,Nseg);
  //PrintArray((const double **)Px,Nseg,beta4p);

  //-----
  int    wc, t0=(int)time(NULL);
  double ee, ees =0.0;  // error entropy
  double arl,arls=0.0;  // ave run-length
  double gcb,gcbs=0.0;  // windowed GC balance (ave of abs)
  unsigned long *RLC = new unsigned long [RLCmax+1];
  unsigned long *BLC = new unsigned long [WinSize*2+1];
  ClearVect(RLC,RLCmax+1);
  ClearVect(BLC,WinSize*2+1);

  unsigned long *** CntXY = new unsigned long ** [HardDecThNum];
  for(int i=0;i<HardDecThNum;i++){
    CntXY[i] = new unsigned long * [beta4p];
    for(int j=0;j<beta4p;j++){
      CntXY[i][j] = new unsigned long [beta4p+1]; // +erasure
      ClearVect(CntXY[i][j],beta4p+1);
    } // for j
  } // for i
  
  //-----
  for(wc=1;wc<=WCmax;wc++){
    ICB->GenRndInfo(ICin,Nseg);
    ICB->Encode(ICout,ICin,Nseg);
    N2 = CH->transmit(IDin,ICout);
    FBA->dbgSetX(ICout); //(dbg)
    FBA->calc(Pout,(const double **)Px,IDin,N2);

    // eval
    for(int i=0;i<Nseg;i++) ICoutI[i] = ConvVectInt(&ICout[i*beta],beta);
    for(int i=0;i<HardDecThNum;i++){
      HardDecErs(HDout,(const double**)Pout,Nseg,beta4p,HardDecThList[i]);
    } // for i
    
    ee  = ErrEntropy(ICoutI,(const double **)Pout,Nseg,beta4p);
    arl = AveRunLength(ICout,Nseg*beta);
    gcb = GCbalanceWin(ICout,Nseg*beta,WinSize);
    ees +=ee;
    arls+=arl;
    gcbs+=gcb;
    CntRunLength(ICout,Nseg*beta,RLC,RLCmax);
    GCbalanceCnt(ICout,Nseg*beta,WinSize,BLC);
    
    //TMP
    //memcpy(IDin,ICout,sizeof(unsigned char)*Nseg*beta); //TMP: noiseless
    //ICB->Decode(IDout,IDin,Nseg);

    if(wc==WCmax || (int)time(NULL)-t0>Tout){
      printf("%04d(%d) (%.2e %.2e) XE:%.2e (%.2e %.2e | %.2e %.2e)\n",
	     wc,N2,ee,ees/wc,ees/(wc*beta*2),arl,(double)arls/wc,gcb,(double)gcbs/wc);
      t0=(int)time(NULL);
    }
    //printf(" RLC: "); PrintVect(RLC,RLCmax+1); printf("\n");
    //dbg_print(ICin,ICout,IDin,IDout,Nseg,beta);
    //printf("%04d %d %d\n",wc,HammingDist(ICin,IDout,Nseg),N2);
    //for(int i=0;i<Nseg;i++) printf(" %04d: Pxout[%02X]=%.2e\n",i,ICoutI[i],Pout[i][ICoutI[i]]);
  } // for wc
  //---
  PrintVect(RLC,RLCmax+1,   "RL cnt: ","\n");
  PrintVect(BLC,WinSize*2+1,"BL cnt: ","\n"); 
  
  delete ICB;
  delete CH;
  delete FBA;
  delete [] ICin;
  delete [] ICout;
  delete [] ICoutI;
  delete [] IDin;
  //delete [] IDout;
  delete [] HDout;
  for(int i=0;i<Nseg;i++){
    delete [] Px[i];
    delete [] Pout[i];
  }
  delete [] Px;
  delete [] Pout;
  delete [] RLC;
  for(int i=0;i<HardDecThNum;i++){
    for(int j=0;j<beta4p;j++) delete [] CntXY[i][j];
    delete [] CntXY[i];
  } // for i
  delete [] CntXY;
  return 0;
}

//=================================================================================
void HardDecErs(int *HDout, const double **Px, int Nseg, int beta4p, double Pers){
  assert(Pers>0 && Pers<1);
  int idx,v;
  for(idx=0;idx<Nseg;idx++){
    v = ArgMax(Px[idx],beta4p);
    if(Px[idx][v]<Pers) v=beta4p; // erasure
    HDout[idx]=v;
  } // for idx
}

//=================================================================================
void dbg_print(const int *ICin, const unsigned char *ICout,
	       const unsigned char *IDin, const int *IDout, int Nseg, int beta){
  for(int i=0;i<Nseg;i++){
    printf("%04d: %04X->",i,ICin[i]);
    PrintVect4(&ICout[i*beta],beta);
    printf("->");
    PrintVect4(&IDin[i*beta],beta);
    printf("->%04X",IDout[i]);
    printf("\n");
  } // for i
}
