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

#define WCmax 10
#define Tout   0  // output interval 

void dbg_print(const int *ICin, const unsigned char *ICout,
	       const unsigned char *IDin, const int *IDout, int Nseg, int beta);

//=================================================================================
int main(int argc, char *argv[]){
  double Pid, Ps;
  char   *fnI, *fnG;  // inner CB & Gtable
  int    Nseg, seed;
  int    beta, beta4p, nu, N, Dmax, N2;
  
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
  printf("# (Pid,Ps)=(%e,%e) Nseg=%d [%d]\n",Pid,Ps,Nseg,seed);
  printf("# InnerCB: %s\n",fnI);
  printf("# Gtable:  %s\n",fnG);
  assert(Pid>=0.0 && Pid<1.0);
  assert(Ps >=0.0 && Ps <1.0);

  //-----
  class InnerCodebook *ICB = new class InnerCodebook(fnI);
  beta  = ICB->Get_beta();
  beta4p= ICB->Get_beta4p();
  nu    = ICB->Get_nu();
  N     = Nseg*beta;
  Dmax  = (int)ceil((double)N*Pid+5); //???
  //Dmax  = 20;
  printf("# beta=%d(%d) nu=%d N=%d Dmax=%d\n",beta,beta4p,nu,N,Dmax);
  int *B = new int [nu];
  ICB->Get_B(B);
  PrintVect(B,nu,"# B=","\n");
  class NBIDSchannel *CH = new class NBIDSchannel(N,4,Dmax,Pid,Ps);
  class FBalg4      *FBA = new class FBalg4(Nseg,Dmax,fnG);
  assert(FBA->getBeta()==beta && FBA->getPid()==Pid && FBA->getPs()==Ps);
  
  //-----
  int           *ICin  = new int [Nseg];
  unsigned char *ICout = new unsigned char [Nseg*beta];
  int           *ICoutI= new int [Nseg];  // for eval
  unsigned char *IDin  = new unsigned char [Nseg*beta+Dmax];
  int           *IDout = new int [Nseg];
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
  double Ixy,sIxy=0.0;
  
  //-----
  for(wc=1;wc<=WCmax;wc++){
    ICB->GenRndInfo(ICin,Nseg);
    ICB->Encode(ICout,ICin,Nseg);
    N2 = CH->transmit(IDin,ICout);
    FBA->dbgSetX(ICout); //(dbg)
    Ixy = FBA->calcIxy((const double **)Px,IDin,N2,(const unsigned char *)ICout,nu,B);
    sIxy += Ixy;
    
    // eval
    //for(int i=0;i<Nseg;i++) ICoutI[i] = ConvVectInt(&ICout[i*beta],beta);
    
    //TMP
    //memcpy(IDin,ICout,sizeof(unsigned char)*Nseg*beta); //TMP: noiseless
    //ICB->Decode(IDout,IDin,Nseg);

    if(wc==WCmax || (int)time(NULL)-t0>=Tout){
      printf("%04d(%d): (%e %e) %e\n",wc,N2,Ixy,sIxy/wc,sIxy/(wc*beta*2.0));
      t0=(int)time(NULL);
    }
    //dbg_print(ICin,ICout,IDin,IDout,Nseg,beta);
    //printf("%04d %d %d\n",wc,HammingDist(ICin,IDout,Nseg),N2);
    //for(int i=0;i<Nseg;i++) printf(" %04d: Pxout[%02X]=%.2e\n",i,ICoutI[i],Pout[i][ICoutI[i]]);
  } // for wc
  //---
  
  delete ICB;
  delete CH;
  delete FBA;
  delete [] ICin;
  delete [] ICout;
  delete [] ICoutI;
  delete [] IDin;
  delete [] IDout;
  for(int i=0;i<Nseg;i++){
    delete [] Px[i];
    delete [] Pout[i];
  }
  delete [] Px;
  delete [] Pout;
  delete [] B;
  return 0;
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
