#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"
#include "ChannelMatrix.hpp"
#include "InnerCodebook.hpp"
#include "IDSchannel.hpp"
#include "SLFBAdec.hpp"

#define BSIZE 8192
#define OutListSize 3
#define WCmax 100000

void OutputConv(long *DWL, const double **P, int N, int Q);
long PdistArgMaxLong(const double *P, int Q, int LS);
void dbgPrint(const int *IW, const int *DW, 
	      const unsigned char *CW, const unsigned char *RW, const double **Pout, 
	      const int *dbgDR, int Nb, int Nb2, int Nu, int Q);

//================================================================================
int main(int argc, char *argv[]){
  int Rho;         // run-length    [constraint.txt]
  int ell,Delta;   // local-balance [constraint.txt]
  int N;           // block length (symbols)
  int Nb,Nb2;      // block length & recv length (bits)
  int Q, Nu;       // numCW, symbol-len [ICB]
  int seed;
  double Pi,Pd,Ps; // IDS prob
  char *fn;
  char *fncb    = new char [BSIZE];
  char *fnconst = new char [BSIZE];
  char *fncm    = new char [BSIZE];
  if(argc!=7){
    fprintf(stderr,"Usage: %s <ICB_dir> <N> <Pi> <Pd> <Ps> <seed|-1>\n",argv[0]);
    return 1;
  } // if
  fn   =      argv[1];
  N    = atoi(argv[2]);
  Pi   = atof(argv[3]);
  Pd   = atof(argv[4]);
  Ps   = atof(argv[5]);
  seed = atoi(argv[6]);
  if(seed==-1) seed = (int)time(NULL);
  srandom(seed);
  assert(N>0);
  assert(Pi>=0.0 && Pi<0.5);
  assert(Pd>=0.0 && Pd<0.5);
  assert(Ps>=0.0 && Ps<0.5);
  snprintf(fncb,   BSIZE,"%s/cb.txt",        fn);  // inner codebook (in)
  snprintf(fnconst,BSIZE,"%s/constraint.txt",fn);  // constraints (in)
  snprintf(fncm,   BSIZE,"%s/EncCM.bin",     fn);  // encoding channel matrix (in)
  ReadConstraints(fnconst, &Rho, &ell, &Delta);
  class InnerCodebook *ICB = new class InnerCodebook(fncb,Rho,ell,Delta);
  Q    = ICB->Get_numCW();
  Nu   = ICB->Get_Nu();
  Nb   = N*Nu;
  printf("# Q=%d N=%d Nu=%d Nb=%d (Pi,Pd,Ps)=(%e,%e,%e) [%d]\n",Q,N,Nu,Nb,Pi,Pd,Ps,seed);
  printf("# ICB:   %s\n",fncb);
  printf("# Const: %s\n",fnconst);
  printf("# EncCM: %s\n",fncm);
  class ChannelMatrix *ECM = new class ChannelMatrix(fncm);
  class IDSchannel    *CH  = new class IDSchannel(Nb,Pi,Pd,Ps);
  class SLFBAdec      *DEC = new class SLFBAdec(ICB,ECM,CH);
  class ChannelMatrix *DCM = new class ChannelMatrix(Q,(int)pow(Q,OutListSize));
  int        *dbgDR = new int [Nb+1];                         // (dbg)drift vector
  int           *IW = new int [N];                            // information word
  unsigned char *CW = new unsigned char [Nb];                 // codeword
  unsigned char *RW = new unsigned char [Nb + CH->GetDmax()]; // received word
  int           *DW = new int [N];                            // decoded word
  long         *DWL = new long [N];                           // (Pout->list->long)
  double **Pout = new double * [N];
  for(int i=0;i<N;i++) Pout[i] = new double [Q];
  int wc;
  long ec,ecmax=0,es=0;
  
  //-----
  for(wc=1;wc<=WCmax;wc++){
    RandVect(IW,N,0,Q-1);
    ICB->Encode(CW,IW,N);
    Nb2 = CH->transmit(RW,CW);
    DEC->Decode(Pout,RW,Nb2,IW);
    HardDecision(DW,(const double **)Pout,N,Q);
    OutputConv(DWL,(const double **)Pout,N,Q);
    for(int i=0;i<N;i++) DCM->countup(IW[i],DWL[i]);
    ec = HammingDist(IW,DW,N);
    es += ec;
    ecmax = max(ec,ecmax);

    if(wc%1000==0 || wc==WCmax){
      printf("%04d %ld/%ld %ld %e : %e %e %e\n",
	     wc,es,(long)wc*N,ecmax,(double)es/(wc*N), DCM->Hx(), DCM->Hxy(), DCM->Ixy());
    } // if wc 

    //(dbg)
    //CH->GetDR(dbgDR);
    //dbgPrint(IW,DW,CW,RW,(const double**)Pout,dbgDR,Nb,Nb2,Nu,Q);
    //printf("%04d %ld %ld/%ld %e\n",wc,ec,es,(long)wc*N,(double)es/(wc*N));
  } // for wc
  //-----
  //DCM->PrintCnt();

  delete ICB;
  delete ECM;
  delete CH;
  delete DEC;
  delete DCM;
  delete [] dbgDR;
  delete [] IW;
  delete [] CW;
  delete [] RW;
  delete [] DW;
  delete [] DWL;
  delete [] fncb;
  delete [] fnconst;
  delete [] fncm;
  for(int i=0;i<N;i++) delete [] Pout[i];
  delete [] Pout;
  return 0;
}

//================================================================================
void OutputConv(long *DWL, const double **P, int N, int Q){
  for(int i=0;i<N;i++){
    DWL[i] = PdistArgMaxLong(P[i],Q,OutListSize);
    //printf("%03d: %ld\n",i,DWL[i]);
    //PrintVect(P[i],Q," ","\n");
  } // for i
}

//================================================================================
long PdistArgMaxLong(const double *P, int Q, int LS){
  assert(Q>0 && LS>0 && LS<=Q);
  long v, val=0;
  double *PX = new double [Q];
  memcpy(PX,P,sizeof(double)*Q);
  for(int i=0; i<Q; i++) assert( P[i]>=0.0 && P[i]<=1.0 );
  for(int i=0; i<LS; i++){
    v = argmax(PX,Q);
    val = val*Q + v;
    PX[v] = -1.0;
  } // for i
  delete [] PX;
  return val;
}

//================================================================================
void dbgPrint(const int *IW, const int *DW, 
	      const unsigned char *CW, const unsigned char *RW, const double **Pout,
	      const int *dbgDR, int Nb, int Nb2, int Nu, int Q){
  int idx;
  for(int i=0;i<max(Nb,Nb2);i++){
    if(i%Nu==0){
      if(i<Nb){
	idx = i/Nu;
	printf("[%03d] %002d %002d\n",idx,IW[idx],DW[idx]);
	PrintVect(Pout[idx],Q,"","\n");
      } else {
	printf("[---]\n");
      } // if i<Nb
    } // if i%Nu
    //---
    printf("%04d: ",i);
    //---
    if(i<Nb) printf("%u ",CW[i]);
    else     printf("- ");
    //---
    if(i<Nb2) printf("%u ",RW[i]);
    else      printf("- ");
    //---
    if(i<Nb+1) printf("(%+03d) ",dbgDR[i]);
    else       printf("(---) ");
    //---
    printf("\n");
  } // for i
}
