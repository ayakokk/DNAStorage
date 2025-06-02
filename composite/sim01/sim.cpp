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
#include "CSFBAdec.hpp"

#define BSIZE 8192
#define Tint  10 //sec
#define OutListSize 3
#define WCmax 100000
//#define WCmax 1000

void OutputConv(long *DWL, const double **P, int N, int Q);
long PdistArgMaxLong(const double *P, int Q, int LS);
void dbgPrint(const int *IW, const int *DW, 
	      const unsigned char *CW, const unsigned char **RW, const double **Pout, 
	      const int *dbgDR, int Nb, int *Nb2, int Nu, int Q, int Nseq);

//================================================================================
int main(int argc, char *argv[]){
  int Rho;         // run-length    [constraint.txt]
  int ell,Delta;   // local-balance [constraint.txt]
  int N;           // block length (symbols)
  int Nb;          // block length (bits)
  int Q, Nu;       // numCW, symbol-len [ICB]
  int Nseq;        // number of sequences (=r)
  int Ndi;         // number of decoding iterations
  int seed;
  int t0=(int)time(NULL);
  double Pi,Pd,Ps; // IDS prob
  char *fn;
  char *fncb    = new char [BSIZE];
  char *fnconst = new char [BSIZE];
  char *fncm    = new char [BSIZE];
  if(argc!=9){
    fprintf(stderr,"Usage: %s <ICB_dir> <N> <Pi> <Pd> <Ps> <Nseq> <NumDecItr> <seed|-1>\n",argv[0]);
    return 1;
  } // if
  fn   =      argv[1];
  N    = atoi(argv[2]);
  Pi   = atof(argv[3]);
  Pd   = atof(argv[4]);
  Ps   = atof(argv[5]);
  Nseq = atoi(argv[6]);
  Ndi  = atoi(argv[7]);
  seed = atoi(argv[8]);
  if(seed==-1) seed = (int)time(NULL);
  srandom(seed);
  assert(N>0 && Nseq>0);
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
  printf("# Q=%d N=%d Nu=%d Nb=%d (Pi,Pd,Ps)=(%e,%e,%e) Nseq=%d Ndi=%d [%d]\n",Q,N,Nu,Nb,Pi,Pd,Ps,Nseq,Ndi,seed);
  printf("# ICB:   %s\n",fncb);
  printf("# Const: %s\n",fnconst);
  printf("# EncCM: %s\n",fncm);
  class ChannelMatrix *ECM = new class ChannelMatrix(fncm);
  class IDSchannel    *CH  = new class IDSchannel(Nb,Pi,Pd,Ps);
  class CSFBAdec      *DEC = new class CSFBAdec(ICB,ECM,CH,Nseq,Ndi);
  class ChannelMatrix *DCM = new class ChannelMatrix(Q,(int)pow(Q,OutListSize));
  int         *dbgDR = new int [Nb+1];                         // (dbg)drift vector
  int            *IW = new int [N];                            // information word
  unsigned char  *CW = new unsigned char [Nb];                 // codeword
  int            *DW = new int [N];                            // decoded word
  int           *Nb2 = new int [Nseq];                         // received word length (bits)
  long          *DWL = new long [N];                           // (Pout->list->long)
  double      **Pout = new double * [N];                       // output PMF [N][Q]
  unsigned char **RW = new unsigned char * [Nseq];             // received word [Nseq][N+Dmax]
  for(int i=0;i<N;   i++) Pout[i] = new double [Q];          
  for(int i=0;i<Nseq;i++) RW[i]   = new unsigned char [Nb + CH->GetDmax()];
  int wc;
  long ec,ecmax=0,es=0;
  
  //-----
  for(wc=1;wc<=WCmax;wc++){
    RandVect(IW,N,0,Q-1);
    ICB->Encode(CW,IW,N);
    for(int i=0;i<Nseq;i++) Nb2[i] = CH->transmit(RW[i],CW);
    //Nb2 = CH->transmit(RW,CW);
    DEC->Decode(Pout,(const unsigned char **)RW,Nb2,IW);
    //DEC->Decode(Pout,RW,Nb2,IW);
    HardDecision(DW,(const double **)Pout,N,Q);
    OutputConv(DWL,(const double **)Pout,N,Q);
    for(int i=0;i<N;i++) DCM->countup(IW[i],DWL[i]);
    ec = HammingDist(IW,DW,N);
    es += ec;
    ecmax = max(ec,ecmax);

    if( time(NULL)-t0 > Tint || wc==WCmax){
      printf("%04d %ld/%ld %ld %e : %e %e %e\n",
	     wc,es,(long)wc*N,ecmax,(double)es/(wc*N), DCM->Hx(), DCM->Hxy(), DCM->Ixy());
      t0 = (int)time(NULL);
    } // if wc 

    //(dbg)
    //CH->GetDR(dbgDR);
    //dbgPrint(IW,DW,CW,(const unsigned char **)RW,(const double**)Pout,dbgDR,Nb,Nb2,Nu,Q,Nseq);
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
  delete [] DW;
  delete [] DWL;
  delete [] Nb2;
  delete [] fncb;
  delete [] fnconst;
  delete [] fncm;
  for(int i=0;i<N;   i++) delete [] Pout[i];
  for(int i=0;i<Nseq;i++) delete [] RW[i];
  delete [] Pout;
  delete [] RW;
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
	      const unsigned char *CW, const unsigned char **RW, const double **Pout,
	      const int *dbgDR, int Nb, int *Nb2, int Nu, int Q, int Nseq){
  int idx;
  int Nb2max = max(Nb2,Nseq);
  for(int i=0;i<max(Nb,Nb2max);i++){
    if(i%Nu==0){
      if(i<Nb){
	idx = i/Nu;
	printf("[%03d] %02d %02d\n",idx,IW[idx],DW[idx]);
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
    for(int j=0;j<Nseq;j++){
      if(i<Nb2[j]) printf("%u",RW[j][i]);
      else         printf("-");
    } // for j
    printf(" ");
    //---
    if(i<Nb+1) printf("(%+03d) ",dbgDR[i]);
    else       printf("(---) ");
    //---
    printf("\n");
  } // for i
}
