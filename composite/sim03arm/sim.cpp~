#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"
#include "ChannelMatrix.hpp"
#include "InnerCodebook.hpp"
#include "CIDSchannel.hpp"
#include "CSFBAdec.hpp"
#include "CPScodec.hpp"
#include "ConcatAB.hpp"

#define BSIZE       8192
#define Tint        10 //sec
#define OutListSize 3
#define Drange      2.0   // Dmin=-Drange*Nb*Pd Dmax=Drange*Nb*Pi  
#define WCmax 100000
//#define WCmax 1

void OutputConv(long  *DWL, const double **P, int Ns, int Q);
void OutputConv(long **DWL, const double **P, int Ns, int Nu, int Q);
long PdistArgMaxLong(const double *P, int Q, int LS);
void dbgPrint(const int *IWA, const int **IWB, const int *DWA, const int **DWB,  
	      const unsigned char *CWA, const int **CCW,
	      const unsigned char **RW, const double **PoutA, const double **PoutB, 
	      const unsigned char **dbgXV, const int **dbgDV, const double **dbgB2, 
	      int Nb, int *Nb2, int Nu, int QA, int QB2, int Nseq);

//================================================================================
int main(int argc, char *argv[]){
  int Rho;         // run-length    [constraint.txt]
  int ell,Delta;   // local-balance [constraint.txt]
  int Ns;          // block length (symbols)
  int Nb;          // block length (bits)
  int QA, Nu;      // numCW-A, symbol-len [ICB]
  int Nseq;        // compo: number of sequences (=r)
  int Kc, QB, QB2; // compo: resolution (=k), Num. symbol, QB2=QB/2;
  int Ndi;         // number of decoding iterations
  int Dmin,Dmax;   // drift min & max
  int seed;
  int t0=(int)time(NULL);
  double Pi,Pd,Ps; // IDS prob
  char *fn, *fncss;
  char *fncb    = new char [BSIZE];
  char *fnconst = new char [BSIZE];
  char *fncm    = new char [BSIZE];
  if(argc!=10){
    fprintf(stderr,"Usage: %s <ICB_dir> <CSS.txt> <Ns> <Pi> <Pd> <Ps> <Nseq> <NumDecItr> <seed|-1>\n",argv[0]);
    return 1;
  } // if
  fn   =      argv[1];
  fncss=      argv[2];
  Ns   = atoi(argv[3]);
  Pi   = atof(argv[4]);
  Pd   = atof(argv[5]);
  Ps   = atof(argv[6]);
  Nseq = atoi(argv[7]);
  Ndi  = atoi(argv[8]);
  seed = atoi(argv[9]);
  if(seed==-1) seed = (int)time(NULL);
  srandom(seed);
  assert(Ns>0 && Nseq>0);
  assert(Pi>=0.0 && Pi<0.5);
  assert(Pd>=0.0 && Pd<0.5);
  assert(Ps>=0.0 && Ps<0.5);
  snprintf(fncb,   BSIZE,"%s/cb.txt",        fn);  // inner codebook (in)
  snprintf(fnconst,BSIZE,"%s/constraint.txt",fn);  // constraints (in)
  snprintf(fncm,   BSIZE,"%s/EncCM.bin",     fn);  // encoding channel matrix (in)
  ReadConstraints(fnconst, &Rho, &ell, &Delta);
  class InnerCodebook *ICB = new class InnerCodebook(fncb,Rho,ell,Delta);
  QA   = ICB->Get_numCW();
  Nu   = ICB->Get_Nu();
  Nb   = Ns*Nu;
  Dmin = (int)floor(-(double)Drange*Nb*Pd);
  Dmax = (int)ceil ( (double)Drange*Nb*Pi);
  printf("# QA=%d Ns=%d Nu=%d Nb=%d (Pi,Pd,Ps)=(%e,%e,%e) Nseq=%d Ndi=%d [%d]\n",
	 QA,Ns,Nu,Nb,Pi,Pd,Ps,Nseq,Ndi,seed);
  printf("# ICB:   %s\n",fncb);
  printf("# Const: %s\n",fnconst);
  printf("# EncCM: %s\n",fncm);
  class ChannelMatrix *ECM  = new class ChannelMatrix(fncm);
  class CSFBAdec      *DECa = new class CSFBAdec(   Nb,Pi,Pd,Ps,Dmin,Dmax, ICB,ECM, Nseq, Ndi);
  class CPScodec      *DECb = new class CPScodec(   Nb,Pi,Pd,Ps,Dmin,Dmax, Nu, Nseq,fncss);
  Kc  = DECb->Get_Kc();
  QB  = DECb->Get_QB();
  QB2 = QB/2; 
  printf("# Kc=%d QB=%d QB2=%d\n",Kc,QB,QB2);
  assert( QB%2==0 );
  class CIDSchannel   *CCH  = new class CIDSchannel(Nb,Pi,Pd,Ps,Dmin,Dmax, Kc,Nseq);
  class ChannelMatrix *DCMa = new class ChannelMatrix(QA, (int)pow(QA, OutListSize));
  class ChannelMatrix *DCMb = new class ChannelMatrix(QB2,(int)pow(QB2,OutListSize));
  class ConcatAB      *CAB  = new class ConcatAB(ICB, DECa, DECb);
  int              *Nb2 = new int       [Nseq];                   // received word length (bits)
  int              *IWA = new int       [Ns];                     // information word A [Ns]
  int             **IWB = new int *     [Ns];                     // information word B [Ns][Nu]
  unsigned char    *CWA = new unsigned char [Nb];                 // codeword A [Nb]
  int             **CCW = new int *     [Nb];                     // composite codeword [Nb][4]
  unsigned char    **RW = new unsigned char * [Nseq];             // received word [Nseq][Nb+Dmax]
  unsigned char   **RW2 = new unsigned char * [Nseq];             //   (binarized -> DECa)
  double        **PoutA = new double *  [Ns];                     // output PMF A [Ns][QA]
  double        **PoutB = new double *  [Nb];                     // output PMF B [Nb][QB2]
  int              *DWA = new int       [Ns];                     // decoded word A [Ns]
  int             **DWB = new int *     [Ns];                     // decoded word B [Ns][Nu]
  long            *DWAL = new long      [Ns];                     // (PoutA->list->long) [Ns]
  long           **DWBL = new long *    [Ns];                     // (PoutB->list->long) [Ns][Nu]
  int           **dbgDV = new int *     [Nseq];                   // (dbg)drift vector [Nseq][Nb+1]
  unsigned char **dbgXV = new unsigned char * [Nseq];             // (dbg)channel input [Nseq][Nb]
  double        **dbgB2 = new double *  [Nb]; 
  for(int i=0; i<Ns;  i++) IWB[i]   = new int [Nu];
  for(int i=0; i<Nb;  i++) CCW[i]   = new int [4]; 
  for(int i=0; i<Nseq;i++) RW[i]    = new unsigned char [Nb + CCH->GetDmax()];
  for(int i=0; i<Nseq;i++) RW2[i]   = new unsigned char [Nb + CCH->GetDmax()];
  for(int i=0; i<Ns;  i++) PoutA[i] = new double [QA];
  for(int i=0; i<Nb;  i++) PoutB[i] = new double [QB2];
  for(int i=0; i<Ns;  i++) DWB[i]   = new int    [Nu];
  for(int i=0; i<Ns;  i++) DWBL[i]  = new long   [Nu];
  for(int i=0; i<Nseq;i++) dbgDV[i] = new int    [Nb+1];
  for(int i=0; i<Nseq;i++) dbgXV[i] = new unsigned char [Nb];
  for(int i=0; i<Nb;  i++) dbgB2[i] = new double [2];
  int wc;
  long eca, ecamax=0, esa=0;
  long ecb, ecbmax=0, esb=0;
  
  //-----
  for(wc=1;wc<=WCmax;wc++){
    // 1. random info
    RandVect(IWA,Ns,0,QA-1);
    for(int i=0;i<Ns;i++) RandVect(IWB[i],Nu,0,QB2-1);
    // 2. encode
    ICB->Encode(CWA,IWA,Ns);                      // encode-A
    DECb->Encode(CCW,(const int **)IWB,CWA);      // encode-B
    // 3. transmit
    CCH->transmit(RW,Nb2,(const int**)CCW);
    CCH->funcFd(RW2,(const unsigned char **)RW);  // binarize for dec-A
    // 4. decode-A 
    DECa->Decode(PoutA, (const unsigned char **)RW2, Nb2,IWA);
    HardDecision(DWA, (const double **)PoutA, Ns, QA);
    // 5. decode-B 
    CAB->AtoB();
    DECb->Decode(PoutB, (const unsigned char **)RW, Nb2, (const int **)IWB);
    HardDecision(DWB, (const double **)PoutB, Ns, Nu, QB2);
    // (count)
    OutputConv(DWAL, (const double **)PoutA, Ns, QA);
    OutputConv(DWBL, (const double **)PoutB, Ns, Nu, QB2);
    for(int i=0;i<Ns;i++){
      DCMa->countup(IWA[i],DWAL[i]);
      for(int j=0; j<Nu; j++) DCMb->countup(IWB[i][j], DWBL[i][j]);
    } // for i
    eca = HammingDist(IWA, DWA, Ns);
    ecb = HammingDist((const int **)IWB, (const int **)DWB, Ns, Nu);
    esa += eca;
    esb += ecb;
    ecamax = max(eca,ecamax);
    ecbmax = max(ecb,ecbmax);

    if( time(NULL)-t0 > Tint || wc==WCmax){
      printf("%04d : %ld/%ld=%e %ld : %ld/%ld=%e %ld : %e %e %e : %e %e %e\n", wc,
	     esa, (long)wc*Ns,    (double)esa/(wc*Ns),    ecamax, 
	     esb, (long)wc*Ns*Nu, (double)esb/(wc*Ns*Nu), ecbmax,  
	     DCMa->Hx(), DCMa->Hxy(), DCMa->Ixy(),
	     DCMb->Hx(), DCMb->Hxy(), DCMb->Ixy());
      t0 = (int)time(NULL);
    } // if wc 

    //(dbg)
    // CCH->GetDV(dbgDV);
    // CCH->GetXV(dbgXV);
    // CAB->Get_PrB2(dbgB2);
    // dbgPrint(IWA,(const int **)IWB,DWA,(const int **)DWB,CWA,(const int **)CCW,(const unsigned char **)RW,
    // 	     (const double**)PoutA, (const double**)PoutB,
    // 	     (const unsigned char **)dbgXV,(const int **)dbgDV,(const double **)dbgB2,Nb,Nb2,Nu,QA,QB2,Nseq);
  } // for wc
  //-----
  //DCMa->PrintCnt();

  delete ICB;
  delete ECM;
  delete CCH;
  delete DECa;
  delete DECb;
  delete DCMa;
  delete DCMb;
  delete CAB;
  delete [] IWA;
  delete [] CWA;
  delete [] DWA;
  delete [] DWAL;
  delete [] Nb2;
  delete [] fncb;
  delete [] fnconst;
  delete [] fncm;
  for(int i=0; i<Ns;  i++) delete [] IWB[i];
  for(int i=0; i<Nb;  i++) delete [] CCW[i];
  for(int i=0; i<Ns;  i++) delete [] PoutA[i];
  for(int i=0; i<Nb;  i++) delete [] PoutB[i];
  for(int i=0; i<Nseq;i++) delete [] RW[i];
  for(int i=0; i<Nseq;i++) delete [] RW2[i];
  for(int i=0; i<Ns;  i++) delete [] DWB[i];
  for(int i=0; i<Ns;  i++) delete [] DWBL[i];
  for(int i=0; i<Nseq;i++) delete [] dbgDV[i];
  for(int i=0; i<Nseq;i++) delete [] dbgXV[i];
  for(int i=0; i<Nb;  i++) delete [] dbgB2[i];
  delete [] IWB;
  delete [] CCW;
  delete [] PoutA;
  delete [] PoutB;
  delete [] RW;
  delete [] RW2;
  delete [] DWB;
  delete [] DWBL;
  delete [] dbgDV;
  delete [] dbgXV;
  delete [] dbgB2;
  return 0;
}

//================================================================================
void OutputConv(long *DWL, const double **P, int Ns, int Q){
  for(int i=0;i<Ns;i++){
    DWL[i] = PdistArgMaxLong(P[i], Q, OutListSize);
    //printf("%03d: %ld\n",i,DWL[i]);
    //PrintVect(P[i],Q," ","\n");
  } // for i
}

//================================================================================
void OutputConv(long **DWL, const double **P, int Ns, int Nu, int Q){
  for(int i=0; i<Ns; i++){
    for(int j=0; j<Nu; j++){
      DWL[i][j] = PdistArgMaxLong(P[i*Nu+j], Q, OutListSize);
    } // for j
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
void dbgPrint(const int *IWA, const int **IWB, const int *DWA, const int **DWB,  
	      const unsigned char *CWA, const int **CCW,
	      const unsigned char **RW, const double **PoutA, const double **PoutB,
	      const unsigned char **dbgXV, const int **dbgDV, const double **dbgB2,
	      int Nb, int *Nb2, int Nu, int QA, int QB2, int Nseq){
  int idx;
  int Nb2max = max(Nb2,Nseq);
  for(int i=0;i<max(Nb,Nb2max);i++){
    if(i%Nu==0){
      if(i<Nb){
	idx = i/Nu;
	printf("[%03d] %02d %02d\n",idx,IWA[idx],DWA[idx]);
	PrintVect(PoutA[idx],QA, "",       "\n");
      } else {
	printf("[---]\n");
      } // if i<Nb
    } // if i%Nu
    //---
    printf("%04d: ",i);
    //---
    if(i<Nb){
      printf("%u (%.2e,%.2e) %d %d (%d %d %d %d)  ",
	     CWA[i],dbgB2[i][0],dbgB2[i][1],IWB[i/Nu][i%Nu],DWB[i/Nu][i%Nu],
	     CCW[i][0],CCW[i][1],CCW[i][2],CCW[i][3]);
      for(int j=0;j<Nseq;j++) printf("%u",dbgXV[j][i]);
      printf(" ");
    } else {
      printf("- (--------,--------) - (- - - -)  ");
      for(int j=0;j<Nseq;j++) printf("-");
      printf(" ");
    } // if
    //---
    for(int j=0;j<Nseq;j++){
      if(i<Nb2[j]) printf("%u",RW[j][i]);
      else         printf("-");
    } // for j
    printf(" ");
    //---
    if(i<Nb+1){
      for(int j=0;j<Nseq;j++) printf("%+03d ",dbgDV[j][i]);
    } else {
      for(int j=0;j<Nseq;j++) printf("--- ");
    }
    //---
    if(i<Nb) PrintVect(PoutB[i], QB2, " ","");
    //---
    printf("\n");
  } // for i
}
