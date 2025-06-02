#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"
#include "MutualInfo.hpp"
#include "InnerCodebook.hpp"

#define WCmax 1000000

void WritePx(const char *fn, const double *Px, int Nu2p);

//================================================================================
int main(int argc, char *argv[]){
  int Rho;         // run-length
  int ell,Delta;   // local-balance
  int N;           // block length (symbols)
  int Nb;          // block length (bits)
  int Q, Nu, Nu2p; // ICB
  int seed;
  char *fn, *fnPx;
  if(argc!=8){
    fprintf(stderr,"Usage: %s <ICB.txt> <Rho> <ell> <Delta> <N> <Px_out.txt> <seed|-1>\n",argv[0]);
    return 1;
  } // if
  fn    =      argv[1];
  Rho   = atoi(argv[2]);
  ell   = atoi(argv[3]);
  Delta = atoi(argv[4]);
  N     = atoi(argv[5]);
  fnPx  =      argv[6];
  seed  = atoi(argv[7]);
  if(seed==-1) seed = (int)time(NULL);
  srandom(seed);
  class InnerCodebook *ICB = new class InnerCodebook(fn,Rho,ell,Delta);
  Q    = ICB->Get_numCW();
  Nu   = ICB->Get_Nu();
  Nu2p = ICB->Get_Nu2p();
  Nb   = N*Nu;
  printf("# Q=%d Nu=%d(%d) Nb=%d [%d]\n",Q,Nu,Nu2p,Nb,seed);
  printf("# Input ICB: %s\n",fn);
  printf("# Output Px: %s\n",fnPx);
  class MutualInfo *MI0 = new class MutualInfo(Q,Nu2p);
  int *IV = new int [N];                       // information word
  unsigned char *CV = new unsigned char [Nb];  // codeword
  double *Px = new double [Nu2p];              // channel prior
  int wc;
  //-----
  for(wc=1;wc<=WCmax;wc++){
    RandVect(IV,N,0,Q-1);
    ICB->Encode(CV,IV,N);
    // count
    for(int i=0;i<N;i++){
      if(i*Nu>=ell) MI0->countup(IV[i],(int)VectToLong(&CV[i*Nu],Nu));
    } // for i
    // print
    if(wc%100000==0 || wc==WCmax-1){
      printf("wc=%07d H(D)=%e H(D|X)=%e I(D;X)=%e\n",wc,MI0->Hx(),MI0->Hxy(),MI0->Ixy());
    }
    //(dbg)
    //printf("wc=%d\n",wc);
    //MI0->PrintCnt();
    //PrintVectX(IV,N, "IV\n","\n");
    //PrintVectB(CV,Nb,"CV\n","\n");
  }
  //-----
  //MI0->PrintCnt();
  //printf("H(D)=%e H(D|X)=%e I(D;X)=%e\n",MI0->Hx(),MI0->Hxy(),MI0->Ixy());
  MI0->GetPy(Px);
  WritePx(fnPx,Px,Nu2p);
  
  delete ICB;
  delete MI0;
  delete [] IV;
  delete [] CV;
  return 0;
}

//================================================================================
void WritePx(const char *fn, const double *Px, int Nu2p){
  FILE *fp;
  if((fp=fopen(fn,"w"))==NULL){
    fprintf(stderr,"Cannot open %s\n",fn);
    exit(1);
  } // if
  fprintf(fp,"%d\n",Nu2p);
  for(int i=0;i<Nu2p;i++){
    if(Px[i]>0) fprintf(fp,"%d %e\n",i,Px[i]);
  } // for i
  fclose(fp);
}
