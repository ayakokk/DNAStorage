#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "bmatrix.hpp"
#include "func.hpp"
#include "InnerCodebook.hpp"

//=================================================================================
int main(int argc, char *argv[]){
  int  beta,B,nu,omg,lmd,seed;
  char *fn;
  if(argc!=8){
    fprintf(stderr,"Usage: %s <beta> <B> <nu> <omg> <lmd> <out.bmat> <seed|-1>\n",argv[0]);
    return 1;
  }
  beta = atoi(argv[1]);
  B    = atoi(argv[2]);
  nu   = atoi(argv[3]);
  omg  = atoi(argv[4]);
  lmd  = atoi(argv[5]);
  fn   =      argv[6];
  seed = atoi(argv[7]);
  if(seed==-1) seed = (int)time(NULL);
  srandom(seed);
  printf("# beta=%d B=%d nu=%d omg=%d lmd=%d [%d]\n",beta,B,nu,omg,lmd,seed);
  printf("# output: %s\n",fn);
  assert(beta>0 && beta<(int)sizeof(int)*2);
  assert(B>0 && B<2*beta);
  assert(nu>0 && omg>0 && lmd>1);
  
  class InnerCodebook *ICB = new class InnerCodebook(beta,B,nu);

  // forbidden set
  int beta4p = ICB->Get_beta4p();
  int B2p    = ICB->Get_B2p();
  int gcb,mrl;
  int Fsize;
  unsigned char *V = new unsigned char [beta];
  class bmatrix *F = new class bmatrix(beta4p,1); // forbidden
  F->clear();
  for(int i=0;i<beta4p;i++){
    ConvIntVect(i,V,beta);
    gcb = GCbalance(V,beta);
    mrl = MaxRunLength(V,beta);
    if(abs(gcb)>=omg || mrl>=lmd) F->setV(i,0,1);
    // printf("%04X: ",i);
    // PrintVect4(V,beta);
    // printf(": %+d %d: %d\n",gcb,mrl,F->getV(i,0));
  } // for i
  Fsize = F->HammingWeight();
  printf("# num forbidden patterns: %d\n",Fsize);
  assert(beta4p-Fsize>=B2p);

  // select CW (random)
  class bmatrix *CW = new class bmatrix(beta4p,1);
  for(int i=0;i<nu;i++){
    CW->clear();
    RandomSel(CW,B2p,F);
    ICB->CWMset(CW,i);
  } // for i
  
  // output
  ICB->GenMap();
  ICB->check();
  ICB->CWMwrite(fn);
  
  // (dbg)
  //ICB->dump();
  /*
  for(int i=0;i<beta4p;i++){
    printf("%04X[%d] %d\n",i,F->getV(i,0),CW->getV(i,0));
  }
  ICB->CWMwrite("tmp/tmp.bin");
  class InnerCodebook *ICB2 = new class InnerCodebook("tmp/tmp.bin");
  delete ICB2;
  */
  
  delete [] V;
  delete F;
  delete CW;
  delete ICB;
  return 0;
}
