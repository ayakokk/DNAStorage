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
  int  mint, beta, mval;  // beta=marker_length
  int  beta4p,nu;
  char *fn;
  if(argc!=5){
    fprintf(stderr,"Usage: %s <interval> <marker_len> <marker_val> <output.bin>\n",argv[0]);
    return 1;
  }
  mint = atoi(argv[1]);
  beta = atoi(argv[2]);
  mval = atoi(argv[3]);
  fn   =      argv[4];
  assert(beta>0);
  assert(mint>0 && mint%beta==0);
  beta4p = (int)pow(4,beta);
  nu = mint/beta + 1;   // info + marker
  printf("# marker: interval=%d length=%d val=%d\n",mint,beta,mval);
  printf("# beta4p=%d nu=%d\n",beta4p,nu);
  printf("# output: %s\n",fn);
  assert(mval>=0 && mval<beta4p);
  
  int *B = new int [nu];
  for(int i=0;i<nu-1;i++) B[i] = beta*2; // info
  B[nu-1] = 0; // marker
  PrintVect(B,nu,"# B=","\n");

  class InnerCodebook *ICB = new class InnerCodebook(beta,B,nu);
  ICB->CWMclear();
  for(int cb=0;cb<nu-1;cb++){
    for(int i=0;i<beta4p;i++) ICB->CWMset(i,cb,1);
  } // for cb
  ICB->CWMset(mval,nu-1,1);
  ICB->GenMap();
  //ICB->dump();
  assert(ICB->check());
  ICB->CWMwrite(fn);
  
  delete ICB;
  delete [] B;
}
