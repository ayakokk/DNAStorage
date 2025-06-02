#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "InnerCodec.hpp"

//================================================================================
int main(int argc, char *argv[]){
  int ell, rlk, eps;
  if(argc!=4){
    fprintf(stderr,"Usage: %s <ell> <k> <eps>\n",argv[0]);
    return 1;
  } //if
  ell = atoi(argv[1]);
  rlk = atoi(argv[2]);
  eps = atoi(argv[3]);
  printf("# ell(win size):%d  rlk(max RL):%d  eps(balance):%d\n",ell,rlk,eps);
  assert(ell>0);
  assert(rlk>0 && rlk<ell);
  assert(eps>=0 && eps<ell);

  class InnerCodec *ICD = new class InnerCodec(ell,rlk,eps);
  
  delete ICD;
  return 0;
}
