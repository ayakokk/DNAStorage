#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "InnerCodebook.hpp"

//================================================================================
int main(int argc, char *argv[]){
  int Rho;       // run-length
  int ell,Delta; // local-balance
  char *fn;
  if(argc!=5){
    fprintf(stderr,"Usage: %s <ICB.txt> <Rho> <ell> <Delta>\n",argv[0]);
    return 1;
  } // if
  fn    = argv[1];
  Rho   = atoi(argv[2]);
  ell   = atoi(argv[3]);
  Delta = atoi(argv[4]);
  class InnerCodebook *ICB = new class InnerCodebook(fn,Rho,ell,Delta);

  delete ICB;
  return 0;
}
