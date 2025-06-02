#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"
#include "InnerCodebook.hpp"

//================================================================================
int main(int argc, char *argv[]){
  int Rho;       // run-length
  int ell,Delta; // local-balance
  int N;         // block length (symbols)
  int Nb;        // block length (bits)
  int Q, Nu;     // ICB
  int seed;
  char *fn;
  if(argc!=7){
    fprintf(stderr,"Usage: %s <ICB.txt> <Rho> <ell> <Delta> <N> <seed|-1>\n",argv[0]);
    return 1;
  } // if
  fn    = argv[1];
  Rho   = atoi(argv[2]);
  ell   = atoi(argv[3]);
  Delta = atoi(argv[4]);
  N     = atoi(argv[5]);
  seed  = atoi(argv[6]);
  if(seed==-1) seed = (int)time(NULL);
  srandom(seed);
  class InnerCodebook *ICB = new class InnerCodebook(fn,Rho,ell,Delta);
  Q  = ICB->Get_numCW();
  Nu = ICB->Get_Nu();
  Nb = N*Nu;
  printf("# Q=%d Nu=%d Nb=%d\n",Q,Nu,Nb);
  int *IV = new int [N];
  unsigned char *CV = new unsigned char [Nb];

  //-----
  RandVect(IV,N,0,Q-1);
  ICB->Encode(CV,IV,N);

  //(dbg)
  PrintVectX(IV,N, "IV\n","\n");
  PrintVectB(CV,Nb,"CV\n","\n");
  //-----
  
  delete ICB;
  delete [] IV;
  delete [] CV;
  return 0;
}
