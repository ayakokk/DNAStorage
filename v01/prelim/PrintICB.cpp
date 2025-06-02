#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "bmatrix.hpp"
#include "InnerCodebook.hpp"

//=================================================================================
int main(int argc, char *argv[]){
  if(argc!=2){
    fprintf(stderr,"Usage: %s <codebook.bmat>\n",argv[0]);
    return 1;
  }
  class InnerCodebook *ICB = new class InnerCodebook(argv[1]);
  ICB->dump();
  delete ICB;
  return 0;
}
