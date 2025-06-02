#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "InnerCodebook.hpp"

//================================================================================
int main(int argc, char *argv[]){
  if(argc!=2){
    fprintf(stderr,"Usage: %s <ICB.txt>\n",argv[0]);
    return 1;
  } // if
  class InnerCodebook *ICB = new class InnerCodebook(argv[1]);

  delete ICB;
  return 0;
}
