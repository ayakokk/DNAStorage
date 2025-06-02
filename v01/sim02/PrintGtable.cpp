#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "Gtable.hpp"

#define BSIZE 4096

//=================================================================================
int main(int argc, char *argv[]){
  char *fn;
  if(argc!=2){
    fprintf(stderr,"Usage: %s <Gtable.bin>\n",argv[0]);
    return 1;
  }
  fn = argv[1];
  printf("# read from: %s\n",fn);
  class Gtable *GT = new class Gtable(fn);
  GT->dump();
  delete GT;
  return 0;
}
