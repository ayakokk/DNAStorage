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
  int    beta;
  double Pid,Ps,Pth;
  char   *fn;
  char   *fn2 = new char [BSIZE];
  if(argc!=6){
    fprintf(stderr,"Usage: %s <beta> <Pid> <Ps> <Pth> <output.bin>\n",argv[0]);
    return 1;
  }
  beta = atoi(argv[1]);
  Pid  = atof(argv[2]);
  Ps   = atof(argv[3]);
  Pth  = atof(argv[4]);
  fn   =      argv[5];
  printf("# beta=%d (Pid,Ps)=(%e,%e) Pth=%e\n",beta,Pid,Ps,Pth);
  printf("# output: %s\n",fn);

  // generate
  class Gtable *GT = new class Gtable(beta,Pid,Ps,Pth);
  GT->write(fn);

  // read and re-write (check)
  snprintf(fn2,BSIZE,"%s.chk",fn);
  class Gtable *GT2 = new class Gtable(fn);
  GT2->write(fn2);
  
  delete GT;
  delete GT2;
  delete [] fn2;
  return 0;
}
