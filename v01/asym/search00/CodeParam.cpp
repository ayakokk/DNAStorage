#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"
#include "CodeParam.hpp"

#define BSIZE 8192


//================================================================================
//================================================================================
//================================================================================

//================================================================================
CodeParam::CodeParam(const char *fn){
  FILE *fp;
  char *buf = new char [BSIZE];
  
  printf("# CodeParam: (input) %s\n",fn);
  if((fp=fopen(fn,"r"))==NULL){
    fprintf(stderr,"Cannot open %s\n",fn);
    exit(1);
  } // if

  // L1: beta nu
  assert(fgets(buf,BSIZE,fp)!=NULL);
  beta = atoi(strtok(buf, " ,\t"));
  nu   = atoi(strtok(NULL," ,\t\n"));
  assert(beta>0 && nu>0);

  // L2: B
  B = new int [nu];
  assert(fgets(buf,BSIZE,fp)!=NULL);
  B[0] = atoi(strtok(buf," ,\t\n"));
  for(int i=1;i<nu;i++) B[i] = atoi(strtok(NULL," ,\t\n"));

  
  // ---
  printf("# CodeParam: beta=%d nu=%d\n",beta,nu);
  PrintVect(B,nu,"# CodeParam: B=","\n");
  
  fclose(fp);
  delete [] buf;
}

//================================================================================
CodeParam::~CodeParam(){
  delete [] B;
  printf("# CodeParam: deleted\n");
}


//================================================================================
int    CodeParam::get_beta(){return beta;}
int    CodeParam::get_nu(){  return nu;}
void   CodeParam::get_B(int *_B){memcpy(_B,B,nu*sizeof(int));}
