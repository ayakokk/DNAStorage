#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"

#define BSIZE 8192
#define NuMax ((int)sizeof(int)*8-1)

bool CheckCond(int c);
void WriteCB(const bool *CW, const char *fn);

int Nu,Nu2p;
int rho,delta;
int wt;      // Nu/2 

//================================================================================
int main(int argc, char *argv[]){
  char *fn;
  if(argc!=5){
    fprintf(stderr,"Usage: %s <Nu> <rho> <delta> <cbout.txt>\n",argv[0]);
    return 1;
  }
  Nu    = atoi(argv[1]);
  rho   = atoi(argv[2]);
  delta = atoi(argv[3]);
  fn    = argv[4];
  Nu2p  = (int)pow(2,Nu);
  wt    = Nu/2;
  printf("# Nu=%d N2p=%d wt=%d rho=%d delta=%d (NuMax=%d)\n",Nu,Nu2p,wt,rho,delta,NuMax);
  printf("# output: %s\n",fn);
  assert(Nu>0 && Nu<=NuMax && Nu%2==0); 

  bool *CW = new bool [Nu2p];
  for(int c=0;c<Nu2p;c++) CW[c] = CheckCond(c);
  WriteCB(CW,fn);
  
  delete [] CW;
  return 0;
}

//================================================================================
bool CheckCond(int c){
  int w,len;
  bool ret = true;
  unsigned char *V = new unsigned char [Nu];
  LongToVect(V,c,Nu);

  // (1) w(x)=Nu/2
  if(HammingWeight(V,Nu)!=wt) ret = false;

  // (2) partial weight (even)
  if(ret){
    for(int i=0;i<Nu;i+=2){
      len = Nu-i;
      w = HammingWeight(&V[i],len);
      w = (2*w) - len;
      if(abs(w)>2*delta){
	ret = false;
	break;
      } // if abs(w)
    } // for i
  } // if ret
  
  // (3) partial weight (odd)
  if(ret){
    for(int i=1;i<Nu;i+=2){
      len = Nu-i;
      w = HammingWeight(&V[i],len);
      w = (2*w) - len;
      if(abs(w)>2*delta-1){
	ret = false;
	break;
      } // if abs(w)
    } // for i
  } // if ret

  // (4) max runlength <= rho
  if(ret){
    len = MaxRunLength(V,Nu);
    if(len>rho) ret = false;
    //printf("%d\n",len);
  } // if ret

  // (5) last runlength <= rho-1
  if(ret){
    len = RunLength(V,Nu-1,Nu);
    if(len>rho-1) ret=false;
    //printf("%d\n",len);
  } // if ret
  
  //(dbg)
  // if(ret){
  //   printf("%04X ",c);
  //   PrintVectB(V,Nu,"","\n");
  // } // if ret
  
  delete [] V;
  return ret;
}

//================================================================================
void WriteCB(const bool *CW, const char *fn){
  int  numCW;
  unsigned char *V = new unsigned char [Nu];
  FILE *fp;
  if((fp=fopen(fn,"w"))==NULL){
    fprintf(stderr,"Cannot open %s\n",fn);
    exit(1);
  } // if
  
  // L1: nu numCW
  numCW = HammingWeight(CW,Nu2p);
  fprintf(fp,"%d %d\n",Nu,numCW);

  // L2+: codewords
  for(int c=0;c<Nu2p;c++){
    if(CW[c]){
      LongToVect(V,c,Nu);
      for(int i=0;i<Nu;i++) fputc(V[i]+'0',fp);
      fprintf(fp,"\n");
      // (dbg)
      //PrintVectB(V,Nu,"","\n");
    } // if CW
  } // for c
  
  fclose(fp);
  delete [] V;
}
