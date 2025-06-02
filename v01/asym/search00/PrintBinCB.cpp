#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "bmatrix.hpp"
#include "func.hpp"

//================================================================================
int main(int argc, char *argv[]){
  if(argc!=2){
    fprintf(stderr,"Usage: %s <BinCB.bin>\n",argv[0]);
    return 1;
  }
  class bmatrix *CB = new class bmatrix(argv[1],true);

  unsigned int l2p = (unsigned int)CB->getM();
  unsigned int b2p = (unsigned int)CB->getN();
  unsigned int s,v;
  int lmd = (int)log2(l2p);
  int beta= (int)log2(b2p);
  unsigned char *S = new unsigned char [lmd];
  unsigned char *V = new unsigned char [beta];
  int cnt;
  
  printf("l2p=%u b2p=%u : lmd=%d beta=%d\n",l2p,b2p,lmd,beta);
  assert(l2p==(unsigned int)pow(2,lmd));
  assert(b2p==(unsigned int)pow(2,beta));
  for(s=0;s<l2p;s++){
    ConvIntV2(S,s,lmd);
    printf("%04X:",s);
    PrintVect2(S,lmd,"",": ");
    cnt=0;
    for(v=0;v<b2p;v++){
      if(CB->getV(s,v)==1){
	ConvIntV2(V,v,beta);
	PrintVect2(V,beta,""," ");
	cnt++;
      } // if
    } // for v
    printf("[%d]\n",cnt);
  } // for s

  delete [] S;
  delete [] V;
  delete CB;
  return 0;
}
