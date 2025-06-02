#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"

//================================================================================
int main(int argc, char *argv[]){
  int len;
  unsigned int u,u2p,u4p;
  if(argc!=2){
    fprintf(stderr,"Usage: %s <len>\n",argv[0]);
    return 1;
  }
  len = atoi(argv[1]);
  u2p = (unsigned long)pow(2,len);
  u4p = (unsigned long)pow(4,len);
  printf("# len=%d u2p=%u u4p=%u\n",len,u2p,u4p);
  unsigned char *V = new unsigned char [len];
  
  // binary
  for(u=0;u<u2p;u++){
    ConvIntV2(V,u,len);
    assert(ConvV2Int(V,len)==u);
  } // for u
  printf("binary: checked\n");

  // 4-ary
  for(u=0;u<u4p;u++){
    ConvIntV4(V,u,len);
    assert(ConvV4Int(V,len)==u);
  } // for u
  printf("4-ary:  checked\n");
  
  delete [] V;
  return 0;
}
