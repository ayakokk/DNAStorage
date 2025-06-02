#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "bmatrix.hpp"
#include "func.hpp"


//=================================================================================
int main(int argc, char *argv[]){
  int  ell,w,delta,t0;
  int  gcb,mrl;
  bool flg;
  unsigned int w4p,u,cnt=0;
  if(argc!=4){
    fprintf(stderr,"Usage: %s <ell> <w> <delta>\n",argv[0]);
    return 1;
  }
  ell  = atoi(argv[1]);
  w    = atoi(argv[2]);
  delta= atoi(argv[3]);
  w4p  = (unsigned int)pow(4,w);
  printf("# ell=%d w=%d(%u) delta=%d\n",ell,w,w4p,delta);
  assert(ell>0);
  assert(w>0 && (long unsigned int)w*2 <= sizeof(unsigned int)*8); //(bits)
  assert(delta>=0);

  unsigned char *V = new unsigned char [w];
  t0 = (int)time(NULL);
  
  for(u=0;u<w4p;u++){
    ConvIntVect(u,V,w);
    gcb = GCbalance(V,w);
    mrl = MaxRunLength(V,w);
    flg = ((int)abs(gcb)<=delta) && (mrl<=ell);
    if(flg) cnt++;
    
    //-----
    if((int)time(NULL)-t0>10){
      printf("%08X: ",u);
      PrintVect(V,w);
      printf(": %+d\t%d\t%d\t(%u/%u)\n",gcb,mrl,flg,cnt,u);
    } // if
  } // for u
  printf("%u/%u\n",cnt,u);
  printf("rate: %e[bit]\n",log2(cnt)/w);
  
  delete [] V;
  return 0;
}
