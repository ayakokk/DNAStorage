#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

int binom(int n, int k);

//======================================================================
int main(int argc, char *argv[]){
  int ell;
  if(argc!=2){
    fprintf(stderr,"Usage: %s <WinSize>\n",argv[0]);
    return 1;
  }
  ell = atoi(argv[1]);
  printf("# ell=%d\n",ell);

  double p,s=0.0;
  for(int w=0;w<=ell;w++){
    p = binom(ell,w)/pow(2,ell);
    s += p*fabs((double)ell-2.0*w);
    //printf("w=%d p=%e\n",w,p);
  } // 
  printf("absolute GCB: %e\n",s);
  return 0;
}


//======================================================================
int binom(int n, int k){
  assert(n>=k && k>=0);
  if(k==0 || k==n){
    return 1;
  } else {
    if(2*k>n) k=n-k;
    return binom(n-1,k) + binom(n-1,k-1);
  }
}
