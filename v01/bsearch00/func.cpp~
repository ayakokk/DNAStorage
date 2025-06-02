#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "func.hpp"

//============================================================================
void PrintVect(const unsigned char *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d ",V[i]);
  printf("%s",post);
}

//============================================================================
int binom(int n, int k){
  assert(n>=k);
  if(n==k || k==0) return 1;
  if(2*k>n) k=n-k; 
  return binom(n-1,k) + binom(n-1,k-1);
}
