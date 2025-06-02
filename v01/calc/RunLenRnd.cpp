#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

//======================================================================
int main(int argc, char *argv[]){
  int q,n;
  if(argc!=3){
    fprintf(stderr,"Usage: %s <q> <n>\n",argv[0]);
    return 1;
  }
  q = atoi(argv[1]);
  n = atoi(argv[2]);
  printf("# q=%d n=%d\n",q,n);
  assert(q>=2);
  assert(n>0);

  // calc-1
  double e1,s=0.0;
  for(int i=n;i>0;i--) s+=(double)i/pow(q,i);
  e1 = (double)(q-1)*s;

  // calc-2
  double e2;
  e2 = (double)q/(q-1)*(1.0-1.0/pow(q,n))-((double)n/pow(q,n));
  
  printf("e1=%e e2=%e\n",e1,e2);
  return 0;
}
