#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"


//================================================================================
void RandVect(int *V, int len, int Vmin, int Vmax){
  assert(Vmin<Vmax);
  int W = Vmax-Vmin+1;
  for(int i=0;i<len;i++){
    V[i] = (int)floor(((double)random()/RAND_MAX)*(double)W) + Vmin;
    assert(V[i]>=Vmin && V[i]<=Vmax);
  } // for i
}

//================================================================================
void PrintVect(const int *V, int len, const char *pre, const char *post){
  int w=40;
  printf("%s",pre);
  for(int i=0;i<len;i++){
    if(i%w==0) printf("%04d: ",i);
    printf("%d ",V[i]);
    if(i%w==w-1) printf("\n");
  } // for i
  printf("\n");
  printf("%s",post);
}

//================================================================================
void PrintVect(const double *V, int len, const char *pre, const char *post){
  int w=40;
  printf("%s",pre);
  for(int i=0;i<len;i++){
    if(i%w==0) printf("%04d: ",i);
    printf("%.2e ",V[i]);
    if(i%w==w-1) printf("\n");
  } // for i
  printf("\n");
  printf("%s",post);
}

//================================================================================
void PrintVect1(const double *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++){
    printf("%.2e ",V[i]);
  } // for i
  printf("%s",post);
}

//================================================================================
void PxUnif(double **Px, int N, int Q){
  assert(N>0 && Q>0);
  double p = (double)1.0/Q;
  for(int i=0;i<N;i++){
    for(int v=0;v<Q;v++){
      Px[i][v] = p;
    } // for v
  } // for i
}

//================================================================================
int max(int a, int b){return (a>=b)? a : b ;}
int min(int a, int b){return (a<=b)? a : b ;}
