#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "func.hpp"

#define BSIZE 8192

//================================================================================
int max(int a, int b){return (a>b)? a : b;}
int min(int a, int b){return (a<b)? a : b;}

//================================================================================
int argmax(const double *V, int len){
  assert(len>0);
  int ret=0;
  for(int i=1;i<len;i++){
    if(V[i]>V[ret]) ret=i;
  } // for i
  return ret;
}

//================================================================================
int argmin(const double *V, int len){
  assert(len>0);
  int ret=0;
  for(int i=1;i<len;i++){
    if(V[i]<V[ret]) ret=i;
  } // for i
  return ret;
}

//================================================================================
double max(const double *V, int len){return V[ argmax(V,len) ]; }
double min(const double *V, int len){return V[ argmin(V,len) ]; }

//================================================================================
void PrintVect(const double *V, int len, const char *pre, const char *post){
  int w=20;
  printf("%s",pre);
  for(int i=0;i<len;i++){
    printf("%02d:%.2e ",i,V[i]);
    if(i%w==w-1) printf("\n");
  } // for i
  printf("%s",post);
}

//================================================================================
void PrintVectX(const int *V, int len, const char *pre, const char *post){
  int w=20;
  printf("%s",pre);
  for(int i=0;i<len;i++){
    if(i%w==0) printf("%04d: ",i);
    printf("%04X ",V[i]);
    if(i%w==w-1) printf("\n");
  } // for i
  printf("%s",post);
}

//================================================================================
void PrintVectB(const unsigned char *V, int len, const char *pre, const char *post){
  int w=100;
  printf("%s",pre);
  for(int i=0;i<len;i++){
    if(i%w==0) printf("%04d: ",i);
    printf("%u",V[i]);
    if(i%w==w-1) printf("\n");
  } // for i
  printf("%s",post);
}

//================================================================================
void RandVect(int *V, int len, int Vmin, int Vmax){
  assert(len>0);
  assert(Vmin<Vmax);
  int A = Vmax-Vmin+1;
  for(int i=0;i<len;i++){
    V[i] = random()%A;
    V[i] += Vmin;
  } // for i
}

//================================================================================
long VectToLong(const unsigned char *V, int len){
  assert(len>0);
  long val = 0;
  for(int i=0;i<len;i++){
    assert(V[i]==0 || V[i]==1);
    val <<= 1;
    if(V[i]==1) val |= 0x1;
  } // for i
  return val;
}

//================================================================================
void LongToVect(unsigned char *V, long val, int len){
  assert(len>0 && val>=0);
  long mask = 0x1 << (len-1);
  for(int i=0;i<len;i++){
    V[i] = ( (val & mask)==0 )? 0 : 1;
    mask >>= 1;
  } // for i
}

//================================================================================
void ReadConstraints(const char *fn, int *Rho, int *ell, int *Delta){
  FILE *fp;
  char *buf = new char [BSIZE];
  if((fp=fopen(fn,"r"))==NULL){
    fprintf(stderr,"Cannot open %s\n",fn);
    exit(1);
  } // if
  assert( fgets(buf,BSIZE,fp)!=NULL );
  (*Rho)   = atoi(strtok(buf, " \t"));
  (*ell)   = atoi(strtok(NULL," \t"));
  (*Delta) = atoi(strtok(NULL," \t\n"));
  fclose(fp);
  assert( (*Rho)  >0 );
  assert( (*ell)  >0 );
  assert( (*Delta)>0 );
}

//================================================================================
void HardDecision(int *V, const double **P, int N, int Q){
  assert( N>0 && Q>0 );
  for(int i=0;i<N;i++) V[i] = argmax( P[i], Q );
}

//================================================================================
int HammingDist(const int *V0, const int *V1, int len){
  assert(len>0);
  int cnt=0;
  for(int i=0;i<len;i++){
    if(V0[i]!=V1[i]) cnt++;
  } // for i
  return cnt;
}
