#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "EnumVect.hpp"


//============================================================================
bool EnumVect::nextV(int pos){
  assert(pos>=0 && pos<N);
  if(V[pos]<NE[pos]-1){
    V[pos]++;
    return true;
  } else {
    if(pos==N-1){
      return false;
    } else {
      V[pos]=0;
      return nextV(pos+1);
    } // if pos
  } // if V[pos]
}

//============================================================================
void EnumVect::PrintVect(const int *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d ",V[i]);
  printf("%s",post);
}


//============================================================================
//============================================================================
//============================================================================

//============================================================================
EnumVect::EnumVect(int _N, const int *_NE){
  assert(_N>0);
  N = _N;
  NE= new int [N];
  V = new int [N];
  memcpy(NE,_NE,sizeof(int)*N);
  initV();
  PrintVect(NE,N,"# EnumVect: ","\n");
}

//============================================================================
EnumVect::~EnumVect(){
  delete [] NE;
  delete [] V;
  printf("# EnumVect: deleted\n");
}

//============================================================================
void EnumVect::initV(){
  for(int i=0;i<N;i++) V[i]=0;
}

//============================================================================
void EnumVect::getV(int *X){
  memcpy(X,V,sizeof(int)*N);
}

//============================================================================
bool EnumVect::nextV(){
  return nextV(0);
}
