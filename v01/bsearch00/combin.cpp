#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "combin.hpp"


//=====================================================================
int combin::HammingWeight(const unsigned char *X, int len){
  int i,w=0;
  for(i=0;i<len;i++)
    if(X[i]!=0) w++;
  return w;
}

//=====================================================================
//=====================================================================
//=====================================================================

//=====================================================================
combin::combin(unsigned int _N, unsigned int _K){
  N = _N;
  K = _K;
  assert(N>=K);
  Num = (unsigned long)boost::math::binomial_coefficient<double>(N,K);
  //printf("# combin: N=%d K=%d Num=%lu\n",N,K,Num);

  V = new unsigned char [N];
  initV();
}

//=====================================================================
combin::~combin(){
  delete [] V;
  //printf("# combin: deleted\n");
}

//=====================================================================
unsigned long combin::getNum(){ return Num; }

//=====================================================================
void combin::initV(){
  for(int i=0;i<N;i++) V[i]=0;
  for(int i=0;i<K;i++) V[i]=1;
}

//=====================================================================
bool combin::nextV(){
  int i,j,cnt;
  
  for(i=0;i<K;i++){
    cnt=-1;
    for(j=0;j<N;j++){
      if(V[j]==1) cnt++;
      if(cnt==i) break;
    } // for j
    assert(j<N);
    assert(V[j]==1);
    if(j==N-1) return false;
    if(j<N-1){
      if(V[j+1]==0){
	// move to right
	V[j]=0;
	V[j+1]=1;
	assert(HammingWeight(V,N)==K);
	return true;
      } else {
	// reset to left-end
	V[j]=0;
	V[i]=1;
      }
    } // if j<N-1 
  } // for i

  assert(K==0);
  return false;
}

//=====================================================================
void combin::getV(unsigned char *X){
  memcpy(X,V,sizeof(unsigned char)*N);
}

//=====================================================================
void combin::printV(){
  for(int i=0;i<N;i++) printf("%u",V[i]);
  printf("\n");
}
