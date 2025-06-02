#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "NBIDSchannel.hpp"


//=================================================================================
void NBIDSchannel::ClearVect(unsigned char *V, int len){
  for(int i=0;i<len;i++) V[i]=0;
}

//=================================================================================
int NBIDSchannel::Rnd01(double p1){
  assert(p1>=0.0 && p1<=1.0);
  if( ((double)random()/RAND_MAX) < p1) return 1;
  else                                  return 0;
}

//=================================================================================
int NBIDSchannel::NextD(int d){
  assert(d>=-Dmax && d<=Dmax);
  if(Rnd01(2.0*Pid)==1){
    // ins or del
    if( (random()%2)==0 ){
      // ins
      if(d< Dmax) return d+1;
      else        return d;   // right-end
    } else {
      // del
      if(d>-Dmax) return d-1;
      else        return d;   // left-end
    } // if random()%2
  } else {
    // no ins/del
    return d;
  } // if Pid
}

//=================================================================================
unsigned char NBIDSchannel::substitution(unsigned char z){
  int ze;
  assert(z>=0 && z<Q);
  if(Rnd01(Ps)==1){
    ze = (int)z + ( 1+(random()%(Q-1)) );
    return (unsigned char) (ze%Q);
  } else {
    return z;
  }
}

//=================================================================================
//=================================================================================
//=================================================================================

//=================================================================================
NBIDSchannel::NBIDSchannel(int _N, int _Q, int _Dmax, double _Pid, double _Ps){
  N   = _N;
  Q   = _Q;
  Dmax= _Dmax;
  Pid = _Pid;
  Ps  = _Ps;
  printf("# NBIDSchannel: N=%d Q=%d Dmax=%d (Pid,Ps)=(%e,%e)\n",N,Q,Dmax,Pid,Ps);
  assert(N>0);
  assert(Q>0 && Q<=256);
  assert(Dmax>0);
  assert(Pid>=0.0 && Pid<0.5);
  assert(Ps >=0.0 && Ps <1.0); 
  DV = new int [N+1];
}

//=================================================================================
NBIDSchannel::~NBIDSchannel(){
  delete [] DV;
  printf("# NBIDSchannel: deleted\n");
}

//=================================================================================
int NBIDSchannel::transmit(unsigned char *Y, const unsigned char *X){
  unsigned char *Z = new unsigned char [N+Dmax];
  ClearVect(Y,N+Dmax);
  ClearVect(Z,N+Dmax);
  
  // drift vector -----
  DV[0]=0;
  for(int i=0;i<N;i++) DV[i+1] = NextD(DV[i]);

  // set Z -----
  for(int i=0;i<N;i++){
    for(int j=DV[i];j<=DV[i+1];j++) Z[i+j] = X[i];
  } // for i

  // set Y -----
  for(int i=0;i<N+DV[N];i++)
    Y[i] = substitution(Z[i]);
  
  //(dbg)
  /*
  for(int i=0;i<N+Dmax;i++){
    if(i<N) printf("%04d (%+03d) %X %X %X (%d)\n",i,DV[i],X[i],Z[i],Y[i],(Y[i]==Z[i]));
    else    printf("%04d (---) - %X %X (%d)\n",   i,           Z[i],Y[i],(Y[i]==Z[i]));
  } // for i
  */
  
  delete [] Z;
  return N+DV[N];
}

