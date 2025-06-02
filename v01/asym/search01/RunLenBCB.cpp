#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"
#include "bmatrix.hpp"

int beta,nu;
int *L0,*L1, *C0,*C1, *R0,*R1;  // max RL
unsigned int beta2p;
class bmatrix *BCB;

void SetLCR();
int  RunLen(const unsigned char *U, int pos0, int pos1, unsigned char v);  // pos0 -> pos1 
int  MaxRL( const unsigned char *U, int len, unsigned char v);
int  RunLenCB(int cb, unsigned char v);

//================================================================================
int main(int argc, char *argv[]){
  char *fn;
  if(argc!=2){
    fprintf(stderr,"Usage: %s <BCB.bin>\n",argv[0]);
    return 1;
  } // if
  fn = argv[1];
  printf("# BCB: %s\n",fn);
  BCB   = new class bmatrix(fn,true); 
  beta2p= BCB->getM();
  nu    = BCB->getN();
  beta  = (int)log2(beta2p);
  printf("# beta=%d beta2p=%u nu=%d\n",beta,beta2p,nu);
  assert(beta2p==(unsigned int)pow(2,beta));

  //---L0,L1, C0,C1, R0,R1: max run lengths
  L0 = new int [nu]; 
  L1 = new int [nu]; 
  C0 = new int [nu]; 
  C1 = new int [nu]; 
  R0 = new int [nu]; 
  R1 = new int [nu]; 
  SetLCR();

  //---calc
  printf("[Run length]\n");
  for(int cb=0;cb<nu;cb++){
    printf("CB%d: (0)%02d (1)%02d\n",cb,RunLenCB(cb,0),RunLenCB(cb,1));
  } // for cb
  
  //--- delete
  delete [] L0;
  delete [] L1;
  delete [] C0;
  delete [] C1;
  delete [] R0;
  delete [] R1;
  delete BCB;
  return 0;
}

//================================================================================
void SetLCR(){
  int l0,l1, c0,c1, r0,r1;
  unsigned int u;
  unsigned char *U = new unsigned char [beta];
  // --- init
  for(int cb=0;cb<nu;cb++){
    L0[cb]=0;
    L1[cb]=0;
    C0[cb]=0;
    C1[cb]=0;
    R0[cb]=0;
    R1[cb]=1;
  } // for cb
  
  // --- set
  for(u=0;u<beta2p;u++){
    ConvIntV2(U,u,beta);
    l0 = RunLen(U,0,beta-1,0);
    l1 = RunLen(U,0,beta-1,1);
    r0 = RunLen(U,beta-1,0,0);
    r1 = RunLen(U,beta-1,0,1);
    c0 = MaxRL(U,beta,0);
    c1 = MaxRL(U,beta,1);
    for(int cb=0;cb<nu;cb++){
      if(BCB->getV(u,cb)==1){
	if(l0>L0[cb]) L0[cb]=l0;
	if(l1>L1[cb]) L1[cb]=l1;
	if(r0>R0[cb]) R0[cb]=r0;
	if(r1>R1[cb]) R1[cb]=r1;
	if(c0>C0[cb]) C0[cb]=c0;
	if(c1>C1[cb]) C1[cb]=c1;
      } // if BCB
    } // for cb
    
    //(dbg)
    //PrintVect2(U,beta,""," ");
    //printf("%d:%d %d:%d %d:%d\n",l0,l1,c0,c1,r0,r1);
  } // for u

  //(dbg)
  for(int cb=0;cb<nu;cb++){
    printf("%d: %d:%d %d:%d %d:%d\n",cb,L0[cb],L1[cb],C0[cb],C1[cb],R0[cb],R1[cb]);
  } // for cb
  
  delete [] U;
}

//================================================================================
int RunLen(const unsigned char *U, int pos0, int pos1, unsigned char v){
  assert(pos0>=0 && pos1>=0);
  assert(v==0 || v==1);
  int len=0;
  if(pos0<=pos1){
    for(int i=pos0;i<=pos1;i++){
      if(U[i]!=v) break;
      len++;
    } // for i
  } else {
    for(int i=pos0;i>=pos1;i--){
      if(U[i]!=v) break;
      len++;
    } // for i
  } // if pos0  
  return len;
}

//================================================================================
int MaxRL( const unsigned char *U, int len, unsigned char v){
  assert(len>0);
  assert(v==0 || v==1);
  int l,lmax=0;
  for(int i=0;i<len;i++){
    l = RunLen(U,i,len-1,v);
    if(l>lmax) lmax=l;
  } // for i
  return lmax;
}

//================================================================================
int RunLenCB(int cb, unsigned char v){
  assert(cb>=0 && cb<nu);
  assert(v==0 || v==1);
  int i, len0, len1;
  if(v==0){
    len0 = C0[cb];
    len1 = R0[cb];
    for(i=1;i<cb;i++){
      if(L0[(cb+i)%nu]<beta) break;
      len1 += beta;
    }
    len1 += L0[(cb+i)%nu]; 
  } else {
    len0 = C1[cb];
    len1 = R1[cb];
    for(i=1;i<cb;i++){
      if(L1[(cb+i)%nu]<beta) break;
      len1 += beta;
    }
    len1 += L1[(cb+i)%nu]; 
  } 
  return max(len0,len1);
}
