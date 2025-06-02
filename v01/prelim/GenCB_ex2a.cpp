#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "bmatrix.hpp"
#include "func.hpp"
#include "InnerCodebook.hpp"

//=================================================================================
int main(int argc, char *argv[]){
  int beta=4, B=7, nu=2;         // fixed param
  int beta4p=(int)pow(4,beta);
  char *fn;
  
  if(argc!=2){
    fprintf(stderr,"Usage: %s <ICBout.bin>\n",argv[0]);
    return 1;
  }
  fn = argv[1];
  printf("# beta=%d(%d) B=%d nu=%d (fixed param)\n",beta,beta4p,B,nu);
  printf("# output: %s\n",fn);
  class InnerCodebook *ICB = new class InnerCodebook(beta,B,nu);

  unsigned int B0[] = { 3, 5, 6, 9,10,12, 2, 4};
  unsigned int B1[] = { 3, 5, 6, 9,10,12,13,11};
  unsigned char *X0 = new unsigned char [beta];  // binary: B0 
  unsigned char *X1 = new unsigned char [beta];  // binary: B1
  unsigned char *Y  = new unsigned char [beta];  // binary: 0-15
  unsigned char *V0 = new unsigned char [beta];  // 4ary: X0Y
  unsigned char *V1 = new unsigned char [beta];  // 4ary: X1Y
  int v0,v1;
  
  for(unsigned int i=0;i<8;i++){
    ConvToBinVect(X0,B0[i],beta);
    ConvToBinVect(X1,B1[i],beta);
    printf(" X0:"); PrintVect(X0,beta); 
    printf(" X1:"); PrintVect(X1,beta); printf("\n"); 
    for(unsigned int j=0;j<16;j++){
      ConvToBinVect(Y,j,beta);
      for(int k=0;k<beta;k++){
	V0[k] = X0[k]*2 + Y[k];
	V1[k] = X1[k]*2 + Y[k];
      } // for k
      v0 = ConvVectInt(V0,beta);
      v1 = ConvVectInt(V1,beta);
      ICB->CWMset(v0,0,1);
      ICB->CWMset(v1,1,1);
      //printf("(%02u %02u) (%02u %02u)\n",B0[i],j,B1[i],j);
      //printf("   "); PrintVect4(V0,beta); printf("(%03d) ",v0); PrintVect4(V1,beta); printf("(%03d)\n",v1);
    } // for j
  } // for i
  ICB->GenMap();
  ICB->check();
  ICB->CWMwrite(fn);
  
  delete ICB;
  delete [] X0;
  delete [] X1;
  delete [] Y;
  delete [] V0;
  delete [] V1;
  return 0;
}
