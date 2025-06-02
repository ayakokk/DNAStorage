#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "bmatrix.hpp"
#include "func.hpp"
#include "InnerCodebook.hpp"

#define BSIZE 4096

int beta, nu;
int *B;
unsigned int beta4p, beta2p;
unsigned int *numBW;
unsigned int **BX;

void ReadBaseFile(const char *fn);

//=================================================================================
int main(int argc, char *argv[]){
  char *fnIn, *fnOut;
  if(argc!=3){
    fprintf(stderr,"Usage: %s <baseICB.txt> <ICBout.bin>\n",argv[0]);
    return 1;
  }
  fnIn = argv[1];
  fnOut= argv[2];
  ReadBaseFile(fnIn);
  beta4p=(int)pow(4,beta);
  beta2p=(int)pow(2,beta);
  printf("# beta=%d(%u,%u) nu=%d\n",beta,beta4p,beta2p,nu);
  PrintVect(B,    nu,"# B=    ","\n");
  PrintVect(numBW,nu,"# numBW=","\n");
  printf("# input:  %s\n",fnIn);
  printf("# output: %s\n",fnOut);
  
  class InnerCodebook *ICB = new class InnerCodebook(beta,B,nu);
  unsigned char *X = new unsigned char [beta];  // binary: BX 
  unsigned char *Y = new unsigned char [beta];  // binary: 0-beta2p
  unsigned char *V = new unsigned char [beta];  // 4ary: XY
  int v;

  for(int cb=0;cb<nu;cb++){
    printf("cb=%d\n",cb);
    for(unsigned int i=0;i<numBW[cb];i++){
      ConvToBinVect(X,BX[cb][i],beta);
      printf(" X:"); PrintVect(X,beta); printf("\n"); 
      for(unsigned int j=0;j<beta2p;j++){
	ConvToBinVect(Y,j,beta);
	for(int k=0;k<beta;k++) V[k] = X[k]*2 + Y[k];
	v = ConvVectInt(V,beta);
	ICB->CWMset(v,cb,1);
	//printf("(%02u %02u)\n",BX[cb][i],j);
	//printf("   "); PrintVect4(V,beta); printf("(%03d)\n",v);
      } // for j
    } // for i
  } // for cb
  ICB->GenMap();
  assert(ICB->check());
  ICB->CWMwrite(fnOut);
  
  delete ICB;
  delete [] X;
  delete [] Y;
  delete [] V;
  delete [] B;
  delete [] numBW;
  return 0;
}

//=================================================================================
void ReadBaseFile(const char *fn){
  FILE *fp;
  char *buf = new char [BSIZE];
  
  if((fp=fopen(fn,"r"))==NULL){
    fprintf(stderr,"Cannot open %s\n",fn);
  }
  // 1st line
  assert(fgets(buf,BSIZE,fp)!=NULL);
  beta = atoi(strtok(buf, " \t"));
  nu   = atoi(strtok(NULL," \t\n"));
  B    = new int [nu];
  numBW= new unsigned int [nu];

  // 2nd line
  assert(fgets(buf,BSIZE,fp)!=NULL);
  B[0] = atoi(strtok(buf," \t"));
  for(int i=1;i<nu;i++) B[i] = atoi(strtok(NULL," \t\n"));

  BX = new unsigned int * [nu];
  for(int i=0;i<nu;i++){
    numBW[i]= (int)pow(2,B[i]-beta);
    BX[i]= new unsigned int [numBW[i]];
  } // for i

  // 3nd+ lines 
  for(int i=0;i<nu;i++){
    assert(fgets(buf,BSIZE,fp)!=NULL);
    BX[i][0] = atoi(strtok(buf," \t"));
    for(unsigned int j=1;j<numBW[i];j++){
      BX[i][j] = atoi(strtok(NULL," \t\n"));
    } // for j
  } // for i
  
  fclose(fp);
  delete [] buf;
}
