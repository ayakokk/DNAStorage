#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "BaseFile.hpp"

#define BSIZE 8192


//============================================================================
void BaseFile::PrintVect(const int *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d ",V[i]);
  printf("%s",post);
}


//============================================================================
//============================================================================
//============================================================================

//============================================================================
BaseFile::BaseFile(const char *fn){
  FILE *fp;
  char *buf = new char [BSIZE];
  
  if((fp=fopen(fn,"r"))==NULL){
    fprintf(stderr,"Cannot open %s\n",fn);
    exit(1);
  }

  // L1: beta nu
  assert( fgets(buf,BSIZE,fp) != NULL);
  beta = atoi(strtok(buf, " ,\t"));
  nu   = atoi(strtok(NULL," ,\t\n"));
  assert(beta>0 && nu>0);
  
  // L2: B
  B =  new int [nu];
  B2p= new int [nu];
  assert( fgets(buf,BSIZE,fp) != NULL);
  B[0] = atoi(strtok(buf," ,\t\n"));
  for(int i=1;i<nu;i++) B[i] = atoi(strtok(NULL," ,\t\n"));
  for(int i=0;i<nu;i++){
    assert(B[i]>=beta && B[i]<=2*beta);
    B2p[i] = (int)pow(2,B[i]-beta);
  } // for i

  // L3+: CW
  CW = new int * [nu];
  for(int idx=0;idx<nu;idx++){
    CW[idx] = new int [B2p[idx]];
    assert( fgets(buf,BSIZE,fp) != NULL);
    CW[idx][0] = atoi(strtok(buf," ,\t\n"));
    for(int i=1;i<B2p[idx];i++)
      CW[idx][i] = atoi(strtok(NULL," ,\t\n"));
  } // for idx

  fclose(fp);
  delete [] buf;

  /*
  printf("# BaseFile: %s [beta=%d nu=%d]\n",fn,beta,nu);
  PrintVect(B,  nu,"B:","\n");
  PrintVect(B2p,nu,"B:","\n");
  for(int i=0;i<nu;i++) PrintVect(CW[i],B2p[i],"CW:","\n");
  */
}

//============================================================================
BaseFile::BaseFile(int _beta, int _nu, const int *_B){
  beta = _beta;
  nu   = _nu;
  B  = new int [nu];
  B2p= new int [nu];
  memcpy(B,_B,sizeof(int)*nu);
  assert(beta>0 && nu>0);
  for(int i=0;i<nu;i++){
    assert(B[i]>=beta && B[i]<=2*beta);
    B2p[i] = (int)pow(2,B[i]-beta);
  } // for i

  CW = new int * [nu];
  for(int i=0;i<nu;i++) CW[i] = new int [B2p[i]];
}

//============================================================================
BaseFile::~BaseFile(){
  delete [] B;
  delete [] B2p;
  for(int i=0;i<nu;i++) delete [] CW[i];
  delete [] CW;
}


//============================================================================
int  BaseFile::get_beta(){return beta; }
int  BaseFile::get_nu(){  return nu;   }
void BaseFile::get_B(  int *_B){  memcpy(_B,  B,  sizeof(int)*nu); }
void BaseFile::get_B2p(int *_B2p){memcpy(_B2p,B2p,sizeof(int)*nu); }

//============================================================================
void BaseFile::get_CW(int idx, int *_CW){
  assert(idx>=0 && idx<nu);
  memcpy(_CW,CW[idx],sizeof(int)*B2p[idx]);
}

//============================================================================
void BaseFile::set_CW(int idx, const int *_CW){
  assert(idx>=0 && idx<nu);
  memcpy(CW[idx],_CW,sizeof(int)*B2p[idx]);
}

//============================================================================
void BaseFile::write(const char *fn){
  FILE *fp;

  if((fp=fopen(fn,"w"))==NULL){
    fprintf(stderr,"Cannot open %s\n",fn);
    exit(1);
  }

  // L1: beta nu
  fprintf(fp,"%d %d\n",beta,nu);

  // L2: B
  for(int i=0;i<nu-1;i++) fprintf(fp,"%d ",B[i]);
  fprintf(fp,"%d\n",B[nu-1]);

  // L3+: CW[]
  for(int idx=0;idx<nu;idx++){
    for(int i=0;i<B2p[idx]-1;i++) fprintf(fp,"%d ",CW[idx][i]);
    fprintf(fp,"%d\n",CW[idx][B2p[idx]-1]);
  } // for idx
  
  fclose(fp);
}
