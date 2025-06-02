#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "InnerCodebook.hpp"

#define BSIZE 4096

//================================================================================
void InnerCodebook::PrintVect(const unsigned char *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%u",V[i]);
  printf("%s",post);
}

//================================================================================
bool InnerCodebook::IsEqual(const unsigned char *V0, const unsigned char *V1, int len){
  assert(len>0);
  for(int i=0;i<len;i++){
    assert(V0[i]==0 || V0[i]==1);
    assert(V1[i]==0 || V1[i]==1);
    if(V0[i]!=V1[i]) return false;
  } // for i
  return true;
}

//================================================================================
bool InnerCodebook::IsInvEqual(const unsigned char *V0, const unsigned char *V1, int len){
  assert(len>0);
  for(int i=0;i<len;i++){
    assert(V0[i]==0 || V0[i]==1);
    assert(V1[i]==0 || V1[i]==1);
    if(V0[i]==V1[i]) return false;
  } // for i
  return true;
}

//================================================================================
void InnerCodebook::ReadFile(const char *fn){
  FILE *fp;
  char *buf = new char [BSIZE];
  if((fp=fopen(fn,"r"))==NULL){
    fprintf(stderr,"Cannot open %s\n",fn);
    exit(1);
  } // if
  // L1
  assert(fgets(buf,BSIZE,fp)!=NULL);
  N     = atoi(strtok(buf, " \t"));
  numCW = atoi(strtok(NULL," \t\n")); 
  assert(N>0 && numCW>0);
  CW = new unsigned char * [numCW];
  for(int i=0;i<numCW;i++) CW[i] = new unsigned char [N];
  // L2+
  for(int i=0;i<numCW;i++){
    assert(fgets(buf,BSIZE,fp)!=NULL);
    for(int j=0;j<N;j++){
      CW[i][j] = buf[j]-'0';
      assert(CW[i][j]==0 || CW[i][j]==1);
    } // for j
  } // for i
  fclose(fp);
  delete [] buf;
}

//================================================================================
void InnerCodebook::SetFlg(){
  bool flg;
  // FlgUnique
  FlgUnique = true;
  for(int i=0;i<numCW;i++){
    for(int j=i+1;j<numCW;j++){
      if(IsEqual(CW[i],CW[j],N)){
	FlgUnique = false;
	break;
      } // if
    } // for j
    if(!FlgUnique) break;
  } // for i
  // FlgInvertible
  FlgInvertible = true;
  for(int i=0;i<numCW;i++){
    flg = false;
    for(int j=0;j<numCW;j++){
      if(IsInvEqual(CW[i],CW[j],N)){
	flg = true;
	break;
      } // if
    } // for j
    if(!flg){
      FlgInvertible = false;
      break;
    } // if
  } // for i
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
InnerCodebook::InnerCodebook(const char *fn){
  ReadFile(fn);
  SetFlg();
  printf("# InnerCodebook: input %s\n",fn);
  printf("# InnerCodebook: N=%d numCW=%d FlgUnique=%d FlgInvertivle=%d\n"
	 ,N,numCW,FlgUnique,FlgInvertible);

  //(dbg)
  PrintCodebook();
}

//================================================================================
InnerCodebook::~InnerCodebook(){
  for(int i=0;i<numCW;i++) delete [] CW[i];
  delete [] CW;
  printf("# InnerCodebook: deleted\n");
}

//================================================================================
void InnerCodebook::PrintCodebook(){
  for(int i=0;i<numCW;i++){
    printf("%03d: ",i);
    PrintVect(CW[i],N,"","\n");
  } // for i
}

//================================================================================
bool InnerCodebook::GetFlgUnique(){    return FlgUnique;}
bool InnerCodebook::GetFlgInvertible(){return FlgInvertible;}
