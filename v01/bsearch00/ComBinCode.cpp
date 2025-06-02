#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <assert.h>

#include "func2.hpp"
#include "BaseFile.hpp"
#include "EnumVect.hpp"

#define BSIZE 8192

void GetNumFiles(int *NumFiles, const char **dn, int nu);
void GetFileNames(char ***fn, const char **dn, int nu);
void Generate(const char ***fn, const int *sel, int nu, const char *odir, int cnt);
void PrintFileNames(const char ***fn, const int *NumFiles, int nu);

//============================================================================
int main(int argc, char *argv[]){
  int nu,cnt;
  char *odir;
  if(argc<=3){
    fprintf(stderr,"Usage: %s <outdir> <dir1> <dir2> ...\n",argv[0]);
    return 1;
  }
  nu = argc-2;
  odir = argv[1];
  printf("# nu=%d\n",nu);
  printf("# output dir: %s\n",odir);
  
  // set file names -----
  char **dn = new char * [nu];
  for(int i=0;i<nu;i++){
    dn[i] = argv[i+2];
    printf("# dn[%d]: %s\n",i,dn[i]);
  } // for i

  int *NumFiles = new int [nu];
  GetNumFiles(NumFiles,(const char **)dn,nu);
  PrintVect(NumFiles,nu,"# NumFiles: ","\n");
  char ***fn = new char ** [nu];
  for(int i=0;i<nu;i++){
    fn[i] = new char * [NumFiles[i]];
    for(int j=0;j<NumFiles[i];j++) fn[i][j] = new char [BSIZE];
  } // for i
  GetFileNames(fn,(const char **)dn,nu);
  PrintFileNames((const char ***)fn, (const int *)NumFiles, nu);

  // generate -----
  class EnumVect *EV = new class EnumVect(nu,NumFiles);
  int *sel = new int [nu];
  cnt=0;
  do {
    EV->getV(sel);
    Generate((const char ***)fn, (const int *)sel, nu, odir, cnt);
    cnt++;
    //PrintVect(sel,nu,"","\n");
  } while( EV->nextV() );

  
  // delete -----
  delete EV;
  delete [] sel;
  for(int i=0;i<nu;i++){
    for(int j=0;j<NumFiles[i];j++) delete [] fn[i][j];
    delete [] fn[i];
  } // for i
  delete [] dn;
  delete [] NumFiles;
  delete [] fn;
  return 0;
}


//============================================================================
void GetNumFiles(int *NumFiles, const char **dn, int nu){
  int cnt;
  DIR *dp;
  struct dirent *entry;

  for(int i=0;i<nu;i++){
    if((dp = opendir(dn[i]))==NULL){
      fprintf(stderr,"Cannot open %s\n",dn[i]);
      exit(1);
    }
    cnt=0;
    while((entry=readdir(dp))!=NULL){
      if(entry->d_type == DT_REG) cnt++;
    } // while
    NumFiles[i]=cnt;
    closedir(dp);
  } // for i
}

//============================================================================
void GetFileNames(char ***fn, const char **dn, int nu){
  int cnt;
  DIR *dp;
  struct dirent *entry;

  for(int i=0;i<nu;i++){
    if((dp = opendir(dn[i]))==NULL){
      fprintf(stderr,"Cannot open %s\n",dn[i]);
      exit(1);
    }
    cnt=0;
    while((entry=readdir(dp))!=NULL){
      if(entry->d_type == DT_REG){
	sprintf(fn[i][cnt],"%s/%s",dn[i],entry->d_name);
	cnt++;
      }
    } // while
    closedir(dp);
  } // for i
}


//============================================================================
void Generate(const char ***fn, const int *sel, int nu, const char *odir, int cnt){
  int  beta=0;
  int  *B   = new int [nu];
  int  *B2p = new int [nu];
  int  **CW = new int * [nu];
  char *ofn = new char [BSIZE];
  class BaseFile *BFin;

  // read nu files
  for(int idx=0;idx<nu;idx++){
    BFin = new class BaseFile(fn[idx][sel[idx]]);
    if(idx==0) beta = BFin->get_beta();
    else       assert(BFin->get_beta()==beta);
    assert(BFin->get_nu()==1);
    BFin->get_B(&B[idx]);
    BFin->get_B2p(&B2p[idx]);
    CW[idx] = new int [B2p[idx]];
    BFin->get_CW(0,CW[idx]);
    
    delete BFin;
  } // for idx

  // write
  sprintf(ofn,"%s/%08d.txt",odir,cnt);
  class BaseFile *BFout = new class BaseFile(beta,nu,B);
  for(int i=0;i<nu;i++) BFout->set_CW(i,CW[i]);
  BFout->write(ofn);
  delete BFout;

  for(int i=0;i<nu;i++) printf("%s ",fn[i][sel[i]]);
  PrintVect(B,  nu,"(B) "," ");
  PrintVect(B2p,nu,"(B2p) "," : ");
  for(int i=0;i<nu;i++) PrintVect(CW[i],B2p[i]," ","\t");
  printf("%s\n",ofn);

  for(int i=0;i<nu;i++) delete [] CW[i];
  delete [] CW;
  delete [] B;
  delete [] B2p;
  delete [] ofn;
}

//============================================================================
void PrintFileNames(const char ***fn, const int *NumFiles, int nu){
  for(int i=0;i<nu;i++){
    printf("# [dir %d]\n",i);
    for(int j=0;j<NumFiles[i];j++) printf("#  %s\n",fn[i][j]);
  } // for i
}
