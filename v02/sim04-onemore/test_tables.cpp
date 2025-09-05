#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"
#include "ChannelMatrix.hpp"
#include "InnerCodebook.hpp"
#include "IDSchannel.hpp"
#include "SLFBAdec.hpp"

#define BSIZE 8192

int main(int argc, char *argv[]){
  int Rho;         // run-length    [constraint.txt]
  int ell,Delta;   // local-balance [constraint.txt]
  int N;           // block length (symbols)
  int Nb;      // block length & recv length (bits)
  int Q, Nu;       // numCW, symbol-len [ICB]
  double Pi,Pd,Ps; // IDS prob
  char *fn;
  char *fncb    = new char [BSIZE];
  char *fnconst = new char [BSIZE];
  char *fncm    = new char [BSIZE];
  
  if(argc!=6){
    fprintf(stderr,"Usage: %s <ICB_dir> <N> <Pi> <Pd> <Ps>\\n",argv[0]);
    fprintf(stderr,"Example: %s ICB/ex01/ 3 0.1 0.1 0.1\\n",argv[0]);
    return 1;
  }
  
  fn   =      argv[1];
  N    = atoi(argv[2]);
  Pi   = atof(argv[3]);
  Pd   = atof(argv[4]);
  Ps   = atof(argv[5]);
  
  assert(N>0);
  assert(Pi>=0.0 && Pi<0.5);
  assert(Pd>=0.0 && Pd<0.5);
  assert(Ps>=0.0 && Ps<0.5);
  
  snprintf(fncb,   BSIZE,"%s/cb.txt",        fn);  // inner codebook (in)
  snprintf(fnconst,BSIZE,"%s/constraint.txt",fn);  // constraints (in)
  snprintf(fncm,   BSIZE,"%s/EncCM.bin",     fn);  // encoding channel matrix (in)
  
  printf("# SLFBAdec sim01 Table Export Test\\n");
  printf("# ICB:   %s\\n",fncb);
  printf("# Const: %s\\n",fnconst);
  printf("# EncCM: %s\\n",fncm);
  printf("# Parameters: N=%d Pi=%.6f Pd=%.6f Ps=%.6f\\n",N,Pi,Pd,Ps);
  
  ReadConstraints(fnconst, &Rho, &ell, &Delta);
  class InnerCodebook *ICB = new class InnerCodebook(fncb,Rho,ell,Delta);
  Q    = ICB->Get_numCW();
  Nu   = ICB->Get_Nu();
  Nb   = N*Nu;
  printf("# Q=%d Nu=%d Nb=%d\\n",Q,Nu,Nb);
  
  class ChannelMatrix *ECM = new class ChannelMatrix(fncm);
  class IDSchannel    *CH  = new class IDSchannel(Nb,Pi,Pd,Ps);
  class SLFBAdec      *DEC = new class SLFBAdec(ICB,ECM,CH);
  
  // テーブル出力実行
  printf("\\n# Exporting probability tables...\\n");
  DEC->exportAllTables("sim01_tables");
  
  printf("\\n# Table export completed!\\n");
  printf("# Generated files:\\n");
  printf("#   sim01_tables/GD_table.txt - Drift probability table\\n");
  printf("#   sim01_tables/GX_table.txt - Channel output probability table\\n");
  printf("#   sim01_tables/sim01_tables_summary.txt - Summary information\\n");
  
  // クリーンアップ
  delete ICB;
  delete ECM;
  delete CH;
  delete DEC;
  delete [] fncb;
  delete [] fnconst;
  delete [] fncm;
  
  return 0;
}