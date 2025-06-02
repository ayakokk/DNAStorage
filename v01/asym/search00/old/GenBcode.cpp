#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "bmatrix.hpp"
#include "func.hpp"

void GenTR(int *TR0, int *TR1, int lmd, int Lr, int Lb, int wL, int wH);
void RemoveInvalidPath(int *TR0, int *TR1, int lmd);
bool CheckTR(const int *TR0, const int *TR1, int lmd);
void PrintTR(const int *TR0, const int *TR1, int lmd);
void GenCB(class bmatrix *CB, const int *TR0, const int *TR1, int lmd, int beta);
void GenCBrec(unsigned char *CW, const int *TR0, const int *TR1,
	      class bmatrix *CB, int src, int stat, int pos, int beta);
void PrintCB(class bmatrix *CB);

//================================================================================
int main(int argc, char *argv[]){
  int    beta,Lr,Lb,lmd,seed;
  int    wL,wH;     // LGCB min max
  double eps;
  unsigned int l2p; // 2^lmd
  unsigned int b2p; // 2^beta
  char   *fn;
  if(argc!=7){
    fprintf(stderr,"Usage: %s <beta> <Lr> <Lb> <eps> <Bcode.bin> <seed|-1>\n",argv[0]);
    return 1;
  }
  beta = atoi(argv[1]);
  Lr   = atoi(argv[2]);
  Lb   = atoi(argv[3]);
  eps  = atof(argv[4]);
  fn   =      argv[5];
  seed = atoi(argv[6]);
  if(seed==-1) seed = (int)time(NULL);
  srandom(seed);
  lmd  = max(Lr,Lb)-1;
  l2p  = (unsigned int)pow(2,lmd);
  b2p  = (unsigned int)pow(2,beta);
  wL   = (int)ceil( Lb*(0.5-eps));
  wH   = (int)floor(Lb*(0.5+eps));
  printf("# CodeParam: beta=%d b2p=%d Lr=%d (Lb,eps)=(%d,%e) [%d]\n",beta,b2p,Lr,Lb,eps,seed);
  printf("# CodeParam: lmd=%d l2p=%u (wL,wH)=(%d,%d)\n",lmd,l2p,wL,wH);
  printf("# Output: %s\n",fn);
  assert(beta>0);
  assert(Lr>0 && Lb>0);
  assert(eps>0 && eps<=0.5);

  // generate trellis
  int *TR0 = new int [l2p];
  int *TR1 = new int [l2p];
  GenTR(TR0,TR1,lmd,Lr,Lb,wL,wH);
  //PrintTR(TR0,TR1,lmd); //(dbg)

  // remove invalid paths
  RemoveInvalidPath(TR0,TR1,lmd);
  assert(CheckTR(TR0,TR1,lmd));
  //PrintTR(TR0,TR1,lmd); //(dbg)

  // generate codebook
  class bmatrix *CB = new class bmatrix(l2p,b2p); // [stat][cw_flg]
  GenCB(CB,TR0,TR1,lmd,beta);

  // write
  CB->write(fn);
  PrintCB(CB);
  
  // delete
  delete [] TR0;
  delete [] TR1;
  delete CB;
  return 0;
} 

//================================================================================
void GenTR(int *TR0, int *TR1, int lmd, int Lr, int Lb, int wL, int wH){
  int  RL,HW;
  bool flg;
  unsigned int l2p = (unsigned int)pow(2,lmd);
  unsigned char *S  = new unsigned char [lmd+1];
  // init
  for(unsigned int s=0;s<l2p;s++){
    TR0[s]=-1;
    TR1[s]=-1;
  } // for s
  // set 
  for(unsigned int s=0;s<l2p*2;s++){
    ConvIntV2(S,s,lmd+1);
    RL = MaxRunLength( S,lmd+1);
    HW = HammingWeight(S,lmd+1);
    flg= (RL<=Lr) && (HW>=wL) && (HW<=wH);  // satisfy constraints
    if(flg){
      if(s%2==0) TR0[s/2] = s%l2p;
      else       TR1[s/2] = s%l2p;
    } // if flg
    //(dbg)
    //PrintVect2(S,lmd+1,"",":");
    //printf("%02d %02d %d\n",RL,HW,flg);
  } // for s
  delete [] S;
}

//================================================================================
void RemoveInvalidPath(int *TR0, int *TR1, int lmd){
  unsigned int l2p = (unsigned int)pow(2,lmd);
  bool upf;
  bool *vf = new bool [l2p];
  SetVal(vf,true,l2p);

  while(1){
    upf = false;
    // set vf
    for(unsigned int s=0;s<l2p;s++){
      if(TR0[s]<0 && TR1[s]<0 && vf[s]){
	vf[s]=false;
	upf = true;
      }
    } // for s
    // set TR
    for(unsigned int s=0;s<l2p;s++){
      if(TR0[s]>=0){
	if(!vf[TR0[s]]){
	  TR0[s]=-1;
	  upf = true;
	  printf("  del TR0: %u\n",s);
	} // if !vf
      } // if TR0
      if(TR1[s]>=0){
	if(!vf[TR1[s]]){
	  TR1[s]=-1;
	  upf = true;
	  printf("  del TR1: %u\n",s);
	} // if !vf
      } // if TR0    
    } // for s
    if(!upf) break;
  } // while(1)
  
  delete [] vf;
}

//================================================================================
bool CheckTR(const int *TR0, const int *TR1, int lmd){
  unsigned int l2p = (unsigned int)pow(2,lmd);
  int s0,s1;
  for(unsigned int s=0;s<l2p;s++){
    s0 = TR0[s];
    s1 = TR1[s];
    if(s0>=0){
      if(TR0[s0]<0 && TR1[s0]<0) return false;
    }
    if(s1>=0){
      if(TR0[s1]<0 && TR1[s1]<0) return false;
    }    
  } // for s
  return true;
}

//================================================================================
void PrintTR(const int *TR0, const int *TR1, int lmd){
  unsigned int l2p = (unsigned int)pow(2,lmd);
  unsigned char *S = new unsigned char [lmd];
  
  for(unsigned int s=0;s<l2p;s++){
    printf("%04d: %04d %04d: ",s,TR0[s],TR1[s]);
    ConvIntV2(S,s,lmd);
    PrintVect2(S,lmd,""," ");
    if(TR0[s]>=0){
      ConvIntV2(S,(unsigned int)TR0[s],lmd);
      PrintVect2(S,lmd,""," ");
    } else {
      PrintSP(lmd+1);
    }
    if(TR1[s]>=0){
      ConvIntV2(S,(unsigned int)TR1[s],lmd);
      PrintVect2(S,lmd,""," ");
    } else {
      PrintSP(lmd+1);
    }    
    printf("\n");
  } // s
  
  delete [] S;
}

//================================================================================
void GenCB(class bmatrix *CB, const int *TR0, const int *TR1, int lmd, int beta){
  unsigned int l2p = (unsigned int)pow(2,lmd);
  unsigned int b2p = (unsigned int)pow(2,beta);
  unsigned int s;
  unsigned char *CW = new unsigned char [beta];
  assert(CB->getM()==(int)l2p);
  assert(CB->getN()==(int)b2p);
  CB->clear();
  for(s=0;s<l2p;s++){
    GenCBrec(CW,TR0,TR1,CB,s,s,beta-1,beta);
  } // for s
  delete [] CW;
}

//================================================================================
void GenCBrec(unsigned char *CW, const int *TR0, const int *TR1,
	      class bmatrix *CB, int src, int stat, int pos, int beta){
  //printf("[pos=%d]\n",pos);
  if(pos==beta-1){
    // --------------- pos==beta-1
    assert(src==stat);
    if(TR0[stat]>=0){
      if(stat%2==1 || TR1[stat]<0){
	CW[pos]=0;
	GenCBrec(CW,TR0,TR1,CB,src,TR0[stat],pos-1,beta);
      } // if stat%2
    } // if TR0
    if(TR1[stat]>=0){
      if(stat%2==0 || TR0[stat]<0){
	CW[pos]=1;
	GenCBrec(CW,TR0,TR1,CB,src,TR1[stat],pos-1,beta);
      } // if stat%2
    } // if TR1
  } else if(pos>0){
    // --------------- pos==1...beta-2
    if(TR0[stat]>=0){
      CW[pos]=0;
      GenCBrec(CW,TR0,TR1,CB,src,TR0[stat],pos-1,beta);
    } // if TR0
    if(TR1[stat]>=0){
      CW[pos]=1;
      GenCBrec(CW,TR0,TR1,CB,src,TR1[stat],pos-1,beta);
    } // if TR1    
  } else {  
    // --------------- pos==0
    assert(pos==0);
    unsigned int c;
    if(TR0[stat]>=0){
      CW[pos]=0;
      c = ConvV2Int(CW,beta);
      CB->setV(src,c,1);
      //PrintVect2(CW,beta,"","\n");
    } // if TR0
    if(TR1[stat]>=0){
      CW[pos]=1;
      c = ConvV2Int(CW,beta);
      CB->setV(src,c,1);
      //PrintVect2(CW,beta,"","\n");
    } // if TR1
  } // if pos
}

//================================================================================
void PrintCB(class bmatrix *CB){
  unsigned int l2p = (unsigned int)CB->getM();
  unsigned int b2p = (unsigned int)CB->getN();
  unsigned int s,v;
  int lmd = (int)log2(l2p);
  int beta= (int)log2(b2p);
  unsigned char *S = new unsigned char [lmd];
  unsigned char *V = new unsigned char [beta];
  
  printf("l2p=%u b2p=%u : lmd=%d beta=%d\n",l2p,b2p,lmd,beta);
  assert(l2p==(unsigned int)pow(2,lmd));
  assert(b2p==(unsigned int)pow(2,beta));
  for(s=0;s<l2p;s++){
    ConvIntV2(S,s,lmd);
    printf("%04X:",s);
    PrintVect2(S,lmd,"",": ");
    for(v=0;v<b2p;v++){
      if(CB->getV(s,v)==1){
	ConvIntV2(V,v,beta);
	PrintVect2(V,beta,""," ");
      } // if
    } // for v
    printf("\n");
  } // for s

  delete [] S;
  delete [] V;
}
