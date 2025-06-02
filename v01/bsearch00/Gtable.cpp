#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include <string>
#include <cstring>
#include <fstream>
#include <iostream>

#include "bmatrix.hpp"
#include "func.hpp"
#include "Gtable.hpp"


//=================================================================================
double Gtable::Psub(unsigned char x, unsigned char y){
  assert(x<4 && y<4);
  if(x==y) return 1.0-Ps;
  else     return Ps/3.0;
}

//=================================================================================
double Gtable::CalcProb(const unsigned char *X, int aX, int bX,
			const unsigned char *Y, int aY, int bY){
  assert(aX>=0 && aX< bX);
  assert(aY>=0 && aY<=bY);
  int d0 = aY-aX;
  int d1 = bY-bX;
  double p=0.0;
  //printf("aX=%d bX=%d aY=%d bY=%d\n",aX,bX,aY,bY);
  
  if(aX==bX-1){ // 1 symbol of X
    if(     d1==d0   && aY  <bY) p=Ptrn* Psub(X[aX],Y[aY]);                       // trans
    else if(d1==d0+1 && aY+1<bY) p=Pid * Psub(X[aX],Y[aY]) * Psub(X[aX],Y[aY+1]); // ins
    else if(d1==d0-1           ) p=Pid;                                           // del
  } else {      // recursion
    if(aY<=bY-1) p+= (Ptrn * Psub(X[aX],Y[aY]) *                       CalcProb(X,aX+1,bX,Y,aY+1,bY)); // trans
    if(aY<=bY-2) p+= (Pid  * Psub(X[aX],Y[aY]) * Psub(X[aX],Y[aY+1]) * CalcProb(X,aX+1,bX,Y,aY+2,bY)); // ins
    p+=              (Pid  *                                           CalcProb(X,aX+1,bX,Y,aY  ,bY)); // del
  } //if alpha==beta-1
  
  //printf("aX=%d bX=%d aY=%d bY=%d: %e\n",aX,bX,aY,bY,p);
  return p;
}

//=================================================================================
//=================================================================================
//=================================================================================

//=================================================================================
Gtable::Gtable(int _beta, double _Pid, double _Ps, double _Pth){
  beta = _beta;
  Pid  = _Pid;
  Ps   = _Ps;
  Ptrn = 1.0 - Pid*2.0;
  Pth  = _Pth;
  printf("# Gtable: beta=%d (Pid,Ps,Ptrn)=(%e,%e,%e) Pth=%e\n",beta,Pid,Ps,Ptrn,Pth);
  assert(beta*2*2 < (int)sizeof(int)*8); // bit width
  beta4p = (int)pow(4,beta);
  beta4pp= (int)pow(4,beta*2);
  printf("# Gtable: beta4p=%d beta4pp=%d\n",beta4p,beta4pp);
    
  // LIST -----
  LIST = new std::list<Gtable_xp>* [beta*2+1];
  for(int i=0;i<=2*beta;i++) LIST[i] = new std::list<Gtable_xp> [(int)pow(4,i)];
  SetLIST();
  
  //dump();
}

//=================================================================================
Gtable::Gtable(const char *fn){
  struct Gtable_xp xp;
  int   x,num,cnt;
  float p;
  std::ifstream fin;
  fin.open(fn, std::ios::in|std::ios::binary);
  if(!fin){
    std::cout << "Cannot open " << fn << std::endl;
    exit(1);
  }
  fin.read((char *)&beta,sizeof(int));
  fin.read((char *)&Pid, sizeof(double));
  fin.read((char *)&Ps,  sizeof(double));
  fin.read((char *)&Pth, sizeof(double));
  Ptrn = 1.0 - Pid*2.0;
  printf("# Gtable: beta=%d (Pid,Ps,Ptrn)=(%e,%e,%e) Pth=%e\n",beta,Pid,Ps,Ptrn,Pth);
  beta4p = (int)pow(4,beta);
  beta4pp= (int)pow(4,beta*2);
  printf("# Gtable: beta4p=%d beta4pp=%d\n",beta4p,beta4pp);

  // LIST -----
  LIST = new std::list<Gtable_xp>* [beta*2+1];
  for(int i=0;i<=2*beta;i++) LIST[i] = new std::list<Gtable_xp> [(int)pow(4,i)];
  for(int len=0;len<=beta*2;len++){
    cnt=0;
    for(int y=0;y<(int)pow(4,len);y++){
      fin.read((char *)&num,sizeof(int));
      for(int j=0;j<num;j++){
	fin.read((char *)&x,sizeof(int));
	fin.read((char *)&p,sizeof(float));
	xp.x = x;
	xp.p = p;
	LIST[len][y].push_back(xp);
	cnt++;
      } // for j
      /*
      num = LIST[len][y].size();
      fout.write((char *)&num,sizeof(int));
      for(auto XP=LIST[len][y].begin();XP!=LIST[len][y].end();++XP){
	x = (*XP).x;
	p = (*XP).p;
	fout.write((char *)&x,sizeof(int));
	fout.write((char *)&p,sizeof(float));
      } // for XP
      */
    } // for y
    printf("# Gtable::Gtable len=%d cnt=%d\n",len,cnt);
  } // for len

  fin.close();
}

//=================================================================================
Gtable::~Gtable(){
  for(int i=0;i<beta*2+1;i++) delete [] LIST[i];
  delete [] LIST;
  printf("# Gtable: deleted\n");
}

//=================================================================================
void Gtable::SetLIST(){
  int    cnt;
  double p;
  unsigned char *X = new unsigned char [beta];
  unsigned char *Y = new unsigned char [beta*2];
  struct Gtable_xp xp;
  
  // ConvIntVect(100,X,beta);
  // ConvIntVect(100,Y,beta*2);
  // CalcProb(X,0,beta, Y,0,beta*2);

  for(int beta2=0;beta2<=beta*2;beta2++){
    cnt=0;
    for(int y=0;y<(int)pow(4,beta2);y++){
      ConvIntVect(y,Y,beta2);
      //PrintVect4(Y,beta2); printf("\n");
      for(int x=0;x<beta4p;x++){
	ConvIntVect(x,X,beta);
	p = CalcProb(X,0,beta, Y,0,beta2);
	if(p>=Pth){
	  xp.x = x;
	  xp.p = p;
	  LIST[beta2][y].push_back(xp);
	  cnt++;
	  //printf(" >"); PrintVect4(X,beta); printf(": %e\n",p);
	}
	//if(p>=1.0e-5){printf(" >"); PrintVect4(X,beta); printf(": %e\n",p);}
      } // for x
      LIST[beta2][y].sort();
    } // for y
    printf("# Gtable::SetLIST beta2=%d cnt=%d\n",beta2,cnt);
  } // for beta2
  
  delete [] X;
  delete [] Y;
}

//=================================================================================
void Gtable::write(const char *fn){
  int   num,x;
  float p;
  std::ofstream fout;
  fout.open(fn, std::ios::out|std::ios::binary|std::ios::trunc);
  if(!fout){
    std::cout << "Cannot open " << fn << std::endl;
    exit(1);
  }
  fout.write((char *)&beta,sizeof(int));
  fout.write((char *)&Pid, sizeof(double));
  fout.write((char *)&Ps,  sizeof(double));
  fout.write((char *)&Pth, sizeof(double));
  for(int len=0;len<=beta*2;len++){
    for(int y=0;y<(int)pow(4,len);y++){
      num = LIST[len][y].size();
      fout.write((char *)&num,sizeof(int));
      for(auto XP=LIST[len][y].begin();XP!=LIST[len][y].end();++XP){
	x = (*XP).x;
	p = (*XP).p;
	fout.write((char *)&x,sizeof(int));
	fout.write((char *)&p,sizeof(float));
      } // for XP
    } // for y
  } // for len
  fout.close();
}

//=================================================================================
int    Gtable::getBeta(){return beta;};
double Gtable::getPid(){ return Pid;};
double Gtable::getPs(){  return Ps;};
double Gtable::getPth(){ return Pth;};

//=================================================================================
void Gtable::dump(){
  unsigned char *X = new unsigned char [beta];
  unsigned char *Y = new unsigned char [beta*2];

  // Param
  printf("beta=%d (Pid,Ps)=(%e,%e) Pth=%e\n",beta,Pid,Ps,Pth);
  
  // LIST
  for(int len=0;len<=beta*2;len++){
    printf("[%d]\n",len);
    for(int y=0;y<(int)pow(4,len);y++){
      ConvIntVect(y,Y,len);
      printf(" %04X ",y);
      PrintVect4(Y,len);
      printf("[%03lu] ",LIST[len][y].size());
      //LIST[len][y].sort();
      for(auto XP=LIST[len][y].begin();XP!=LIST[len][y].end();++XP){
	ConvIntVect((*XP).x,X,beta);
	PrintVect4(X,beta);
	printf(" %.2e ",(*XP).p);
	//printf("(%04X)%.2e ",(*XP).x,(*XP).p);
      } // for XP
      printf("\n");
    } // for y
  } // for len
  delete [] X;
  delete [] Y;
}
