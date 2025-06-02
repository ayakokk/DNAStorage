#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"
#include "StatICB.hpp"


//=================================================================================
void StatICB::InitHdCnt(){
  for(int i=0;i<nu;i++)
    for(int x=0;x<beta4p;x++)
      for(int y=0;y<(beta4p*lmd);y++)
	  HdCnt[i][x][y]=0;
}

//=================================================================================
double StatICB::CalcIxy(const unsigned long **C, int numX, int numY){
  assert(numX>0 && numY>0);
  double *Px = new double [numX];
  double *Py = new double [numY];
  double *Q  = new double [numX];
  double *Hxy= new double [numY];
  double hx,Ixy;
  unsigned long *Cx = new unsigned long [numX];
  unsigned long *Cy = new unsigned long [numY];
  unsigned long sum=0;

  ClearVect(Cx,numX);
  ClearVect(Cy,numY);
  
  // count
  for(int i=0;i<numX;i++){
    for(int j=0;j<numY;j++){
     sum+=C[i][j];
     Cx[i]+=C[i][j];
     Cy[j]+=C[i][j];
    } // for i
  } // for j
  assert(sum>0);
  
  //Px,Py
  for(int i=0;i<numX;i++) Px[i] = (double)Cx[i]/sum;
  for(int j=0;j<numY;j++) Py[j] = (double)Cy[j]/sum;
  hx = Entropy(Px,numX);
  //hy = Entropy(Py,numY);

  // H(X|Y=y)
  for(int j=0;j<numY;j++){
    if(Cy[j]>0){
      for(int i=0;i<numX;i++) Q[i] = (double)C[i][j]/Cy[j];
      Hxy[j] = Entropy(Q,numX);
    } else {
      Hxy[j] = 0.0; 
    } // if Cy[j]
  } // for j

  // I(X;Y)
  Ixy=0.0;
  for(int j=0;j<numY;j++) Ixy+= Py[j]*Hxy[j];
  Ixy = hx-Ixy;
  
  //(dbg)
  // PrintVect(Px, numX,"Px: ","\n");
  // PrintVect(Py, numY,"Py: ","\n");
  // PrintVect(Hxy,numY,"Hxy:","\n");
  // printf("H(X)=%e H(Y)=%e Ixy=%e\n",hx,hy,Ixy);
  
  delete [] Px;
  delete [] Py;
  delete [] Q;
  delete [] Hxy;
  delete [] Cx;
  delete [] Cy;
  return Ixy;
}

//=================================================================================
void StatICB::SetY(int *Y, const double **P, int Nseg){
  for(int i=0;i<Nseg;i++){
    Y[i] = ArgMax(P[i],beta4p);
  } // for i
}

//=================================================================================
void StatICB::SetHp(int *Hp,const double **P, int Nseg){
  int x;
  for(int i=0;i<Nseg;i++){
    x = (int)floor( (Entropy(P[i],beta4p)*(double)lmd)/B[i%nu] );
    if(x>=lmd) x=lmd-1;
    if(x<0)    x=0;
    Hp[i] = x;
  } // for i
}

//=================================================================================
double StatICB::Entropy(const double *P, int len){
  double h=0.0;
  for(int i=0;i<len;i++){
    assert(P[i]>=0);
    if(P[i]>0) h -= P[i]*log2(P[i]);
  } // for i
  return h;
}

//=================================================================================
int StatICB::ArgMax(const double *V, int len){
  int i,r=0;
  for(i=1;i<len;i++)
    if(V[i]>V[r]) r=i;
  return r;
}

//=================================================================================
void StatICB::ClearVect(unsigned long *V, int len){
  for(int i=0;i<len;i++) V[i]=0;
}

//=================================================================================
//=================================================================================
//=================================================================================

//=================================================================================
StatICB::StatICB(int _nu, int _beta, const int *_B, int _lmd){
  nu   = _nu;
  beta = _beta;
  lmd  = _lmd;
  B    = new int [nu];
  memcpy(B,_B,sizeof(int)*nu);
  Bsum = Sum(B,nu);
  
  beta4p  = (int)pow(4,beta);
  printf("# StatICB: nu=%d beta=%d beta4p=%d Bsum=%d lmd=%d\n",nu,beta,beta4p,Bsum,lmd);
  PrintVect(B,nu,"# StatICB: B=","\n");
  assert(nu>0);
  assert(beta4p>0);

  // HdCnt
  HdCnt = new unsigned long ** [nu];
  for(int i=0;i<nu;i++){
    HdCnt[i] = new unsigned long * [beta4p];
    for(int j=0;j<beta4p;j++) HdCnt[i][j] = new unsigned long [beta4p*lmd];
  } // for i
  InitHdCnt();
}

//=================================================================================
StatICB::~StatICB(){
  for(int i=0;i<nu;i++){
    for(int j=0;j<beta4p;j++) delete [] HdCnt[i][j];
    delete [] HdCnt[i];
  } // for i
  delete [] HdCnt;
  delete [] B;
  printf("# StatICB: deleted\n");
}

//=================================================================================
void StatICB::count(const int *X, const double **P, int Nseg){
  int *Y = new int [Nseg];
  int *Hp= new int [Nseg];
  SetY( Y, P,Nseg);
  SetHp(Hp,P,Nseg);
  for(int i=0;i<Nseg;i++){
      assert(X[i] >=0 && X[i] <beta4p);
      assert(Y[i] >=0 && Y[i] <beta4p);
      assert(Hp[i]>=0 && Hp[i]<lmd);
      HdCnt[i%nu][X[i]][Hp[i]*beta4p+Y[i]]++;
  } // for i

  //(dbg)
  /*
  for(int i=0;i<Nseg;i++)
    printf("%04d: %03d %03d %03d\n",i,X[i],Y[i],Hp[i]);
  */
  
  delete [] Y;
  delete [] Hp;
}

//=================================================================================
void StatICB::PrintIxy(){
  double v, sum=0;
  
  for(int i=0;i<nu;i++){
    v = CalcIxy((const unsigned long **)HdCnt[i],beta4p,beta4p*lmd);
    sum += v;
    printf("%e ",v);
  } // for j
  printf(": %e\n",sum/Bsum);

  //CalcIxy((const unsigned long**)HdCnt[0][0],beta4p,beta4p+1);//TMP
}

//=================================================================================
void StatICB::dump(){
  for(int i=0;i<nu;i++){
    printf("[%d]\n",i);
    for(int x=0;x<beta4p;x++){
      printf("%03d: ",x);
      for(int y=0;y<beta4p*lmd;y++){
	printf("%04lu ",HdCnt[i][x][y]);
	if(y%beta4p==beta4p-1) printf("| ");
      } // for y
      printf("\n");
    } // for x
  } // for i
}
