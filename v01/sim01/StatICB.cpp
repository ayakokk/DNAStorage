#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "func.hpp"
#include "StatICB.hpp"

#define _HdThNum  20    // >=1: num of thresholds / 0: no erasure 
#define _HdThStep 0.05 


//=================================================================================
void StatICB::InitHdCnt(){
  for(int i=0;i<HdThNum;i++)
    for(int j=0;j<nu;j++)
      for(int x=0;x<beta4p;x++)
	for(int y=0;y<beta4p+1;y++)
	  HdCnt[i][j][x][y]=0;
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
void StatICB::SetY(int *Y, const double **P, double Pth, int Nseg){
  for(int i=0;i<Nseg;i++){
    Y[i] = ArgMax(P[i],beta4p);
    if(P[i][Y[i]]<Pth) Y[i]=beta4p; // erasure
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
StatICB::StatICB(int _nu, int _beta, const int *_B){
  nu   = _nu;
  beta = _beta;
  B    = new int [nu];
  memcpy(B,_B,sizeof(int)*nu);
  Bsum=Sum(B,nu);
  
  beta4p  = (int)pow(4,beta);
  HdThNum = max(_HdThNum,1);
  HdThStep= _HdThStep;
  printf("# StatICB: nu=%d beta=%d beta4p=%d Bsum=%d HdThNum=%d HdThStep=%.2e\n",
	 nu,beta,beta4p,Bsum,HdThNum,HdThStep);
  PrintVect(B,nu,"# StatICB: B=","\n");
  assert(nu>0);
  assert(beta4p>0);

  // hard decision threshold list
  HdThList = new double [HdThNum];
  if(_HdThNum<=0){
    HdThList[0]=1.0e-10;
  } else {
    for(int i=0;i<HdThNum;i++){
      HdThList[i] = 1.0-((double)(i+1)*HdThStep);
      if(HdThList[i]<=0) HdThList[i]=1.0e-10;
    } // for i
  }
  PrintVect(HdThList,HdThNum,"# StatICB: HdThList ","\n");
  //assert(HdThList[HdThNum-1]>0);

  // HdCnt
  HdCnt = new unsigned long *** [HdThNum];
  for(int i=0;i<HdThNum;i++){
    HdCnt[i] = new unsigned long ** [nu];
    for(int j=0;j<nu;j++){
      HdCnt[i][j] = new unsigned long * [beta4p];
      for(int k=0;k<beta4p;k++) HdCnt[i][j][k] = new unsigned long [beta4p+1];
    } // for j
  } // for i
  InitHdCnt();


}

//=================================================================================
StatICB::~StatICB(){
  delete [] HdThList;
  for(int i=0;i<HdThNum;i++){
    for(int j=0;j<nu;j++){
      for(int k=0;k<beta4p;k++) delete [] HdCnt[i][j][k];
      delete [] HdCnt[i][j];
    } // for j
    delete [] HdCnt[i];
  } // for i
  delete [] HdCnt;
  delete [] B;
  printf("# StatICB: deleted\n");
}

//=================================================================================
void StatICB::count(const int *X, const double **P, int Nseg){
  int *Y = new int [Nseg];
  for(int ith=0;ith<HdThNum;ith++){
    SetY(Y,P,HdThList[ith],Nseg);
    for(int i=0;i<Nseg;i++){
      assert(X[i]>=0 && X[i]<beta4p);
      assert(Y[i]>=0 && Y[i]<beta4p+1);
      HdCnt[ith][i%nu][X[i]][Y[i]]++;
    } // for i
    
    //(dbg)
    // printf("[ith=%d HdTh=%.2e]\n",ith,HdThList[ith]);
    // for(int i=0;i<Nseg;i++) printf("%03d: %03d %03d\n",i,X[i],Y[i]);
  } // for ith
  delete [] Y;
}

//=================================================================================
void StatICB::PrintIxy(){
  double v, sum;
  
  for(int i=0;i<HdThNum;i++){
    printf("[%02d] Th=%e: ",i,HdThList[i]);
    sum = 0;
    for(int j=0;j<nu;j++){
      v = CalcIxy((const unsigned long **)HdCnt[i][j],beta4p,beta4p+1);
      sum += v;
      printf("%e ",v);
    } // for j
    printf(": %e\n",sum/Bsum);
  } // for ith
  
  //CalcIxy((const unsigned long**)HdCnt[0][0],beta4p,beta4p+1);//TMP
}

//=================================================================================
void StatICB::dump(){
  for(int i=0;i<HdThNum;i++){
    for(int j=0;j<nu;j++){
      printf("Pers=%.2e (%d)\n",HdThList[i],j);
      for(int x=0;x<beta4p;x++){
	printf("%03d: ",x);
	for(int y=0;y<beta4p+1;y++) printf("%04lu ",HdCnt[i][j][x][y]);
	printf("\n");
      } // for x
    } // for j
  } // for i
}
