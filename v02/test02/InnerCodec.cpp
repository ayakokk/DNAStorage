#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "InnerCodec.hpp"
#include "InnerCodecFenc.cpp"

//================================================================================
int InnerCodec::BitGet(int val, int pos){
  assert(pos>=0 && pos<(int)sizeof(int)*8-1);
  int mask = (0x1 << pos);
  return ((val & mask)==0)? 0 : 1;
}

//================================================================================
void InnerCodec::BitSet(int *val, int pos, int sval){
  assert(pos>=0 && pos<(int)sizeof(int)*8-1);
  assert(sval==0 || sval==1);
  int mask = (0x1 << pos);
  if(sval==0) (*val) = (*val) & (~mask);  // set 0
  else        (*val) = (*val) | mask;     // set 1
}

//================================================================================
void InnerCodec::Normalize(double *V, int len){
  assert(len>0);
  double s = 0;
  for(int i=0;i<len;i++){
    assert(V[i]>=0.0);
    s += V[i];
  } // for i
  assert( s>0 );
  for(int i=0;i<len;i++) V[i]/=s;
}

//================================================================================
void InnerCodec::Normalize(double ***V, int len0, int len1, int len2){
  assert(len0>0 && len1>0 && len2>0);
  double s = 0;
  for(int i0=0;i0<len0;i0++){
    for(int i1=0;i1<len1;i1++){
      for(int i2=0;i2<len2;i2++){
	assert(V[i0][i1][i2]>=0.0);
	s += V[i0][i1][i2];
      } // for i2
    } // for i1
  } // for i0
  assert( s>0 );
  for(int i0=0;i0<len0;i0++){
    for(int i1=0;i1<len1;i1++){
      for(int i2=0;i2<len2;i2++){
	V[i0][i1][i2] /= s;
      } // for i2
    } // for i1
  } // for i0
}

//================================================================================
void InnerCodec::ClearVect(double *V, int len){
  for(int i=0;i<len;i++) V[i]=0.0;
}

//================================================================================
void InnerCodec::PrintVect(const int *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d ",V[i]);
  printf("%s",post);
}

//================================================================================
void InnerCodec::PrintVect(const double *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%.2e ",V[i]);
  printf("%s",post);
}

//================================================================================
double InnerCodec::SumVect(const double *V, int len){
  double s=0.0;
  for(int i=0;i<len;i++) s+=V[i];
  return s;
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
int InnerCodec::CalcL(int p, int q){
  assert(p>=0 && p<NumSTR);
  assert(q>=0 && q<NumSTB);
  int L = 0;
  for(int a=0;a<4;a++){
    if(STR[p][a]>=0 && STB[q][a]>=0) BitSet(&L,a,1);
  } // for a
  assert(L>=0 && L<=15);
  return L;
}

//================================================================================
double InnerCodec::FuncPd(int d0, int d1){
  // assert(d0>=Dmin && d0<=Dmax);
  // assert(d1>=Dmin && d1<=Dmax);
  if( (d0<Dmin) || (d0>Dmax) || (d1<Dmin) || (d1>Dmax) ) return 0;
  if(abs(d0-d1)>=2) return 0;
  if(d0==Dmin){
    if(d1==d0)        return 1.0-Pi;
    else if(d1==d0+1) return Pi;
    else              return 0;
  } else if(d0==Dmax) {
    if(d1==d0)        return 1.0-Pd;
    else if(d1==d0+1) return 0;
    else              return Pd;
  } else {
    if(d1==d0)        return 1.0-(Pi+Pd);
    else if(d1==d0+1) return Pi;
    else              return Pd;
  } // if d0
}

//================================================================================
double InnerCodec::FuncPpqx(int p0, int q0, int x, int p1, int q1){
  assert(p0>=0 && p0<NumSTR);
  assert(p1>=0 && p1<NumSTR);
  assert(q0>=0 && q0<NumSTB);
  assert(q1>=0 && q1<NumSTB);
  assert(x>=0  && x<4);
  if(STR[p0][x]==p1 && STB[q0][x]==q1) return 1.0;
  return 0.0;
}

//================================================================================
void InnerCodec::SetSTR(){
  int st;
  STR = new int * [NumSTR];
  for(int i=0;i<NumSTR;i++) STR[i] = new int [4];
  // initial state
  for(int a=0;a<4;a++) STR[4*RunK][a] = a*RunK;
  // other states
  for(int a=0;a<4;a++){
    for(int i=0;i<RunK;i++){
      st = a*RunK+i;
      for(int a1=0;a1<4;a1++){
	if(a1==a){
	  STR[st][a1] = (i==RunK-1)? -1 : st+1;
	} else {
	  STR[st][a1] = RunK*a1;
	} // if a1==a
      } // for a1
    } // for i
  } // for a
}

//================================================================================
void InnerCodec::SetSTB(){
  STB = new int * [NumSTB];
  for(int i=0;i<NumSTB;i++) STB[i] = new int [4];
  // inital state
  for(int a=0;a<4;a++) STB[2*Delta+1][a] = (a<=1)? -1+Delta : +1+Delta;
  // boundary states
  for(int a=0;a<4;a++) STB[-Delta+Delta][a] = (a<=1)? -1 : -Delta+1+Delta;  // left-end
  for(int a=0;a<4;a++) STB[ Delta+Delta][a] = (a<=1)? Delta-1+Delta : -1;   // right-end
  // other states
  for(int b=-Delta+1;b<Delta;b++){
    for(int a=0;a<4;a++) STB[b+Delta][a] = (a<=1)? b-1+Delta : b+1+Delta;
  } // for b
}

//================================================================================
void InnerCodec::SetGT(){
  // int L; // dest set
  // int z; // Fenc(x)
  double p0,p1;
  //----- GT1
  GT1 = new double *** [4];
  for(int y=0;y<4;y++){
    GT1[y] = new double ** [4];
    for(int x=0;x<4;x++){
      GT1[y][x] = new double * [NumSTR];
      for(int p=0; p<NumSTR;p++){
	GT1[y][x][p] = new double [NumSTB];
	for(int q=0;q<NumSTB;q++){
	  // L = CalcL(p,q);
	  // z = Fenc[L][x];
	  // assert(z>=-1 && z<4);
	  // if(z==-1) GT1[y][x][p][q] = 0.0;  // (invalid)
	  // else      GT1[y][x][p][q] = (y==z)? 1.0-Ps : Ps/3.0;
	  GT1[y][x][p][q] = (y==x)? 1.0-Ps : Ps/3.0;
	} // for q
      } // for p
    } // for x
  } // for y
  //----- GT2
  GT2 = new double **** [4];
  for(int y0=0;y0<4;y0++){
    GT2[y0] = new double *** [4];
    for(int y1=0;y1<4;y1++){
      GT2[y0][y1] = new double ** [4];
      for(int x=0;x<4;x++){
	GT2[y0][y1][x] = new double * [NumSTR];
	for(int p=0; p<NumSTR;p++){
	  GT2[y0][y1][x][p] = new double [NumSTB];
	  for(int q=0;q<NumSTB;q++){
	    // L = CalcL(p,q);
	    // z = Fenc[L][x];
	    // assert(z>=-1 && z<4);
	    // if(z==-1){
	    //   GT2[y0][y1][x][p][q] = 0.0;  // (invalid)
	    // } else {
	    //   p0 = (y0==z)? 1.0-Ps : Ps/3.0;
	    //   p1 = (y1==z)? 1.0-Ps : Ps/3.0;
	    //   GT2[y0][y1][x][p][q] = p0*p1;
	    // } // if z
	    p0 = (y0==x)? 1.0-Ps : Ps/3.0;
	    p1 = (y1==x)? 1.0-Ps : Ps/3.0;
	    GT2[y0][y1][x][p][q] = p0*p1;
	  } // for q
	} // for p
      } // for x
    } // for y1
  } // for y0
}

//================================================================================
void InnerCodec::AllocFG(){
  // PU, PD
  PU = new double * [N];
  PD = new double * [N];
  for(int i=0;i<N;i++){
    PU[i] = new double [4];
    PD[i] = new double [4];
  } // for i
  // PF, PB
  PF = new double *** [N+1];
  PB = new double *** [N+1];
  for(int i=0;i<N+1;i++){
    PF[i] = new double ** [NumSTR];
    PB[i] = new double ** [NumSTR];
    for(int p=0;p<NumSTR;p++){
      PF[i][p] = new double * [NumSTB];
      PB[i][p] = new double * [NumSTB];
      for(int q=0;q<NumSTB;q++){
	PF[i][p][q] = new double [NumD];
	PB[i][p][q] = new double [NumD];
      } // for q
    } // for p
  } // for i
}

//================================================================================
void InnerCodec::DelFG(){
  // PU, PD
  for(int i=0;i<N;i++){
    delete [] PU[i];
    delete [] PD[i];
  } // for i
  delete [] PU;
  delete [] PD;
  // PF, PB
  for(int i=0;i<N+1;i++){
    for(int p=0;p<NumSTR;p++){
      for(int q=0;q<NumSTB;q++){
	delete [] PF[i][p][q];
	delete [] PB[i][p][q];
      } // for q
      delete [] PF[i][p];
      delete [] PB[i][p];
    } // for p
    delete [] PF[i];
    delete [] PB[i];
  } // for i
  delete [] PF;
  delete [] PB;
} 

//================================================================================
void InnerCodec::InitFG(const int *RW, const double **Px, int N2){
  assert(N2>=N+Dmin && N2<=N+Dmax);
  ClearNode();
  // -----Y
  memcpy(Y,RW,sizeof(int)*N2);
  for(int i=0;i<N2;i++) assert(Y[i]>=0 && Y[i]<4);
  // -----PU
  for(int i=0;i<N;i++){
    memcpy(PU[i],Px[i],sizeof(double)*4);
    Normalize(PU[i],4);
  } // for i
  // -----PF[0]
  PF[0][NumSTR-1][NumSTB-1][0-Dmin] = 1.0; // init state
  // ----- PB[N]
  for(int p=0;p<NumSTR;p++){
    for(int q=0;q<NumSTB;q++){
      for(int d=Dmin;d<=Dmax;d++){
	if(d==N2-N && p!=NumSTR-1 && q!=NumSTB-1){
	  PB[N][p][q][d-Dmin] = 1.0/((NumSTR-1)*(NumSTB-1));
	} else {
	  PB[N][p][q][d-Dmin] = 0.0;
	}
      } // for d
    } // for q
  } // for p
  Normalize(PB[N],NumSTR,NumSTB,NumD);
}

//================================================================================
void InnerCodec::ClearNode(int idx){
  assert(idx>=0 && idx<=N);
  if(idx<N){
    ClearVect(PU[idx],4);
    ClearVect(PD[idx],4);
  } // if idx
  for(int p=0;p<NumSTR;p++){
    for(int q=0;q<NumSTB;q++){
      ClearVect(PF[idx][p][q],NumD);
      ClearVect(PB[idx][p][q],NumD);
    } // for q
  } // for p
}

//================================================================================
void InnerCodec::ClearNode(){
  for(int i=0;i<=N;i++) ClearNode(i);
}

//================================================================================
void InnerCodec::Forward(int idx, int N2){
  assert(idx>=0 && idx<N);
  double s,ss,sd,sy;
  for(int p1=0;p1<NumSTR;p1++){
    for(int q1=0;q1<NumSTB;q1++){
      for(int d1=Dmin;d1<=Dmax;d1++){
	if(idx+d1<0 || idx+d1>=N2) continue;
	s = 0;
	for(int d0=d1-1;d0<=d1+1;d0++){
	  if(idx+d0<0 || idx+d0>=N2) continue;
	  sd = FuncPd(d0,d1);
	  ss = 0.0;
	  if(sd==0) continue;
	  for(int x=0;x<4;x++){
	    for(int p0=0;p0<NumSTR;p0++){
	      for(int q0=0;q0<NumSTB;q0++){
		if( FuncPpqx(p0,q0,x,p1,q1)==0.0 ) continue;
		if(d1==d0){
		  // --- trans
		  sy = GT1[Y[idx+d0]][x][p0][q0];
		} else if(d1==d0+1) {
		  // --- ins
		  sy = GT2[Y[idx+d0]][Y[idx+d0+1]][x][p0][q0];
		} else if(d1==d0-1) {
		  // --- del
		  sy = 1.0;
		} else {
		  assert(false);
		} // id d1
		ss += (sy*PF[idx][p0][q0][d0-Dmin]*PU[idx][x]);
		//(dbg)
		/*
		printf("idx=%d (p1,q1)=(%d,%d) (d0,d1)=(%d,%d) (p0,q0)=(%d,%d) x=%d sd=%e sy=%e ss=%e s=%e\n",
		       idx,p1,q1,d0,d1,p0,q0,x,sd,sy,ss,s);
		*/
	      } // for q0
	    } // for p0
	  } // for x
	  s += (ss*sd);
	} // for d0
	PF[idx+1][p1][q1][d1-Dmin] = s;
	// ---
      } // for d1
    } // for q1
  } // for p1
  Normalize(PF[idx+1],NumSTR,NumSTB,NumD);
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
InnerCodec::InnerCodec(int _RunK, int _Delta, double _Pi, double _Pd, double _Ps, int _N, int _Dmin, int _Dmax){
  RunK  = _RunK;
  Delta = _Delta;
  Pi    = _Pi;
  Pd    = _Pd;
  Ps    = _Ps;
  N     = _N;
  Dmin  = _Dmin;
  Dmax  = _Dmax;
  NumSTR = 4*RunK+1;
  NumSTB = 2*Delta+2;
  NumD   = Dmax-Dmin+1;
  Nmax   = N+Dmax;
  printf("# InnerCodec: RunK=%d Delta=%d Pids=(%e,%e,%e) N=%d D=(%d,%d)\n",RunK,Delta,Pi,Pd,Ps,N,Dmin,Dmax);
  printf("# InnerCodec: NumSTR=%d NumSTB=%d NumD=%d Nmax=%d\n",NumSTR,NumSTB,NumD,Nmax);
  assert(RunK>0 && Delta>0);
  assert(Pi>=0.0 && Pi<1.0);
  assert(Pd>=0.0 && Pd<1.0);
  assert(Ps>=0.0 && Ps<1.0);
  assert(N>0);
  assert(Dmin<=0 && Dmax>=0);
  
  SetSTR();
  SetSTB();
  SetGT();
  AllocFG();
  Y = new int [Nmax];
  
  // (dbg)
  // PrintSTR();
  // PrintSTB();
  // PrintFenc();
  PrintGT1();
  // PrintGT2();
}

//================================================================================
InnerCodec::~InnerCodec(){
  DelFG();
  delete [] Y;
  for(int i=0;i<NumSTR;i++) delete [] STR[i];
  for(int i=0;i<NumSTB;i++) delete [] STB[i];
  delete [] STR;
  delete [] STB;
  for(int y=0;y<4;y++){
    for(int x=0;x<4;x++){
      for(int p=0; p<NumSTR;p++){
	delete [] GT1[y][x][p];
      } // for p
      delete [] GT1[y][x];
    } // for x
    delete [] GT1[y];
  } // for y
  delete [] GT1;
  
  printf("# InneCodec: deleted\n");
}

//================================================================================
void InnerCodec::Encode(int *CW, const int *IW){
  int p = NumSTR-1; // STR init
  int q = NumSTB-1; // STB init
  int Z,L;
  for(int i=0;i<N;i++){
    assert(IW[i]>=0 && IW[i]<4);
    L = CalcL(p,q);
    Z = Fenc[L][IW[i]];
    CW[i] = Z;
    //printf("%03d: (%02d,%02d;%02d) %02d->%02d\n",i,p,q,L,IW[i],CW[i]);
    p = STR[p][Z];
    q = STB[q][Z];
    assert(p>=0 && p<NumSTR);
    assert(q>=0 && q<NumSTB);
  } // for i
}

//================================================================================
void InnerCodec::Decode(double **Pout, const int *RW, const double **Px, int N2, const int *DbgCW){
  InitFG(RW,Px,N2);
  for(int idx=0;idx<N;idx++) Forward(idx,N2);

  // Forward(0,N2);
  // Forward(1,N2);
  // Forward(2,N2);
  
  //(dbg)
  PrintNode(DbgCW);
  // PrintNode(0, DbgCW[0]);
  // PrintNode(1, DbgCW[1]);
  // PrintNode(2, DbgCW[2]);
  // PrintNode(3, DbgCW[3]);
}

//================================================================================
void InnerCodec::PrintSTR(){
  printf("# STR\n");
  for(int i=0;i<NumSTR;i++){
    printf("#  %02d: ",i);
    for(int a=0;a<4;a++) printf("%02d ",STR[i][a]);
    printf("\n");
  } // for i
}

//================================================================================
void InnerCodec::PrintSTB(){
  printf("# STB\n");
  for(int i=0;i<NumSTB;i++){
    printf("#  %02d: ",i);
    for(int a=0;a<4;a++) printf("%02d ",STB[i][a]);
    printf("\n");
  } // for i
}

//================================================================================
void InnerCodec::PrintFenc(){
  printf("# Fenc\n");
  for(int i=0;i<16;i++){
    printf("#  %02d: ", i);
    for(int j=0;j<4;j++) printf("%d",BitGet(i,j));
    printf(" ");
    for(int j=0;j<4;j++) printf("%d ",Fenc[i][j]);
    printf("\n");
  } // for i
}

//================================================================================
void InnerCodec::PrintGT1(){
  printf("# GT1:\n");
  for(int p=0;p<NumSTR;p++){
    for(int q=0;q<NumSTB;q++){
      printf("#  p=%02d q=%02d: ",p,q);
      for(int x=0;x<4;x++){
	printf("[%d] ",x);
	for(int y=0;y<4;y++) printf("%.2e ",GT1[y][x][p][q]);
      } // for x
      printf("\n");
    } // for q
  } // for p
}

//================================================================================
void InnerCodec::PrintGT2(){
  printf("# GT2:\n");
  for(int p=0;p<NumSTR;p++){
    for(int q=0;q<NumSTB;q++){
      printf("#  p=%02d q=%02d:\n",p,q);
      for(int x=0;x<4;x++){
	printf("#   [%d] ",x);
	for(int y0=0;y0<4;y0++){
	  for(int y1=0;y1<4;y1++) printf("%.2e ",GT2[y0][y1][x][p][q]);
	  printf("| ");
	} // for y0
	printf("\n");	
      } // for x
    } // for q
  } // for p
}

//================================================================================
void InnerCodec::PrintNode(int idx, int DbgCW){
  assert(idx>=0 && idx<=N);
  printf("[%04d] %d > %d\n",idx,DbgCW,Y[idx]);
  if(idx<N){
    PrintVect(PU[idx],4,"PU: ","\n");
    PrintVect(PD[idx],4,"PD: ","\n");
  } else {
    printf("PU:\nPD:\n");
  } // if idx
  printf("PF:\n");
  for(int p=0;p<NumSTR;p++){
    for(int q=0;q<NumSTB;q++){
      if(SumVect(PF[idx][p][q],NumD)>0){
	printf("%02d:%02d ",p,q);
	PrintVect(PF[idx][p][q],NumD,"","\n");
      } //if
    } // for q
  } // for p
  printf("PB:\n");
  for(int p=0;p<NumSTR;p++){
    for(int q=0;q<NumSTB;q++){
      if(SumVect(PB[idx][p][q],NumD)>0){
	printf("%02d:%02d ",p,q);
	PrintVect(PB[idx][p][q],NumD,"","\n");
      } //if
    } // for q
  } // for p
}

//================================================================================
void InnerCodec::PrintNode(const int *DbgCW){
  int c;
  for(int i=0;i<=N;i++){
    c = (DbgCW==NULL)? -1 : DbgCW[i];
    PrintNode(i,c);
  } // for i
}
