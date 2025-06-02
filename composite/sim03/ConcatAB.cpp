#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "InnerCodebook.hpp"
#include "CSFBAdec.hpp"
#include "CPScodec.hpp"
#include "ConcatAB.hpp"


//================================================================================
void ConcatAB::ClearVect(double *V, int len){
  for(int i=0; i<len; i++) V[i]=0.0;
}

//================================================================================
void ConcatAB::Normalize(double *V, int len){
  assert( len>0 );
  double s=0;
  for(int i=0; i<len; i++){
    assert( V[i] >= 0.0 );
    s += V[i];
  } // for i
  assert( s>0.0 );
  for(int i=0; i<len; i++) V[i] /= s;
}

//================================================================================
void ConcatAB::ConvPrbQB(){
  unsigned char *CW = new unsigned char [Nu];
  double p;
  // init
  for(int i=0; i<Ns; i++) Normalize( PrBQ[i], Qa );
  for(int i=0; i<Nb; i++) ClearVect( PrB2[i], 2  );
  // convert
  for(int idx=0; idx<Ns; idx++){
    for(int v=0; v<Qa; v++){
      ICB->Get_CW(CW,v);
      p = PrBQ[idx][v];
      for(int j=0; j<Nu; j++){
	assert( CW[j]==0 || CW[j]==1 );
	PrB2[idx*Nu+j][ CW[j] ] += p;
      } // for j
    } // for v
  } // for idx
  delete [] CW;
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
ConcatAB::ConcatAB(class InnerCodebook *_ICB, class CSFBAdec *_DECa, class CPScodec *_DECb){
  ICB  = _ICB;
  DECa = _DECa;
  DECb = _DECb;
  Nb   = DECa->Get_Nb();
  Nu   = DECa->Get_Nu();
  Ns   = DECa->Get_Ns();
  Nseq = DECa->Get_Nseq();
  Qa   = DECa->Get_Q();
  Qb   = DECb->Get_QB();
  Dmin = DECa->Get_Dmin();
  Dmax = DECa->Get_Dmax();
  Drng = DECa->Get_Drng();
  printf("# ConcatAB: (Nb,Nu,Ns)=(%d,%d,%d) Nseq=%d (Qa,Qb)=(%d,%d) (Dmin,Dmax,Drng)=(%d,%d,%d)\n",
	 Nb,Nu,Ns,Nseq,Qa,Qb,Dmin,Dmax,Drng);
  assert( DECb->Get_Nb()   == Nb );
  assert( DECb->Get_Nu()   == Nu );
  assert( DECb->Get_Ns()   == Ns );
  assert( DECb->Get_Nseq() == Nseq );
  assert( DECb->Get_Dmin() == Dmin );
  assert( DECb->Get_Dmax() == Dmax );
  assert( DECb->Get_Drng() == Drng );
  assert( ICB->Get_Nu()    == Nu );
  assert( ICB->Get_numCW() == Qa );
  //-----
  PrD = new double ** [Nseq];
  for(int i=0; i<Nseq; i++){
    PrD[i] = new double * [Ns+1];
    for(int j=0; j<Ns+1; j++) PrD[i][j] = new double [Drng];
  } // for i  
  PrB2 = new double * [Nb];
  PrBQ = new double * [Ns];
  for(int i=0; i<Nb; i++) PrB2[i] = new double [2];
  for(int i=0; i<Ns; i++) PrBQ[i] = new double [Qa];
}

//================================================================================
ConcatAB::~ConcatAB(){
  for(int i=0; i<Nseq; i++){
    for(int j=0; j<Ns+1; j++) delete [] PrD[i][j];
    delete [] PrD[i];
  } // for i
  delete [] PrD;
  for(int i=0; i<Nb; i++) delete [] PrB2[i];
  for(int i=0; i<Ns; i++) delete [] PrBQ[i];
  delete [] PrB2;
  delete [] PrBQ;
  printf("# ConcatAB: deleted\n");
}

//================================================================================
void ConcatAB::AtoB(){
  DECa->Get_PFB( PrD );
  DECa->Get_PO( PrBQ );
  ConvPrbQB();    // PrB2 <- PrBQ
  DECb->InitFG();
  DECb->SetFGD( (const double ***)PrD  );
  DECb->SetFGB( (const double  **)PrB2 );
  //(dbg)
  //DECb->PrintFG();
}

//================================================================================
void ConcatAB::Get_PrB2(double **PrB2out){
  for(int i=0; i<Nb; i++) memcpy( PrB2out[i], PrB2[i], sizeof(double)*2 );
}
