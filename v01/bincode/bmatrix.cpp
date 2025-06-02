#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <string>
#include <cstring>
#include <fstream>
#include <iostream>

#include "bmatrix.hpp"


//============================================================================
void bmatrix::GetSize(int *Mx, int *Nx, const char *fn, bool bin){
  std::ifstream fin;
  std::string   buf;
  std::string   dlm = " \t,";
  char *cbuf;

  if(bin){  //----- binary
    fin.open(fn, std::ios::in|std::ios::binary);
    if(!fin){
      std::cout << "Cannot open " << fn << std::endl;
      exit(1);
    }
    fin.read((char *)Mx,sizeof(int));
    fin.read((char *)Nx,sizeof(int));
    fin.close();
  } else {  //----- alist  
    fin.open(fn, std::ios::in);
    if(!fin){
      std::cout << "Cannot open " << fn << std::endl;
      exit(1);
    }
    std::getline(fin,buf);
    cbuf = new char [buf.length()+1];
    strcpy(cbuf,buf.c_str());
    *Nx = atoi(strtok(cbuf," \t,"));
    *Mx = atoi(strtok(NULL," \t,\n"));
    delete [] cbuf;
  }
}

//============================================================================
void bmatrix::bitset(unsigned char *u, int pos, int s){
  // MSB=7,6,5,4,3,2,1,0=LSB
  assert(pos>=0 && pos<8);
  assert(s==0 || s==1);
  unsigned char ptn = 0x01<<pos;
  if(s==1) {
    (*u) = (*u) | ptn; 
  } else {
    ptn = ptn^0xFF;
    (*u) = (*u) & ptn;
  }
}

//============================================================================
int bmatrix::bitget(unsigned char u, int pos){
  // MSB=7,6,5,4,3,2,1,0=LSB
  assert(pos>=0 && pos<8);
  unsigned char ptn = 0x01<<pos;
  return ((u&ptn)==0)? 0 : 1;
}

//============================================================================
bool bmatrix::bittest(){
  int i,x;
  unsigned char a,b;
  for(x=0;x<256;x++){
    a = (unsigned char)x;
    b = 0x00;
    for(i=0;i<8;i++) bitset(&b,i,bitget(a,i));
    //printf("%02X %02X\n",a,b);
    if(a!=b){
      printf("bittest: a=%02X b=%02X\n",a,b);
      return false;
    }
  } // for x
  return true;
}

//============================================================================
int bmatrix::byteweight(unsigned char u, int len){
  assert(len>=0 && len<=8);
  int i,w=0;
  for(i=0;i<len;i++){
    if((u & 0x01)!=0) w++;
    u>>=1;
  }
  return w;
}

//============================================================================
//============================================================================
//============================================================================

//============================================================================
bmatrix::bmatrix(int _M, int _N){
  M = _M;
  N = _N;
  assert(M>0 && N>=0);
  Nb   = (int)ceil((double)N/8.0);
  Nbs  = sizeof(unsigned char)*Nb;
  Nfrac= N-(Nb-1)*8;
  //printf("# bmatrix: %dx%d(Nb=%d Nbs=%d)\n",M,N,Nb,Nbs);
  V = new unsigned char * [M];
  for(int i=0;i<M;i++) V[i] = new unsigned char [Nb];
  clear();
  assert(bittest());
  //printf("# bmatrix: generated\n");
}

//============================================================================
bmatrix::bmatrix(const char *fn, bool bin){
  GetSize(&M,&N,fn,bin);
  assert(M>0 && N>=0);
  Nb = (int)ceil((double)N/8.0);
  Nbs= sizeof(unsigned char)*Nb;  
  //printf("# bmatrix: %dx%d(Nb=%d Nbs=%d)\n",M,N,Nb,Nbs);
  V = new unsigned char * [M];
  for(int i=0;i<M;i++) V[i] = new unsigned char [Nb];
  clear();
  assert(bittest());
  //printf("# bmatrix: generated\n");
  if(bin) read(fn);
  else    read_alist(fn);
}

//============================================================================
bmatrix::~bmatrix(){
  for(int i=0;i<M;i++) delete [] V[i];
  delete [] V;
  //printf("# bmatrix: deleted\n");
}

//============================================================================
int bmatrix::getM(){ return M;}
int bmatrix::getN(){ return N;}
int bmatrix::getNb(){return Nb;}

//============================================================================
void bmatrix::clear(){
  for(int i=0;i<M;i++)
    for(int j=0;j<Nb;j++)
      V[i][j]=0x00;
}

//============================================================================
void bmatrix::randmat(){
  clear();
  for(int i=0;i<M;i++)
    for(int j=0;j<Nb;j++)
      //V[i][j]=(unsigned char)(random() % 256);
      V[i][j]=(unsigned char)floor( ((double)random()/RAND_MAX)*256.0 );
}

//============================================================================
void bmatrix::randvec(int w){
  assert(w>0 && w<=N);
  assert(M==1);
  int pos,wt=0;
  clear();
  while(wt<w){
    pos = (int)floor( ((double)random()/((long)RAND_MAX+1))*N );
    if(getV(0,pos)==0){
      setV(0,pos,1);
      wt++;
    }
  } // while
}

//============================================================================
void bmatrix::print(){
  for(int i=0;i<M;i++){
    printf("%04d: ",i);
    for(int j=0;j<N;j++){
      printf("%d",getV(i,j));
      if(j%8==7) printf(" ");
    } // for j
    printf("\n");
  } // for i
}

//============================================================================
int bmatrix::getV(int r, int c){
  assert(r>=0 && r<M);
  assert(c>=0 && c<N);
  int ptn = 0x00000007;
  int ci  = c>>3;
  int co  = c & ptn;
  return bitget(V[r][ci],co);
}

//============================================================================
void bmatrix::setV(int r, int c, int s){
  assert(r>=0 && r<M);
  assert(c>=0 && c<N);
  assert(s==0 || s==1);
  int ptn = 0x00000007;
  int ci  = c>>3;
  int co  = c & ptn;
  bitset(&V[r][ci],co,s);
}

//============================================================================
void bmatrix::setV(class bmatrix *H){
  assert(H->M==M);
  assert(H->N==N);
  for(int i=0;i<M;i++)
    for(int j=0;j<Nb;j++)
      V[i][j]=H->V[i][j];
}

//============================================================================
void bmatrix::setV(class bmatrix *H, int rt, int ct, int rs, int cs, int lr, int lc){
  int Ms = H->getM();
  int Ns = H->getN();
  assert(lr>=0 && lc>=0);
  assert(rt>=0 && rt+lr<=M);
  assert(ct>=0 && ct+lc<=N);
  assert(rs>=0 && rs+lr<=Ms);
  assert(cs>=0 && cs+lc<=Ns);
  for(int i=0;i<lr;i++)
    for(int j=0;j<lc;j++) setV(rt+i,ct+j,H->getV(rs+i,cs+j));
}

//============================================================================
void bmatrix::invV(int r, int c){
  assert(r>=0 && r<M);
  assert(c>=0 && c<N);
  setV(r,c,(getV(r,c)==0)?1:0);
}

//============================================================================
void bmatrix::getRV(unsigned char *RV, int r){
  assert(r>=0 && r<M);
  memcpy(RV,V[r],Nbs);
}

//============================================================================
void bmatrix::getRV(class bmatrix *RV, int r){
  assert(r>=0 && r<M);
  assert(RV->getM()==1);
  assert(RV->getN()==N);
  memcpy(RV->V[0],V[r],Nbs);
}

//============================================================================
void bmatrix::setRV(int r, const unsigned char *RV){
  assert(r>=0 && r<M);
  memcpy(V[r],RV,Nbs);
}

//============================================================================
void bmatrix::setRV(int r, class bmatrix *RV){
  assert(r>=0 && r<M);
  assert(RV->getM()==1);
  assert(RV->getN()==N);
  memcpy(V[r],RV->V[0],Nbs);
}

//============================================================================
void bmatrix::getRVunp(unsigned char *RV, int r){
  assert(r>=0 && r<M);
  for(int i=0;i<N;i++) RV[i] = getV(r,i);
}

//============================================================================
void bmatrix::setRVunp(int r, const unsigned char *RV){
  assert(r>=0 && r<M);
  for(int i=0;i<N;i++) setV(r,i,RV[i]);
}

//============================================================================
void bmatrix::setIdentity(int r, int c, int l){
  assert(l>0);
  assert(r>=0 && r+l<=M);
  assert(c>=0 && c+l<=N);
  for(int i=0;i<l;i++)
    for(int j=0;j<l;j++)
      setV(r+i,c+j,(i==j)?1:0); 
}

//============================================================================
void bmatrix::RowExchange(int r0, int r1){
  assert(r0>=0 && r0<M);
  assert(r1>=0 && r1<M);
  unsigned char *Vtmp = new unsigned char [Nb];
  memcpy(Vtmp, V[r0],Nbs);
  memcpy(V[r0],V[r1],Nbs);
  memcpy(V[r1],Vtmp, Nbs);
  delete [] Vtmp;
}

//============================================================================
void bmatrix::RowAdd(int r0, int r1){
  assert(r0>=0 && r0<M);
  assert(r1>=0 && r1<M);
  assert(r0!=r1);
  for(int j=0;j<Nb;j++) V[r0][j] = V[r0][j] ^ V[r1][j];
}

//============================================================================
void bmatrix::RowAdd(int r, const unsigned char *RV){
  assert(r>=0 && r<M);
  for(int i=0;i<Nb;i++)
    V[r][i] = V[r][i]^RV[i];
}

//============================================================================
void bmatrix::ColExchange(int c0, int c1){
  int s;
  assert(c0>=0 && c0<N);
  assert(c1>=0 && c1<N);
  for(int i=0;i<M;i++){
    s = getV(i,c0);
    setV(i,c0,getV(i,c1));
    setV(i,c1,s);
  }
}

//============================================================================
void bmatrix::GenCyclicMat(class bmatrix *RV){
  assert(RV->getN()==N);
  assert(RV->getM()==1);
  memcpy(V[0],RV->V[0],Nbs);
  for(int i=1;i<M;i++){
    setV(i,0,getV(i-1,N-1));
    for(int j=1;j<N;j++)
      setV(i,j,getV(i-1,j-1));
  } // for i
}

//============================================================================
void bmatrix::ShiftRight(int r){
  int i,val;
  assert(r>=0 && r<M);
  val = bitget(V[r][Nb-1],7-(Nb*8-N));
  for(i=Nb-1;i>0;i--){
    V[r][i]<<=1;
    bitset(&V[r][i],0,bitget(V[r][i-1],7));
  } // for i
  V[r][0]<<=1;
  bitset(&V[r][0],0,val);
}

//============================================================================
void bmatrix::ShiftLeft(int r){
  int i,val;
  assert(r>=0 && r<M);
  val = bitget(V[r][0],0);
  for(i=0;i<Nb-1;i++){
    V[r][i]>>=1;
    bitset(&V[r][i],7,bitget(V[r][i+1],0));
  } // for i
  V[r][Nb-1]>>=1;
  bitset(&V[r][Nb-1],7-(Nb*8-N),val);
}

//============================================================================
bool bmatrix::isIdentity(){
  if(N!=M) return false;
  for(int i=0;i<M;i++)
    for(int j=0;j<N;j++)
      if((i==j && getV(i,j)==0) || (i!=j && getV(i,j)==1)) return false;
  return true;
}

//============================================================================
bool bmatrix::isZero(){
  for(int i=0;i<M;i++)
    for(int j=0;j<N;j++)
      if(getV(i,j)!=0) return false;
  return true;
}

//============================================================================
bool bmatrix::isRcyclic(){
  for(int i=0;i<M-1;i++){
    for(int j=0;j<N-1;j++)
      if(getV(i,j)!=getV(i+1,j+1)) return false;
    if(getV(i,N-1)!=getV(i+1,0)) return false;
  } // for i
  return true;
}

//============================================================================
bool bmatrix::isLcyclic(){
  for(int i=0;i<M-1;i++){
    for(int j=1;j<N;j++)
      if(getV(i,j)!=getV(i+1,j-1)) return false;
    if(getV(i,0)!=getV(i+1,N-1)) return false;
  } // for i
  return true;
}

//============================================================================
int bmatrix::HammingWeight(){
  int w=0;
  for(int i=0;i<M;i++){
    for(int j=0;j<Nb-1;j++)
      w+=byteweight(V[i][j],8);
    w+=byteweight(V[i][Nb-1],Nfrac);
  }
  return w;
}

//============================================================================
int bmatrix::HammingWeightOld(){
  int w=0;
  for(int i=0;i<M;i++)
    for(int j=0;j<N;j++)
      if(getV(i,j)==1) w++;
  return w;
}

//============================================================================
void bmatrix::setAnd(class bmatrix *H0, class bmatrix *H1){
  assert(H0->getM()==M && H0->getN()==N);
  assert(H1->getM()==M && H1->getN()==N);
  for(int i=0;i<M;i++){
    for(int j=0;j<Nb;j++) V[i][j] = H0->V[i][j] & H1->V[i][j];
  } // for i
}

//============================================================================
void bmatrix::write(const char *fn){
  std::ofstream fout;
  fout.open(fn, std::ios::out|std::ios::binary|std::ios::trunc);
  if(!fout){
    std::cout << "Cannot open " << fn << std::endl;
    exit(1);
  }
  fout.write((char *)&M,sizeof(int));
  fout.write((char *)&N,sizeof(int));
  for(int i=0;i<M;i++)
    for(int j=0;j<Nb;j++)
      fout.write((char *)&V[i][j],sizeof(unsigned char));
  fout.close();
}

//============================================================================
void bmatrix::read(const char *fn){
  int Mx,Nx;
  std::ifstream fin;
  fin.open(fn, std::ios::in|std::ios::binary);
  if(!fin){
    std::cout << "Cannot open " << fn << std::endl;
    exit(1);
  }
  fin.read((char *)&Mx,sizeof(int));
  fin.read((char *)&Nx,sizeof(int));
  assert(Mx==M);
  assert(Nx==N);
  for(int i=0;i<M;i++)
    for(int j=0;j<Nb;j++)
      fin.read((char *)&V[i][j],sizeof(unsigned char));
  fin.close();
}

//============================================================================
void bmatrix::write_alist(const char *fn){
  int i,j;
  int CWmax=0, RWmax=0;
  int *CWlist = new int [N];
  int *RWlist = new int [M];

  // CWlist, CWmax
  for(j=0;j<N;j++){
    CWlist[j]=0;
    for(i=0;i<M;i++)
      if(getV(i,j)==1) CWlist[j]++;
    if(CWlist[j]>CWmax) CWmax=CWlist[j];
  } // for j

  // RWlist, RWmax
  for(i=0;i<M;i++){
    RWlist[i]=0;
    for(j=0;j<N;j++)
      if(getV(i,j)==1) RWlist[i]++;
    if(RWlist[i]>RWmax) RWmax=RWlist[i];
  } // for i
  
  std::ofstream fout;
  fout.open(fn, std::ios::out|std::ios::trunc);
  if(!fout){
    std::cout << "Cannot open " << fn << std::endl;
    exit(1);
  }
  // L1: N M
  fout << N << " " << M << std::endl;
  // L2: CWmax RWmax
  fout << CWmax << " " << RWmax << std::endl;
  // L3: CWlist
  for(j=0;j<N;j++) fout << CWlist[j] << " ";
  fout << std::endl;
  // L4: RWlist
  for(i=0;i<M;i++) fout << RWlist[i] << " ";
  fout << std::endl;
  // columns
  for(j=0;j<N;j++){
    for(i=0;i<M;i++)
      if(getV(i,j)==1) fout << i+1 << " ";
    fout << std::endl;
  }
  // rows
  for(i=0;i<M;i++){
    for(j=0;j<N;j++)
      if(getV(i,j)==1) fout << j+1 << " ";
    fout << std::endl;
  }
  
  fout.close();
  delete [] CWlist;
  delete [] RWlist;
}

//============================================================================
void bmatrix::read_alist( const char *fn){
  int pos,Mx,Nx;
  std::ifstream fin;
  std::string   buf;
  std::string   dlm = " \t,";
  char *cbuf, *ctok;
  
  clear();
  fin.open(fn, std::ios::in);
  if(!fin){
    std::cout << "Cannot open " << fn << std::endl;
    exit(1);
  }
  // L1: N M
  std::getline(fin,buf);
  cbuf = new char [buf.length()+1];
  strcpy(cbuf,buf.c_str());
  Nx = atoi(strtok(cbuf," \t,"));
  Mx = atoi(strtok(NULL," \t,\n"));
  assert(Nx==N && Mx==M);
  delete [] cbuf;
  // L2: CWmax RWmax (unused)
  std::getline(fin,buf);
  // L3: CWlist (unused)
  std::getline(fin,buf);
  // L4: RWlist (unused)
  std::getline(fin,buf);
  // Cols
  for(int j=0;j<N;j++){
    std::getline(fin,buf);
    cbuf = new char [buf.length()+1];
    strcpy(cbuf,buf.c_str());
    ctok = strtok(cbuf," \t,\n");
    pos = atoi(ctok)-1;
    assert(pos>=0 && pos<M);
    setV(pos,j,1);
    while((ctok=strtok(NULL," \t,\n"))!=NULL){
      pos = atoi(ctok)-1;
      assert(pos>=0 && pos<M);
      setV(pos,j,1);
    }
    delete [] cbuf;
  } // for j
  // Rows (consistency check)
  class bmatrix *Hx = new class bmatrix(M,N);
  Hx->clear();
  for(int i=0;i<M;i++){
    std::getline(fin,buf);
    cbuf = new char [buf.length()+1];
    strcpy(cbuf,buf.c_str());
    ctok = strtok(cbuf," \t,\n");
    pos = atoi(ctok)-1;
    assert(pos>=0 && pos<N);
    Hx->setV(i,pos,1);
    while((ctok=strtok(NULL," \t,\n"))!=NULL){
      pos = atoi(ctok)-1;
      assert(pos>=0 && pos<N);
      Hx->setV(i,pos,1);
    }
    delete [] cbuf;
  } // for j
  assert(bmatrix_compare(this,Hx)==0);
  delete Hx;
  
  fin.close();
}

//============================================================================
//============================================================================
//============================================================================


//============================================================================
int bmatrix_compare(class bmatrix *H0, class bmatrix *H1){
  int M = H0->getM();
  int N = H0->getN();
  assert(H1->getM()==M);
  assert(H1->getN()==N);
  
  for(int i=0;i<M;i++){
    for(int j=0;j<N;j++)
      if(H0->getV(i,j) != H1->getV(i,j)) return 1;
  } // for i
  return 0;
}

//============================================================================
bool bmatrix_inverse(class bmatrix *HI, class bmatrix *H){
  int i,j;
  int N = H->getN();
  bool fullrank=true;
  assert(H->getM() ==N);
  assert(HI->getM()==N);
  assert(HI->getN()==N);
  HI->clear();
  
  // set Hx = [ H | I ]
  class bmatrix *Hx = new class bmatrix(N,2*N);
  Hx->clear();
  Hx->setV(H,0,0,0,0,N,N);
  for(i=0;i<N;i++) Hx->setV(i,N+i,1);

  //H->print();
  //Hx->print();
  
  // elementary row operations
  for(i=0;i<N;i++){
    // 1: find nonzero element in the i-th column
    j=i;
    while(Hx->getV(j,i)==0){
      j++;
      if(j==N){
	fullrank=false;
	break;
      }
    } // while
    if(!fullrank) break;

    // 2: row exhange
    Hx->RowExchange(i,j);

    // 3: make the other elements 0
    for(j=0;j<N;j++){
      if(j!=i && Hx->getV(j,i)==1) Hx->RowAdd(j,i);
    } // for j

    //printf("[%d/%d]\n",i,N-1);
    //Hx->print();
    
  } // for i
  
  if(fullrank){
    HI->setV(Hx,0,0,0,N,N,N);
  } else {
    printf("bmatrix_inverse: matrix is not fullrank\n");
  }

  //Hx->print();
  //HI->print();
  
  delete Hx;
  return fullrank;
}

//============================================================================
void bmatrix_transpose(class bmatrix *HT, class bmatrix *H){
  int M=H->getM();
  int N=H->getN();
  assert(HT->getM()==N);
  assert(HT->getN()==M);
  for(int i=0;i<M;i++)
    for(int j=0;j<N;j++)
      HT->setV(j,i,H->getV(i,j));
}

//============================================================================
void bmatrix_mul(class bmatrix *H, class bmatrix *H0, class bmatrix *H1){
  int M = H->getM();
  int N = H->getN();
  int K = H0->getN();
  int Nb= H->getNb();
  int s;
  assert(H0->getM()==M);
  assert(H1->getN()==N);
  assert(H1->getM()==K);
  unsigned char *RV = new unsigned char [Nb];
  
  H->clear();
  for(int i=0;i<M;i++){
    for(int j=0;j<K;j++){
      s = H0->getV(i,j);
      if(s==1){
	H1->getRV(RV,j);
	H->RowAdd(i,RV);
      } // if s==1
    } // for j
  } // for i

  delete [] RV;
}

//============================================================================
void bmatrix_add(class bmatrix *H, class bmatrix *H0, class bmatrix *H1){
  int M = H->getM();
  int N = H->getN();
  int Nb= H->getNb();
  unsigned char *RV  = new unsigned char [Nb];
  unsigned char *RV0 = new unsigned char [Nb];
  unsigned char *RV1 = new unsigned char [Nb];
  assert(H0->getM()==M && H0->getN()==N);
  assert(H1->getM()==M && H1->getN()==N);
  
  for(int i=0;i<M;i++){
    H0->getRV(RV0,i);
    H1->getRV(RV1,i);
    for(int j=0;j<Nb;j++) RV[j] = RV0[j] ^ RV1[j];
    H->setRV(i,RV);
  } // for i

  delete [] RV;
  delete [] RV0;
  delete [] RV1;
}

//============================================================================
void bmatrix_and(class bmatrix *H, class bmatrix *H0, class bmatrix *H1){
  int M = H->getM();
  int N = H->getN();
  int Nb= H->getNb();
  unsigned char *RV  = new unsigned char [Nb];
  unsigned char *RV0 = new unsigned char [Nb];
  unsigned char *RV1 = new unsigned char [Nb];
  assert(H0->getM()==M && H0->getN()==N);
  assert(H1->getM()==M && H1->getN()==N);
  
  for(int i=0;i<M;i++){
    H0->getRV(RV0,i);
    H1->getRV(RV1,i);
    for(int j=0;j<Nb;j++) RV[j] = RV0[j] & RV1[j];
    H->setRV(i,RV);
  } // for i

  delete [] RV;
  delete [] RV0;
  delete [] RV1;
}

//============================================================================
int bmatrix_HammingDist(class bmatrix *H0, class bmatrix *H1){
  int d=0;
  int M=H0->getM();
  int N=H0->getN();
  assert(H1->getM()==M && H1->getN()==N);

  for(int i=0;i<M;i++)
    for(int j=0;j<N;j++)
      if(H0->getV(i,j)!=H1->getV(i,j)) d++;
  
  return d;
}
