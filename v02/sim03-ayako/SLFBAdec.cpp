#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

#include "InnerCodebook.hpp"
#include "ChannelMatrix.hpp"
#include "IDSchannel.hpp"
#include "SLFBAdec.hpp"

#define NuMax 10

//================================================================================
int SLFBAdec::max(int a, int b){return (a>b)? a : b;}
int SLFBAdec::min(int a, int b){return (a<b)? a : b;}

//================================================================================
double SLFBAdec::Psub(unsigned char a, unsigned char b){
  assert( a==0 || a==1 );
  assert( b==0 || b==1 );
  return (a==b)? 1.0-Ps : Ps;
}

//================================================================================
void SLFBAdec::ClearMat(double **X, int M0, int N0){
  for(int i=0;i<M0;i++){
    for(int j=0;j<N0;j++) X[i][j]=0.0;
  } // for i
}

//================================================================================
void SLFBAdec::CopyMat(double **Y, const double **X, int M0, int N0){
  assert(M0>0 && N0>0);
  for(int i=0;i<M0;i++) memcpy(Y[i],X[i],sizeof(double)*N0);
}

//================================================================================
void SLFBAdec::PrintMat(const double **X, int M0, int N0, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<M0;i++){
    printf("%03d: ",i);
    for(int j=0;j<N0;j++) printf("%.2e ",X[i][j]);
    printf("\n");
  } // for i
  printf("%s",post);
}

//================================================================================
void SLFBAdec::MultMat(double **Y, const double **X0, const double **X1, int M0, int N0, int M1, int N1){
  assert(M0>0 && N0>0 && M1>0 && N1>0);
  assert(N0==M1);
  double **Y0 = new double * [M0];
  for(int i=0;i<M0;i++){
    Y0[i] = new double [N1];
    for(int j=0;j<N1;j++) Y0[i][j]=0.0;
  } // for i
  //--- calc
  for(int i=0;i<M0;i++){
    for(int j=0;j<N1;j++){
      for(int k=0;k<N0;k++){
	Y0[i][j] += X0[i][k]*X1[k][j];
      } // for k
    } // for j
  } // for i
  //--- Y <- Y0
  for(int i=0;i<M0;i++) memcpy(Y[i],Y0[i],sizeof(double)*N1);
  //--- del
  for(int i=0;i<M0;i++) delete [] Y0[i];
  delete [] Y0;
}

//================================================================================
void SLFBAdec::PrintVect(const unsigned char *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%u",V[i]);
  printf("%s",post);
}

//================================================================================
void SLFBAdec::PrintVect(const int *V, int len, const char *pre, const char *post){
  printf("%s",pre);
  for(int i=0;i<len;i++) printf("%d ",V[i]);
  printf("%s",post);
}

//================================================================================
void SLFBAdec::PrintVect(const double *V, int len, const char *pre, const char *post){
  int w=20;
  printf("%s",pre);
  for(int i=0;i<len;i++){
    printf("%.2e ",V[i]);
    if(i%w==w-1) printf("\n");
  } // for i
  printf("%s",post);
}

//================================================================================
long SLFBAdec::VectToLong(const unsigned char *V, int len){
  assert(len>=0);
  long val = 0;
  for(int i=0;i<len;i++){
    assert(V[i]==0 || V[i]==1);
    val <<= 1;
    if(V[i]==1) val |= 0x1;
  } // for i
  return val;
}

//================================================================================
void SLFBAdec::LongToVect(unsigned char *V, long val, int len){
  assert(len>0 && val>=0);
  long mask = 0x1 << (len-1);
  for(int i=0;i<len;i++){
    V[i] = ( (val & mask)==0 )? 0 : 1;
    mask >>= 1;
  } // for i
}

//--------------------------------------------------------------------------------
void SLFBAdec::normalize(double *V, int len){
  assert(len>0);
  double s=0;
  for(int i=0;i<len;i++){
    assert( V[i]>=0.0 );
    s += V[i];
  } // for i
  assert( s>0 );
  for(int i=0;i<len;i++) V[i] /= s;
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

//================================================================================
void SLFBAdec::SetGD(){
  // p(d)　ドリフト遷移確率の設定
  double **GDX = new double * [Drng];  // for calc
  GD = new double * [Drng];
  for(int i=0;i<Drng;i++){
    GD[i]  = new double [Drng];
    GDX[i] = new double [Drng];
  } // for i
  // --- GDX (1-bit)
  ClearMat(GDX,Drng,Drng);
  GDX[0][0] = 1.0-Pi;
  GDX[Drng-1][Drng-1] = 1.0-Pd;
  for(int i=0;i<Drng-1;i++) GDX[i][i+1] = Pi;
  for(int i=1;i<Drng;  i++) GDX[i][i-1] = Pd;
  for(int i=1;i<Drng-1;i++) GDX[i][i]   = 1.0-(Pi+Pd);
  // --- GD <- identity mat
  ClearMat(GD,Drng,Drng);
  for(int i=0;i<Drng;i++) GD[i][i] = 1.0;
  // --- calc
  for(int i=0;i<Nu;i++) MultMat(GD,(const double**)GD,(const double **)GDX,Drng,Drng,Drng,Drng);
  // (dbg)
  // PrintMat((const double **)GDX,Drng,Drng,"GDX\n","");
  // PrintMat((const double **)GD, Drng,Drng,"GD\n", "");
  // --- del
  for(int i=0;i<Drng;i++) delete [] GDX[i];
  delete [] GDX;
}

//================================================================================
void SLFBAdec::DelGD(){
  for(int i=0;i<Drng;i++) delete [] GD[i];
  delete [] GD;
}

//================================================================================
void SLFBAdec::SetGX(){
  // 符号語とチャネル出力の関係を事前計算
  // Nu = 送信シンボル長
  // Nu2 = 受診シンボル長
  long ly2p,x;
  unsigned char *X = new unsigned char [Nu];
  GX = new double ** [Nu*2+1];
  for(int ly=0;ly<=Nu*2;ly++){
    if(ly<Nu2min || ly>Nu2max){
      // approximate
      GX[ly]    = new double * [1];
      GX[ly][0] = new double [1];
      GX[ly][0][0] = (ly<Nu)? pow(Pd,Nu-ly) : pow(Pi,ly-Nu); //?
      //printf("GX[%d]=%e\n",ly,GX[ly][0][0]);
    } else {
      // exact
      ly2p = (long)pow(2,ly);
      GX[ly] = new double * [ly2p];
      for(long y=0;y<ly2p;y++){
	GX[ly][y] = new double [Q];
	for(int xi=0;xi<Q;xi++){
	  ICB->Get_CW(X,xi);
	  x = VectToLong(X,Nu);
	  GX[ly][y][xi] = CalcPyx(y,x,ly,Nu);
	  //(dbg)
	  //printf("%e ly=%d y=%ld xi=%d x=%ld ",GX[ly][y][xi],ly,y,xi,x);
	  //PrintVect(X,Nu,"X=","\n");
	} // for x
      } // for y
    } // if ly
  } // for ly
  delete [] X;
}

//================================================================================
void SLFBAdec::DelGX(){
  int i2p;
  for(int i=0;i<=Nu*2;i++){
    if(i<Nu2min || i>Nu2max){
      delete [] GX[i][0];
      delete [] GX[i];
    } else {
      i2p = (long)pow(2,i);
      for(long y=0;y<i2p;y++) delete [] GX[i][y];
      delete [] GX[i];
    } // if i
  } // for i
  delete [] GX;
}

//================================================================================
double SLFBAdec::CalcPyx(long y, long x, int ly, int lx){
  // チャネル確率p(y|x,d)の計算
  assert(lx>0 && ly>=0);
  //printf("y=%ld x=%ld ly=%d lx=%d\n",y,x,ly,lx);
  if(ly==0) return pow(Pd,lx);
  long   x1,y1;
  double ret,qt,qi,qd;
  unsigned char *X = new unsigned char [lx];
  unsigned char *Y = new unsigned char [ly];
  LongToVect(X,x,lx);
  LongToVect(Y,y,ly);
  if(lx==1){
    if(ly==1){
      ret = Psub(X[0],Y[0]) * Pt;
    } else if(ly==2) {
      ret = Psub(X[0],Y[0]) * Psub(X[0],Y[1]) * Pi;
    } else {
      assert(ly>=3); 
      ret = 0.0;
    } // if ly
  } else {
    x1  = VectToLong(&X[1],lx-1);
    // -----trans
    qt = Psub(X[0],Y[0]) * Pt;
    y1 = (ly==1)? 0 : VectToLong(&Y[1],ly-1);
    qt *= CalcPyx(y1,x1,ly-1,lx-1);
    // -----ins
    if(ly<2){
      qi = 0.0;
    } else {
      qi = Psub(X[0],Y[0]) * Psub(X[0],Y[1]) * Pi;
      y1 = (ly==2)? 0 : VectToLong(&Y[2],ly-2);
      qi *= CalcPyx(y1,x1,ly-2,lx-1);	
    } // if ly<2
    // -----del
    qd = Pd;
    y1 = y;
    qd *= CalcPyx(y1,x1,ly,lx-1);
    // ---
    ret = qt + qi + qd;
  } // if lx
  delete [] X;
  delete [] Y;
  //printf("ret=%e\n",ret);
  return ret;
}

//================================================================================
double SLFBAdec::GetGX(int Nu2, long y, long xi){
  assert(Nu2>=0 && Nu2<=2*Nu);
  assert(y >=0);
  assert(xi>=0 && xi<Q);
  if(Nu2<Nu2min || Nu2>Nu2max) return GX[Nu2][0][0 ]; // approx
  else                         return GX[Nu2][y][xi]; 
}

//================================================================================
void SLFBAdec::SetGXNew(){
  // 新しい確率分布の事前計算: P(y|φ_0(u'_i), e_iv, d_iv, d_(i+1)v)
  long ly2p,x;
  unsigned char *X = new unsigned char [Nu];
  GXNew = new double *** [Nu*2+1];
  for(int ly=0;ly<=Nu*2;ly++){
    if(ly<Nu2min || ly>Nu2max){
      // approximate
      GXNew[ly]    = new double ** [1];
      GXNew[ly][0] = new double * [1];
      GXNew[ly][0][0] = new double [4];  // 4つのエラー状態
      // 各エラー状態への適切な確率設定
      if(ly == 0) {
        // 削除のみが可能（受信長0）
        GXNew[ly][0][0][0] = 0.0;              // Match (不可能)
        GXNew[ly][0][0][1] = 0.0;              // Insertion (不可能)
        GXNew[ly][0][0][2] = pow(Pd,Nu);       // Deletion (全削除)
        GXNew[ly][0][0][3] = 0.0;              // Substitution (不可能)
      } else {
        // その他の近似計算
        double base_prob = (ly<Nu)? pow(Pd,Nu-ly) : pow(Pi,ly-Nu);
        GXNew[ly][0][0][0] = base_prob * 0.7;  // Match
        GXNew[ly][0][0][1] = base_prob * 0.1;  // Insertion
        GXNew[ly][0][0][2] = base_prob * 0.1;  // Deletion
        GXNew[ly][0][0][3] = base_prob * 0.1;  // Substitution
      }
    } else {
      // exact
      ly2p = (long)pow(2,ly);
      GXNew[ly] = new double ** [ly2p];
      for(long y=0;y<ly2p;y++){
        GXNew[ly][y] = new double * [Q];
        for(int xi=0;xi<Q;xi++){
          ICB->Get_CW(X,xi);
          x = VectToLong(X,Nu);
          GXNew[ly][y][xi] = new double [4];  // 4つのエラー状態
          // 修正されたCalcPyxNew関数を使用
          // 各エラー状態の確率を個別に計算せず、総和を各状態に分配
          double total_prob = CalcPyxNew(y,x,ly,Nu,0); // 全状態の総確率
          
          // lattice.cppに基づく状態分布（近似）
          double state_weights[4] = {0.7, 0.1, 0.1, 0.1}; // Match, Ins, Del, Sub
          for(int e=0; e<4; e++){
            GXNew[ly][y][xi][e] = total_prob * state_weights[e];
          }
        } // for xi
      } // for y
    } // if ly
  } // for ly
  delete [] X;
}

//================================================================================
void SLFBAdec::DelGXNew(){
  int i2p;
  for(int i=0;i<=Nu*2;i++){
    if(i<Nu2min || i>Nu2max){
      delete [] GXNew[i][0][0];
      delete [] GXNew[i][0];
      delete [] GXNew[i];
    } else {
      i2p = (long)pow(2,i);
      for(long y=0;y<i2p;y++){
        for(int xi=0;xi<Q;xi++){
          delete [] GXNew[i][y][xi];
        }
        delete [] GXNew[i][y];
      }
      delete [] GXNew[i];
    } // if i
  } // for i
  delete [] GXNew;
}

//================================================================================
double SLFBAdec::CalcPyxNew(long y, long x, int ly, int lx, int error_state){
  // lattice.cppに基づく格子計算による正確なチャネル確率p(y|x,e)の計算
  assert(lx>0 && ly>=0);
  assert(error_state>=0 && error_state<4);
  
  unsigned char *X = new unsigned char [lx];
  unsigned char *Y = new unsigned char [ly];
  LongToVect(X,x,lx);
  LongToVect(Y,y,ly);
  
  // 格子計算用の3次元配列 F[n][m][e]
  std::vector<std::vector<std::vector<double>>> F(lx+1, 
    std::vector<std::vector<double>>(ly+1, std::vector<double>(4, 0.0)));
  
  // 格子の初期化
  initializeLattice(F, lx, ly);
  
  // 動的プログラミングによる格子計算
  calculateLatticeDP(F, X, Y, lx, ly);
  
  delete [] X;
  delete [] Y;
  
  // 指定されたエラー状態での終端確率を返す
  return F[lx][ly][error_state];
}

//================================================================================
void SLFBAdec::initializeLattice(std::vector<std::vector<std::vector<double>>>& F, int lx, int ly){
  // 初期状態の設定（lattice.cppと同じ）
  // F[0][0][0] = 1.0でMatch状態から開始
  F[0][0][0] = 1.0; // Match状態から開始
  F[0][0][1] = 0.0; // Substitution
  F[0][0][2] = 0.0; // Deletion  
  F[0][0][3] = 0.0; // Insertion
}

//================================================================================
double SLFBAdec::getTransitionProbability(int prev_state, int curr_state){
  // lattice.cppの遷移確率行列と同じ値
  double P_transition[4][4] = {
    {0.7, 0.1, 0.1, 0.1},  // Match からの遷移
    {0.6, 0.2, 0.1, 0.1},  // Substitution からの遷移  
    {0.5, 0.1, 0.3, 0.1},  // Deletion からの遷移
    {0.5, 0.1, 0.1, 0.3}   // Insertion からの遷移
  };
  
  return P_transition[prev_state][curr_state];
}

//================================================================================
void SLFBAdec::calculateLatticeDP(std::vector<std::vector<std::vector<double>>>& F, 
                                  const unsigned char* X, const unsigned char* Y, int lx, int ly){
  // lattice.cppと同じ格子計算ロジック
  
  for(int n = 0; n <= lx; n++){
    for(int m = 0; m <= ly; m++){
      
      // Match計算（対角線移動: n-1,m-1 -> n,m）
      if(n > 0 && m > 0 && X[n-1] == Y[m-1]){
        for(int prev_e = 0; prev_e < 4; prev_e++){
          double prev_prob = F[n-1][m-1][prev_e];
          double trans_prob = getTransitionProbability(prev_e, 0); // Match=0
          double match_prob = Psub(X[n-1], Y[m-1]) * Pt;
          F[n][m][0] += prev_prob * trans_prob * match_prob;
        }
      }
      
      // Substitution計算（対角線移動: n-1,m-1 -> n,m）
      if(n > 0 && m > 0 && X[n-1] != Y[m-1]){
        for(int prev_e = 0; prev_e < 4; prev_e++){
          double prev_prob = F[n-1][m-1][prev_e];
          double trans_prob = getTransitionProbability(prev_e, 1); // Substitution=1
          double sub_prob = Psub(X[n-1], Y[m-1]) * Pt * (1.0/3.0);
          F[n][m][1] += prev_prob * trans_prob * sub_prob;
        }
      }
      
      // Deletion計算（垂直移動: n-1,m -> n,m）
      if(n > 0){
        for(int prev_e = 0; prev_e < 4; prev_e++){
          double prev_prob = F[n-1][m][prev_e];
          double trans_prob = getTransitionProbability(prev_e, 2); // Deletion=2
          double del_prob = Pd;
          F[n][m][2] += prev_prob * trans_prob * del_prob;
        }
      }
      
      // Insertion計算（水平移動: n,m-1 -> n,m）
      if(m > 0){
        for(int prev_e = 0; prev_e < 4; prev_e++){
          double prev_prob = F[n][m-1][prev_e];
          double trans_prob = getTransitionProbability(prev_e, 3); // Insertion=3
          double ins_prob = Pi * (1.0/4.0); // 4つの塩基への等確率挿入
          F[n][m][3] += prev_prob * trans_prob * ins_prob;
        }
      }
      
    } // for m
  } // for n
}

//================================================================================
void SLFBAdec::debugLatticeCalculation(long y, long x, int ly, int lx, const char* filename){
  // 特定のy,xペアに対する格子計算の詳細なデバッグ出力
  
  unsigned char *X = new unsigned char [lx];
  unsigned char *Y = new unsigned char [ly];
  LongToVect(X,x,lx);
  LongToVect(Y,y,ly);
  
  std::ofstream file(filename);
  if(!file.is_open()){
    printf("Error: Cannot open debug file %s\n", filename);
    delete [] X;
    delete [] Y;
    return;
  }
  
  file << "=== Lattice Calculation Debug ===" << std::endl;
  file << "Input (x): ";
  for(int i=0; i<lx; i++) file << (int)X[i];
  file << " (binary), " << x << " (decimal)" << std::endl;
  file << "Output (y): ";
  for(int i=0; i<ly; i++) file << (int)Y[i];
  file << " (binary), " << y << " (decimal)" << std::endl;
  file << "Length: lx=" << lx << ", ly=" << ly << std::endl;
  file << std::endl;
  
  // 格子計算用の3次元配列 F[n][m][e]
  std::vector<std::vector<std::vector<double>>> F(lx+1, 
    std::vector<std::vector<double>>(ly+1, std::vector<double>(4, 0.0)));
  
  // 格子の初期化
  F[0][0][0] = 1.0; // Match状態から開始
  file << "=== Initialization ===" << std::endl;
  file << "F[0][0][0] = 1.0 (Match state)" << std::endl;
  file << std::endl;
  
  // 遷移確率行列の出力
  file << "=== Transition Probability Matrix ===" << std::endl;
  file << "P_transition[prev_state][curr_state]:" << std::endl;
  file << "        M      S      D      I" << std::endl;
  file << "M:   0.70   0.10   0.10   0.10" << std::endl;
  file << "S:   0.60   0.20   0.10   0.10" << std::endl;
  file << "D:   0.50   0.10   0.30   0.10" << std::endl;
  file << "I:   0.50   0.10   0.10   0.30" << std::endl;
  file << std::endl;
  
  // 動的プログラミングによる格子計算（詳細ログ付き）
  file << "=== Dynamic Programming Calculation ===" << std::endl;
  
  for(int n = 0; n <= lx; n++){
    for(int m = 0; m <= ly; m++){
      
      if(n == 0 && m == 0) continue; // 初期状態はスキップ
      
      file << "--- Processing F[" << n << "][" << m << "] ---" << std::endl;
      if(n > 0 && m > 0){
        file << "Comparing: X[" << (n-1) << "]=" << (int)X[n-1] 
             << " vs Y[" << (m-1) << "]=" << (int)Y[m-1];
        file << " -> " << (X[n-1] == Y[m-1] ? "MATCH" : "MISMATCH") << std::endl;
      }
      
      double old_values[4] = {F[n][m][0], F[n][m][1], F[n][m][2], F[n][m][3]};
      
      // Match計算（対角線移動: n-1,m-1 -> n,m）
      if(n > 0 && m > 0 && X[n-1] == Y[m-1]){
        file << "  Match calculation:" << std::endl;
        for(int prev_e = 0; prev_e < 4; prev_e++){
          double prev_prob = F[n-1][m-1][prev_e];
          double trans_prob = getTransitionProbability(prev_e, 0); // Match=0
          double match_prob = Psub(X[n-1], Y[m-1]) * Pt;
          double contribution = prev_prob * trans_prob * match_prob;
          F[n][m][0] += contribution;
          
          if(prev_prob > 0){
            file << "    From state " << prev_e << ": " << prev_prob 
                 << " * " << trans_prob << " * " << match_prob 
                 << " = " << contribution << std::endl;
          }
        }
      }
      
      // Substitution計算（対角線移動: n-1,m-1 -> n,m）
      if(n > 0 && m > 0 && X[n-1] != Y[m-1]){
        file << "  Substitution calculation:" << std::endl;
        for(int prev_e = 0; prev_e < 4; prev_e++){
          double prev_prob = F[n-1][m-1][prev_e];
          double trans_prob = getTransitionProbability(prev_e, 1); // Substitution=1
          double sub_prob = Psub(X[n-1], Y[m-1]) * Pt * (1.0/3.0);
          double contribution = prev_prob * trans_prob * sub_prob;
          F[n][m][1] += contribution;
          
          if(prev_prob > 0){
            file << "    From state " << prev_e << ": " << prev_prob 
                 << " * " << trans_prob << " * " << sub_prob 
                 << " = " << contribution << std::endl;
          }
        }
      }
      
      // Deletion計算（垂直移動: n-1,m -> n,m）
      if(n > 0){
        file << "  Deletion calculation:" << std::endl;
        for(int prev_e = 0; prev_e < 4; prev_e++){
          double prev_prob = F[n-1][m][prev_e];
          double trans_prob = getTransitionProbability(prev_e, 2); // Deletion=2
          double del_prob = Pd;
          double contribution = prev_prob * trans_prob * del_prob;
          F[n][m][2] += contribution;
          
          if(prev_prob > 0){
            file << "    From state " << prev_e << ": " << prev_prob 
                 << " * " << trans_prob << " * " << del_prob 
                 << " = " << contribution << std::endl;
          }
        }
      }
      
      // Insertion計算（水平移動: n,m-1 -> n,m）
      if(m > 0){
        file << "  Insertion calculation:" << std::endl;
        for(int prev_e = 0; prev_e < 4; prev_e++){
          double prev_prob = F[n][m-1][prev_e];
          double trans_prob = getTransitionProbability(prev_e, 3); // Insertion=3
          double ins_prob = Pi * (1.0/4.0); // 4つの塩基への等確率挿入
          double contribution = prev_prob * trans_prob * ins_prob;
          F[n][m][3] += contribution;
          
          if(prev_prob > 0){
            file << "    From state " << prev_e << ": " << prev_prob 
                 << " * " << trans_prob << " * " << ins_prob 
                 << " = " << contribution << std::endl;
          }
        }
      }
      
      // 結果の出力
      file << "  Results: F[" << n << "][" << m << "] = [";
      for(int e = 0; e < 4; e++){
        if(e > 0) file << ", ";
        file << std::scientific << std::setprecision(6) << F[n][m][e];
      }
      file << "]" << std::endl;
      
      // 変化量の出力
      bool changed = false;
      for(int e = 0; e < 4; e++){
        if(F[n][m][e] != old_values[e]) changed = true;
      }
      if(changed){
        file << "  Changes: [";
        for(int e = 0; e < 4; e++){
          if(e > 0) file << ", ";
          double change = F[n][m][e] - old_values[e];
          file << std::scientific << std::setprecision(6) << change;
        }
        file << "]" << std::endl;
      }
      file << std::endl;
      
    } // for m
  } // for n
  
  // 最終結果
  file << "=== Final Results ===" << std::endl;
  file << "F[" << lx << "][" << ly << "] = [";
  for(int e = 0; e < 4; e++){
    if(e > 0) file << ", ";
    file << std::scientific << std::setprecision(6) << F[lx][ly][e];
  }
  file << "]" << std::endl;
  file << std::endl;
  
  file << "Error state probabilities:" << std::endl;
  file << "Match (0):       " << std::scientific << std::setprecision(6) << F[lx][ly][0] << std::endl;
  file << "Substitution (1): " << std::scientific << std::setprecision(6) << F[lx][ly][1] << std::endl;
  file << "Deletion (2):    " << std::scientific << std::setprecision(6) << F[lx][ly][2] << std::endl;
  file << "Insertion (3):   " << std::scientific << std::setprecision(6) << F[lx][ly][3] << std::endl;
  
  file.close();
  delete [] X;
  delete [] Y;
  
  printf("Lattice debug output saved to: %s\n", filename);
}

//================================================================================
double SLFBAdec::GetGXNew(int Nu2, long y, long xi, int error_state){
  assert(Nu2>=0 && Nu2<=2*Nu);
  assert(y >=0);
  assert(xi>=0 && xi<Q);
  assert(error_state>=0 && error_state<4);
  
  if(Nu2<Nu2min || Nu2>Nu2max) return GXNew[Nu2][0][0][error_state]; // approx
  else                         return GXNew[Nu2][y][xi][error_state]; 
}

//================================================================================
void SLFBAdec::SetGENew(){
  // 新しい確率分布の事前計算: P(e_(i+1)v|e_iv, φ0(u'_(i-⌈k/v⌉+1)), φ0(u'_(i-⌈k/v⌉)), φ0(u'_i), φ0(u'_(i+1)), d_iv, d_(i+1)v)
  long ly2p,x;
  unsigned char *X = new unsigned char [Nu];
  GENew = new double *** [Nu*2+1];
  for(int ly=0;ly<=Nu*2;ly++){
    if(ly<Nu2min || ly>Nu2max){
      // approximate
      GENew[ly]    = new double ** [1];
      GENew[ly][0] = new double * [1];
      GENew[ly][0][0] = new double [4];  // 4つのエラー状態
      // k-merコンテキスト依存の基本確率（簡略版）
      double base_prob = (ly<Nu)? pow(Pd,Nu-ly) : pow(Pi,ly-Nu);
      // 遷移確率行列（eLattice_v7.cppと同じ）
      GENew[ly][0][0][0] = base_prob * 0.7;  // Match状態継続確率
      GENew[ly][0][0][1] = base_prob * 0.1;  // Insertion状態遷移確率  
      GENew[ly][0][0][2] = base_prob * 0.1;  // Deletion状態遷移確率
      GENew[ly][0][0][3] = base_prob * 0.1;  // Substitution状態遷移確率
    } else {
      // exact
      ly2p = (long)pow(2,ly);
      GENew[ly] = new double ** [ly2p];
      for(long y=0;y<ly2p;y++){
        GENew[ly][y] = new double * [Q];
        for(int xi=0;xi<Q;xi++){
          ICB->Get_CW(X,xi);
          x = VectToLong(X,Nu);
          GENew[ly][y][xi] = new double [4];  // 4つのエラー状態
          // k-merコンテキスト依存の確率計算
          double base_prob = CalcPexNew(y,x,ly,Nu);
          // 各エラー状態への遷移確率を設定
          GENew[ly][y][xi][0] = base_prob * 0.7;  // Match状態継続
          GENew[ly][y][xi][1] = base_prob * 0.1;  // Insertion状態遷移
          GENew[ly][y][xi][2] = base_prob * 0.1;  // Deletion状態遷移
          GENew[ly][y][xi][3] = base_prob * 0.1;  // Substitution状態遷移
        } // for xi
      } // for y
    } // if ly
  } // for ly
  delete [] X;
}

//================================================================================
void SLFBAdec::DelGENew(){
  int i2p;
  for(int i=0;i<=Nu*2;i++){
    if(i<Nu2min || i>Nu2max){
      delete [] GENew[i][0][0];
      delete [] GENew[i][0];
      delete [] GENew[i];
    } else {
      i2p = (long)pow(2,i);
      for(long y=0;y<i2p;y++){
        for(int xi=0;xi<Q;xi++){
          delete [] GENew[i][y][xi];
        }
        delete [] GENew[i][y];
      }
      delete [] GENew[i];
    } // if i
  } // for i
  delete [] GENew;
}

//================================================================================
double SLFBAdec::CalcPexNew(long y, long x, int ly, int lx){
  // 実験データベースに基づくk-merコンテキスト依存のエラー状態遷移確率計算
  // P(e_(i+1)v|e_iv, φ0(u'_(i-⌈k/v⌉+1)), φ0(u'_(i-⌈k/v⌉)), φ0(u'_i), φ0(u'_(i+1)), d_iv, d_(i+1)v)
  assert(lx>0 && ly>=0);
  
  if(ly==0) return pow(Pd,lx);
  
  long x1,y1;
  double ret,qt,qi,qd;
  unsigned char *X = new unsigned char [lx];
  unsigned char *Y = new unsigned char [ly];
  LongToVect(X,x,lx);
  LongToVect(Y,y,ly);
  
  if(lx==1){
    if(ly==1){
      // k-merを抽出してテーブルから確率を取得
      std::string current_kmer = extractKmerFromSequence(X, lx, 0, k_mer_length);
      
      // デフォルト確率（実験データがない場合）
      double match_prob = 0.7, sub_prob = 0.1;
      
      // 実験データから確率を取得
      if (transition_probs.find(current_kmer) != transition_probs.end()) {
        // 前のイベントを0(Match)と仮定、より正確には前の状態を追跡する必要がある
        auto& probs = transition_probs[current_kmer][0]; // Match から遷移と仮定
        match_prob = probs[0]; // Match確率
        sub_prob = probs[3];   // Substitution確率
      }
      
      if(X[0]==Y[0]){
        ret = Psub(X[0],Y[0]) * Pt * match_prob;  // Match状態継続
      } else {
        ret = Psub(X[0],Y[0]) * Pt * sub_prob * (1.0/3.0);  // Substitution状態遷移
      }
    } else if(ly==2) {
      // k-merを抽出してテーブルから挿入確率を取得
      std::string current_kmer = extractKmerFromSequence(X, lx, 0, k_mer_length);
      double ins_prob = 0.1; // デフォルト
      
      if (transition_probs.find(current_kmer) != transition_probs.end()) {
        auto& probs = transition_probs[current_kmer][0]; // Match から遷移と仮定
        ins_prob = probs[1]; // Insertion確率
      }
      
      ret = Psub(X[0],Y[0]) * Psub(X[0],Y[1]) * Pi * ins_prob * (1.0/4.0);  // Insertion状態遷移
    } else {
      assert(ly>=3); 
      ret = 0.0;
    }
  } else {
    x1 = VectToLong(&X[1],lx-1);
    
    // k-merを抽出してテーブルから確率を取得
    std::string current_kmer = extractKmerFromSequence(X, lx, 0, k_mer_length);
    
    // デフォルト確率
    double match_prob = 0.7, ins_prob = 0.1, del_prob = 0.1, sub_prob = 0.1;
    
    // 実験データから確率を取得
    if (transition_probs.find(current_kmer) != transition_probs.end()) {
      auto& probs = transition_probs[current_kmer][0]; // Match から遷移と仮定
      match_prob = probs[0]; // Match確率
      ins_prob = probs[1];   // Insertion確率
      del_prob = probs[2];   // Deletion確率
      sub_prob = probs[3];   // Substitution確率
    }
    
    // k-merコンテキストを考慮したMatch/Substitution
    if(X[0]==Y[0]){
      qt = Psub(X[0],Y[0]) * Pt * match_prob;  // Match継続
    } else {
      qt = Psub(X[0],Y[0]) * Pt * sub_prob * (1.0/3.0);  // Substitution遷移
    }
    y1 = (ly==1)? 0 : VectToLong(&Y[1],ly-1);
    qt *= CalcPexNew(y1,x1,ly-1,lx-1);
    
    // k-merコンテキストを考慮したInsertion
    if(ly<2){
      qi = 0.0;
    } else {
      qi = Psub(X[0],Y[0]) * Psub(X[0],Y[1]) * Pi * ins_prob * (1.0/4.0);  // Insertion遷移
      y1 = (ly==2)? 0 : VectToLong(&Y[2],ly-2);
      qi *= CalcPexNew(y1,x1,ly-2,lx-1);
    }
    
    // k-merコンテキストを考慮したDeletion  
    qd = Pd * del_prob;  // Deletion遷移
    y1 = y;
    qd *= CalcPexNew(y1,x1,ly,lx-1);
    
    ret = qt + qi + qd;
  }
  
  delete [] X;
  delete [] Y;
  return ret;
}

//================================================================================
double SLFBAdec::GetGENew(int Nu2, long y, long xi){
  assert(Nu2>=0 && Nu2<=2*Nu);
  assert(y >=0);
  assert(xi>=0 && xi<Q);
  
  // この関数はエラー状態の事後確率の合計を返す（正規化前）
  double total_prob = 0.0;
  
  if(Nu2<Nu2min || Nu2>Nu2max) {
    // approximation: 全エラー状態の確率を合計
    for(int e=0; e<4; e++) {
      total_prob += GENew[Nu2][0][0][e];
    }
  } else {
    // exact: 全エラー状態の確率を合計
    for(int e=0; e<4; e++) {
      total_prob += GENew[Nu2][y][xi][e];
    }
  }
  
  return total_prob;
}

//================================================================================
void SLFBAdec::loadTransitionProbabilities(int k) {
  // 確率テーブルを読み込み (eLattice_v7.cppと同じ方法)
  std::string base_path = "DNArSim-main/simulator/probEdit/k" + std::to_string(k) + "/";
  
  std::vector<std::string> event_files = {
    "KmerYi_prevYiM_RatesAvg.txt", // Match
    "KmerYi_prevYiI_RatesAvg.txt", // Insertion  
    "KmerYi_prevYiD_RatesAvg.txt", // Deletion
    "KmerYi_prevYiS_RatesAvg.txt"  // Substitution
  };
  
  std::vector<std::string> event_names = {"Match", "Insertion", "Deletion", "Substitution"};
  
  for (int prev_event = 0; prev_event < 4; prev_event++) {
    std::string full_path = base_path + event_files[prev_event];
    printf("# Loading: %s (prev event: %s)\n", full_path.c_str(), event_names[prev_event].c_str());
    
    std::ifstream file(full_path);
    if (!file.is_open()) {
      printf("# Warning: Cannot open file %s\n", full_path.c_str());
      continue;
    }
    
    std::string line;
    int line_count = 0;
    
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      std::string kmer;
      double ins, del, sub, err, match;
      
      if (iss >> kmer >> ins >> del >> sub >> err >> match) {
        transition_probs[kmer][prev_event] = {match, ins, del, sub}; // Match, Ins, Del, Sub の順
        line_count++;
      }
    }
    
    printf("# Loaded %d k-mers for %s\n", line_count, event_names[prev_event].c_str());
    file.close();
  }
  
  printf("# Total k-mers in transition table: %zu\n", transition_probs.size());
}

//================================================================================
std::string SLFBAdec::binaryToDNA(const unsigned char* binary_seq, int length) {
  // バイナリペアをDNA塩基に変換するテーブル
  std::map<std::pair<int,int>, char> binary_to_dna = {
    {{0,0}, 'A'}, {{0,1}, 'T'}, {{1,0}, 'G'}, {{1,1}, 'C'}
  };
  
  std::string dna_seq;
  for (int i = 0; i < length; i += 2) {
    if (i + 1 < length) {
      std::pair<int,int> pair = {binary_seq[i], binary_seq[i+1]};
      if (binary_to_dna.find(pair) != binary_to_dna.end()) {
        dna_seq += binary_to_dna[pair];
      } else {
        dna_seq += 'A'; // デフォルト
      }
    }
  }
  return dna_seq;
}

//================================================================================
std::string SLFBAdec::extractKmerFromSequence(const unsigned char* seq, int seq_len, int pos, int k) {
  // バイナリ配列からk-merを抽出してDNA配列に変換
  if (pos < 0 || pos >= seq_len || seq_len < k*2) {
    return std::string(k, 'A'); // デフォルト
  }
  
  // k-mer抽出位置の計算（k文字のDNAに必要な2k個のバイナリビット）
  int start_pos = std::max(0, pos - k*2 + 2);
  int extract_length = std::min(k*2, seq_len - start_pos);
  
  if (extract_length < k*2) {
    return std::string(k, 'A'); // デフォルト
  }
  
  // バイナリからDNAに変換
  std::string dna_kmer = binaryToDNA(seq + start_pos, extract_length);
  
  // k文字に満たない場合は'A'で埋める
  while (dna_kmer.length() < static_cast<size_t>(k)) {
    dna_kmer += 'A';
  }
  
  // k文字を超える場合は切り詰める
  if (dna_kmer.length() > static_cast<size_t>(k)) {
    dna_kmer = dna_kmer.substr(0, k);
  }
  
  return dna_kmer;
}

//================================================================================
void SLFBAdec::exportTransitionProbabilities(const char* filename) {
  printf("# Exporting transition probabilities to %s\n", filename);
  
  std::ofstream file(filename);
  if (!file.is_open()) {
    printf("# Error: Cannot create file %s\n", filename);
    return;
  }
  
  std::vector<std::string> event_names = {"Match", "Insertion", "Deletion", "Substitution"};
  
  // ヘッダー
  file << "# k-mer Transition Probability Table\n";
  file << "# Format: kmer prev_event match_prob ins_prob del_prob sub_prob\n";
  file << "# prev_event: 0=Match, 1=Insertion, 2=Deletion, 3=Substitution\n";
  file << "kmer\tprev_event\tmatch_prob\tins_prob\tdel_prob\tsub_prob\n";
  
  for (const auto& kmer_entry : transition_probs) {
    const std::string& kmer = kmer_entry.first;
    for (const auto& event_entry : kmer_entry.second) {
      int prev_event = event_entry.first;
      const std::vector<double>& probs = event_entry.second;
      
      file << kmer << "\t" << prev_event << "\t";
      file << std::scientific << std::setprecision(6);
      file << probs[0] << "\t" << probs[1] << "\t" << probs[2] << "\t" << probs[3] << "\n";
    }
  }
  
  file.close();
  printf("# Exported %zu k-mers to %s\n", transition_probs.size(), filename);
}

//================================================================================
void SLFBAdec::exportGXNewTable(const char* filename) {
  printf("# Exporting GXNew table to %s\n", filename);
  
  std::ofstream file(filename);
  if (!file.is_open()) {
    printf("# Error: Cannot create file %s\n", filename);
    return;
  }
  
  // ヘッダー
  file << "# GXNew Probability Table: P(y|x,e)\n";
  file << "# Format: Nu2 binary(y) binary(xi) error_state probability\n";
  file << "# error_state: 0=Match, 1=Insertion, 2=Deletion, 3=Substitution\n";
  file << "Nu2\tbinary(y)\tbinary(xi)\terror_state\tprobability\n";
  
  unsigned char *Y = new unsigned char [Nu*2]; // 最大受信長
  unsigned char *X = new unsigned char [Nu];   // 送信コードワード
  
  int count = 0;
  for(int Nu2 = 0; Nu2 <= Nu*2; Nu2++) {
    if(Nu2 < Nu2min || Nu2 > Nu2max) {
      // approximate case
      for(int e = 0; e < 4; e++) {
        file << Nu2 << "\t";
        
        // y=0のバイナリ表現
        if(Nu2 > 0) {
          LongToVect(Y, 0, Nu2);
          for(int i = 0; i < Nu2; i++) file << (int)Y[i];
        } else {
          file << "(empty)";
        }
        file << "\t";
        
        // xi=0のバイナリ表現（コードワード）
        ICB->Get_CW(X, 0);
        for(int i = 0; i < Nu; i++) file << (int)X[i];
        
        file << "\t" << e << "\t" << std::scientific << std::setprecision(6) 
             << GXNew[Nu2][0][0][e] << "\n";
        count++;
      }
    } else {
      // exact case
      long Nu2p = (long)pow(2, Nu2);
      for(long y = 0; y < Nu2p && y < 100; y++) { // 最初の100個のyのみ出力（ファイルサイズ制限）
        // yをバイナリに変換
        LongToVect(Y, y, Nu2);
        
        for(int xi = 0; xi < Q; xi++) {
          // xiのコードワードを取得
          ICB->Get_CW(X, xi);
          
          for(int e = 0; e < 4; e++) {
            file << Nu2 << "\t";
            
            // yのバイナリ表現
            for(int i = 0; i < Nu2; i++) file << (int)Y[i];
            file << "\t";
            
            // xiのバイナリ表現（コードワード）
            for(int i = 0; i < Nu; i++) file << (int)X[i];
            
            file << "\t" << e << "\t" 
                 << std::scientific << std::setprecision(6) 
                 << GXNew[Nu2][y][xi][e] << "\n";
            count++;
          }
        }
      }
    }
  }
  
  delete [] Y;
  delete [] X;
  file.close();
  printf("# Exported %d GXNew entries to %s\n", count, filename);
}

//================================================================================
void SLFBAdec::exportGENewTable(const char* filename) {
  printf("# Exporting GENew table to %s\n", filename);
  
  std::ofstream file(filename);
  if (!file.is_open()) {
    printf("# Error: Cannot create file %s\n", filename);
    return;
  }
  
  // ヘッダー
  file << "# GENew Probability Table: P(e_(i+1)v|...)\n";
  file << "# Format: Nu2 binary(y) binary(xi) error_state probability\n";
  file << "# error_state: 0=Match, 1=Insertion, 2=Deletion, 3=Substitution\n";
  file << "Nu2\tbinary(y)\tbinary(xi)\terror_state\tprobability\n";
  
  unsigned char *Y = new unsigned char [Nu*2]; // 最大受信長
  unsigned char *X = new unsigned char [Nu];   // 送信コードワード
  
  int count = 0;
  for(int Nu2 = 0; Nu2 <= Nu*2; Nu2++) {
    if(Nu2 < Nu2min || Nu2 > Nu2max) {
      // approximate case
      for(int e = 0; e < 4; e++) {
        file << Nu2 << "\t";
        
        // y=0のバイナリ表現
        if(Nu2 > 0) {
          LongToVect(Y, 0, Nu2);
          for(int i = 0; i < Nu2; i++) file << (int)Y[i];
        } else {
          file << "(empty)";
        }
        file << "\t";
        
        // xi=0のバイナリ表現（コードワード）
        ICB->Get_CW(X, 0);
        for(int i = 0; i < Nu; i++) file << (int)X[i];
        
        file << "\t" << e << "\t" << std::scientific << std::setprecision(6) 
             << GENew[Nu2][0][0][e] << "\n";
        count++;
      }
    } else {
      // exact case
      long Nu2p = (long)pow(2, Nu2);
      for(long y = 0; y < Nu2p && y < 100; y++) { // 最初の100個のyのみ出力（ファイルサイズ制限）
        // yをバイナリに変換
        LongToVect(Y, y, Nu2);
        
        for(int xi = 0; xi < Q; xi++) {
          // xiのコードワードを取得
          ICB->Get_CW(X, xi);
          
          for(int e = 0; e < 4; e++) {
            file << Nu2 << "\t";
            
            // yのバイナリ表現
            for(int i = 0; i < Nu2; i++) file << (int)Y[i];
            file << "\t";
            
            // xiのバイナリ表現（コードワード）
            for(int i = 0; i < Nu; i++) file << (int)X[i];
            
            file << "\t" << e << "\t" 
                 << std::scientific << std::setprecision(6) 
                 << GENew[Nu2][y][xi][e] << "\n";
            count++;
          }
        }
      }
    }
  }
  
  delete [] Y;
  delete [] X;
  file.close();
  printf("# Exported %d GENew entries to %s\n", count, filename);
}

//================================================================================
void SLFBAdec::exportAllProbabilityTables(const char* output_dir) {
  printf("# Exporting all probability tables to directory: %s\n", output_dir);
  
  // ディレクトリ作成（簡単な方法）
  std::string mkdir_cmd = "mkdir -p ";
  mkdir_cmd += output_dir;
  system(mkdir_cmd.c_str());
  
  // 各テーブルを個別ファイルに出力
  std::string transition_file = std::string(output_dir) + "/transition_probabilities.txt";
  std::string gxnew_file = std::string(output_dir) + "/GXNew_table.txt";
  std::string genew_file = std::string(output_dir) + "/GENew_table.txt";
  
  exportTransitionProbabilities(transition_file.c_str());
  exportGXNewTable(gxnew_file.c_str());
  exportGENewTable(genew_file.c_str());
  
  // サマリーファイルも作成
  std::string summary_file = std::string(output_dir) + "/probability_tables_summary.txt";
  std::ofstream summary(summary_file);
  if (summary.is_open()) {
    summary << "# SLFBAdec Probability Tables Summary\n";
    summary << "# Generated on: " << __DATE__ << " " << __TIME__ << "\n\n";
    
    summary << "## Parameters\n";
    summary << "Nu (symbol length): " << Nu << "\n";
    summary << "Q (codewords): " << Q << "\n";
    summary << "k-mer length: " << k_mer_length << "\n";
    summary << "Nu2min: " << Nu2min << "\n";
    summary << "Nu2max: " << Nu2max << "\n";
    summary << "Channel parameters: Pi=" << Pi << " Pd=" << Pd << " Ps=" << Ps << " Pt=" << Pt << "\n\n";
    
    summary << "## Generated Files\n";
    summary << "1. transition_probabilities.txt - k-mer dependent transition probabilities\n";
    summary << "2. GXNew_table.txt - Channel output probabilities P(y|x,e)\n";
    summary << "3. GENew_table.txt - Error state transition probabilities P(e_(i+1)v|...)\n\n";
    
    summary << "## k-mer Statistics\n";
    summary << "Total k-mers in transition table: " << transition_probs.size() << "\n";
    
    // k-merごとの統計
    if (!transition_probs.empty()) {
      summary << "\n## Sample k-mers (first 10)\n";
      int sample_count = 0;
      for (const auto& kmer_entry : transition_probs) {
        if (sample_count >= 10) break;
        summary << kmer_entry.first << ": " << kmer_entry.second.size() << " event types\n";
        sample_count++;
      }
    }
    
    summary.close();
    printf("# Created summary file: %s\n", summary_file.c_str());
  }
  
  printf("# All probability tables exported to: %s\n", output_dir);
}

//================================================================================
void SLFBAdec::SetFG(){
  PF = new double * [Ns+1];
  PB = new double * [Ns+1];
  for(int i=0;i<Ns+1;i++){
    PF[i] = new double [Drng];
    PB[i] = new double [Drng];
  } // for i
  PU = new double * [Ns];
  PD = new double * [Ns];
  PI = new double * [Ns];
  PO = new double * [Ns];
  for(int i=0;i<Ns;i++){
    PU[i] = new double [Q];
    PD[i] = new double [Q];
    PI[i] = new double [Q];
    PO[i] = new double [Q];
  } // for i
  Yin = new unsigned char [Nb+Dmax];
}

//================================================================================
void SLFBAdec::DelFG(){
  for(int i=0;i<Ns+1;i++){
    delete [] PF[i];
    delete [] PB[i];
  } // for i
  delete [] PF;
  delete [] PB;
  for(int i=0;i<Ns;i++){
    delete [] PU[i];
    delete [] PD[i];
    delete [] PI[i];
    delete [] PO[i];
  } // for i
  delete [] PU;
  delete [] PD;
  delete [] PI;
  delete [] PO;
  delete [] Yin;
}

//================================================================================
void SLFBAdec::InitFG(const unsigned char *RW, const double **Pin, int Nb2){
  ClearMat(PF,Ns+1,Drng);
  ClearMat(PB,Ns+1,Drng);
  ClearMat(PU,Ns,Q);
  ClearMat(PD,Ns,Q);
  ClearMat(PO,Ns,Q);
  CopyMat(PI,Pin,Ns,Q);
  PF[0 ][0     -Dmin] = 1.0;
  PB[Ns][Nb2-Nb-Dmin] = 1.0;
  memcpy(Yin,RW,sizeof(unsigned char)*Nb2);
}

//================================================================================
void SLFBAdec::CalcPU(int idx){
  // p(u'|u) 符号化確率
  assert(idx>=0 && idx<Ns);
  int Vin,Vup;
  for(Vup=0;Vup<Q;Vup++){
    PU[idx][Vup] = 0.0;
    for(Vin=0;Vin<Q;Vin++){
      PU[idx][Vup] += (PI[idx][Vin] * ECM->GetPxy(Vin,Vup));
    } // for Vin
  } // for Vup
  normalize(PU[idx], Q);
}

//================================================================================
void SLFBAdec::CalcPO(int idx){
  assert(idx>=0 && idx<Ns);
  int Vout,Vdn;
  for(Vout=0;Vout<Q;Vout++){
    PO[idx][Vout] = 0.0;
    for(Vdn=0;Vdn<Q;Vdn++){
      PO[idx][Vout] += (PD[idx][Vdn] * ECM->GetPxy(Vout,Vdn));
    } // for Vdn
  } // for Vout
  normalize(PO[idx], Q);
}

//================================================================================
void SLFBAdec::CalcPF(int idx, int Nb2){
  // Forward確率
  assert(idx>=0 && idx<Ns);
  int  d0,d1;
  int  xi;
  int  Nu2, iL, iR;
  //long x;
  long y;
  double s,ss;
  //unsigned char *X = new unsigned char [Nu];
  for(d1=Dmin;d1<=Dmax;d1++){
    PF[idx+1][d1-Dmin] = 0.0;
    for(d0=Dmin;d0<=Dmax;d0++){
      Nu2 = Nu+d1-d0;
      iL  = idx    *Nu + d0;
      iR  = (idx+1)*Nu + d1 - 1;
      if( (Nu2<0) || (Nu2>2*Nu) || (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
      if( Nu2>=Nu2min && Nu2<=Nu2max ){
	y = VectToLong( &Yin[iL], Nu2 );
	s = 0.0;
	for(xi=0;xi<Q;xi++){
	  //ICB->Get_CW(X,xi);
	  //x  = VectToLong(X,Nu);
	  ss = PU[idx][xi] * GetGX(Nu2,y,xi);
	  s += ss;
	  //(dbg)
	  //printf("%e d1=%d d0=%d Nu2=%d iL=%d iR=%d y=%ld xi=%d x=%ld\n",ss,d1,d0,Nu2,iL,iR,y,xi,x);
	} // for xi
      } else {
	s = GetGX(Nu2,0,0);
      } // if Nu2
      PF[idx+1][d1-Dmin] += (s*PF[idx][d0-Dmin]);
    } // for d0
  } // for d1
  normalize(PF[idx+1],Drng);
  //delete [] X;
}

//================================================================================
void SLFBAdec::CalcPB(int idx, int Nb2){
  // Backward確率
  assert(idx>=0 && idx<Ns);
  int  d0,d1;
  int  xi;
  int  Nu2, iL, iR;
  //long x;
  long y;
  double s,ss;
  //unsigned char *X = new unsigned char [Nu];
  for(d0=Dmin;d0<=Dmax;d0++){
    PB[idx][d0-Dmin] = 0.0;
    for(d1=Dmin;d1<=Dmax;d1++){
      Nu2 = Nu+d1-d0;
      iL  = idx    *Nu + d0;
      iR  = (idx+1)*Nu + d1 - 1;
      if( (Nu2<0) || (Nu2>2*Nu) || (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
      if( Nu2>=Nu2min && Nu2<=Nu2max ){
	y = VectToLong( &Yin[iL], Nu2 );
	s = 0.0;
	for(xi=0;xi<Q;xi++){
	  //ICB->Get_CW(X,xi);
	  //x  = VectToLong(X,Nu);
	  ss = PU[idx][xi] * GetGX(Nu2,y,xi);
	  s += ss;
	  //(dbg)
	  //printf("%e d1=%d d0=%d Nu2=%d iL=%d iR=%d y=%ld xi=%d x=%ld\n",ss,d1,d0,Nu2,iL,iR,y,xi,x);
	} // for xi
      } else {
	s = GetGX(Nu2,0,0);
      } // if Nu2
      PB[idx][d0-Dmin] += (s*PB[idx+1][d1-Dmin]);
    } // for d1
  } // for d0
  normalize(PB[idx],Drng);
  //delete [] X;
}

//================================================================================
void SLFBAdec::CalcPD(int idx, int Nb2){
  //　事後確率
  assert(idx>=0 && idx<Ns);
  int  d0,d1;
  int  xi;
  int  Nu2, iL, iR;
  //long x;
  long y;
  double s,ss;
  //unsigned char *X = new unsigned char [Nu];
  for(xi=0;xi<Q;xi++){
    PD[idx][xi] = 0.0;
    //ICB->Get_CW(X,xi);
    //x  = VectToLong(X,Nu);
    for(d1=Dmin;d1<=Dmax;d1++){
      s = 0.0;
      for(d0=Dmin;d0<=Dmax;d0++){
	Nu2 = Nu+d1-d0;
	iL  = idx    *Nu + d0;
	iR  = (idx+1)*Nu + d1 - 1;
	if( (Nu2<0) || (Nu2>2*Nu) || (iL<0) || (iL>=Nb2) || (iR<-1) || (iR>=Nb2) ) continue;
	y = VectToLong( &Yin[iL], Nu2 );
	ss = PF[idx][d0-Dmin] * GetGX(Nu2,y,xi);
	s += ss;
	//(dbg)
	//printf("%e d1=%d d0=%d Nu2=%d iL=%d iR=%d y=%ld xi=%d x=%ld\n",ss,d1,d0,Nu2,iL,iR,y,xi,x);
      } // for d0
      PD[idx][xi] += (s*PB[idx+1][d1-Dmin]);
    } // for d1
  } // for xi
  normalize(PD[idx],Q);
  //delete [] X;
}

//================================================================================
//================================================================================
//================================================================================

//================================================================================
SLFBAdec::SLFBAdec(class InnerCodebook *_ICB, class ChannelMatrix *_ECM, class IDSchannel *_CH){
  ICB  = _ICB;
  ECM  = _ECM;
  CH   = _CH;
  Nu   = ICB->Get_Nu();
  Q    = ICB->Get_numCW();
  Nb   = CH->GetN();
  Dmax = CH->GetDmax();
  Dmin = CH->GetDmin();
  Pi   = CH->GetPi();
  Pd   = CH->GetPd();
  Ps   = CH->GetPs();
  Pt   = 1.0-Pi-Pd;
  Nu2p = (long) pow(2,Nu);
  Ns   = Nb/Nu;
  Drng = Dmax-Dmin+1;
  Nu2max = Nu + (int)ceil( (double)Nu*Pi ) + 2;
  Nu2min = Nu - (int)ceil( (double)Nu*Pd ) - 2;
  Nu2max = min(Nu2max, Nu*2);
  Nu2min = max(Nu2min, 0   );
  k_mer_length = 4; // デフォルトk=4、必要に応じて設定可能にする
  printf("# SLFBAdec: Ns=%d (Nu;Nu2min,Nu2max)=(%d(%ld);%d,%d) Nb=%d Q=%d\n",Ns,Nu,Nu2p,Nu2min,Nu2max,Nb,Q);
  printf("# SLFBAdec: (Pi,Pd,Ps,Pt)=(%e,%e,%e,%e) (Dmin,Dmax:Drng)=(%d,%d:%d)\n",
	 Pi,Pd,Ps,Pt,Dmin,Dmax,Drng);
  printf("# SLFBAdec: k-mer length=%d\n",k_mer_length);
  printf("# SLFBAdec: loading transition probabilities\n");
  loadTransitionProbabilities(k_mer_length);
  assert(Nb%Nu==0);
  assert(Nu<=NuMax);
  assert(ECM->GetM()==Q && ECM->GetN()==Q);
  assert(Pt>0.0 && Pt<=1.0);
  //----- set tables
  printf("# SLFBAdec: generating GD & GX\n");
  SetGD();
  SetGX();
  printf("# SLFBAdec: generating GXNew\n");
  SetGXNew();
  printf("# SLFBAdec: generating GENew\n");
  SetGENew();
  printf("# SLFBAdec: generating factor graph\n");
  SetFG();
}

//================================================================================
SLFBAdec::~SLFBAdec(){
  DelGD();
  DelGX();
  DelGXNew();
  DelGENew();
  DelFG();
  printf("# SLFBAdec: deleted\n");
}

//================================================================================
void SLFBAdec::Decode(double **Pout, const unsigned char *RW, int Nb2, const int *dbgIW, const double **Pin){
  int idx;
  assert( Nb2>=Nb+Dmin && Nb2<=Nb+Dmax );
  InitFG(RW,Pin,Nb2);
  for(idx=0;   idx<Ns;idx++) CalcPU(idx);
  for(idx=0;   idx<Ns;idx++) CalcPF(idx,Nb2);
  for(idx=Ns-1;idx>=0;idx--) CalcPB(idx,Nb2);
  for(idx=0;   idx<Ns;idx++) CalcPD(idx,Nb2);
  for(idx=0;   idx<Ns;idx++) CalcPO(idx);
  CopyMat(Pout,(const double**)PO,Ns,Q);
  
  // PrintNode(dbgIW);
  // PrintNode(0,   dbgIW[0]);
  // PrintNode(1,   dbgIW[1]);
  // PrintNode(2,   dbgIW[2]);
  // PrintNode(Ns-4,dbgIW[Ns-2]);
  // PrintNode(Ns-3,dbgIW[Ns-2]);
  // PrintNode(Ns-2,dbgIW[Ns-2]);
  // PrintNode(Ns-1,dbgIW[Ns-1]);
}

//================================================================================
void SLFBAdec::Decode(double **Pout, const unsigned char *RW, int Nb2, const int *dbgIW){
  double **Pin = new double * [Ns];
  for(int i=0;i<Ns;i++){
    Pin[i] = new double [Q];
    for(int x=0;x<Q;x++) Pin[i][x] = (double)1.0/Q;
  } // for i
  Decode(Pout,RW,Nb2,dbgIW,(const double **)Pin);
  for(int i=0;i<Ns;i++) delete [] Pin[i];
  delete [] Pin;
}

//================================================================================
void SLFBAdec::PrintNode(int idx, int iw){
  assert(idx>=0 && idx<Ns);
  printf("[%04d] %02d\n",idx,iw);
  PrintVect(PF[idx  ],Drng,"PF: ","\n");
  PrintVect(PB[idx+1],Drng,"PB: ","\n");
  PrintVect(PI[idx  ],Q,   "PI: ","\n");
  PrintVect(PU[idx  ],Q,   "PU: ","\n");
  PrintVect(PD[idx  ],Q,   "PD: ","\n");
  PrintVect(PO[idx  ],Q,   "PO: ","\n");
}

//================================================================================
void SLFBAdec::PrintNode(const int *dbgIW){
  for(int i=0;i<Ns;i++) PrintNode(i,dbgIW[i]);
}

