// スパース行列ライブラリ


#ifndef SPARSE_H
#define SPARSE_H


#include <iostream>
#include <vector>
#include <map>
#include "matrix.H"
namespace nmlib{


/******** クラス定義 ********/

// スパース行列クラス
class Sparse{
public:
  explicit Sparse(size_t r=0, size_t c=0);

  size_t nrow(void) const;  // 行数
  size_t ncol(void) const;  // 列数
  double  operator()(size_t i, size_t j) const;  // 成分参照 s=A(i,j)
  double& operator()(size_t i, size_t j);        // 成分設定 A(i,j)=s

  const std::map<size_t,double>& columns(size_t i) const{ return val[i]; }  // 行iの列j=>値Aij
        std::map<size_t,double>& columns(size_t i)      { return val[i]; }  // 〃(非const版)

private:
  size_t row, col;  // サイズ
  std::vector<std::map<size_t,double> > val;  // 成分 val[i][j]  CRS(Condensed Row Storage)形式
};


// 計算条件
class SparseConf{
public:
  SparseConf(void): tol(1.e-8), loop(100), x0(), lu(), verb(false) {}
  void init_ilu0(const Sparse& a);  // 前処理行列の設定(ILU0)

  double  tol;   // 許容誤差 |Ax-b|/|b|
  int     loop;  // 反復回数上限
  Matrix  x0;    // 初期近似解
  Sparse  lu;    // 前処理行列
  bool    verb;  // 途中経過をcerrに表示
};


/******** ユーティリティ一覧 ********/

// インクリメンタル演算
Sparse& operator+=(Sparse& a, double s);   // A=A+sI
Sparse& operator-=(Sparse& a, double s);   // A=A-sI
Sparse& operator*=(Sparse& a, double s);   // A=A*s
Sparse& operator/=(Sparse& a, double s);   // A=A/s
Sparse& operator+=(Sparse& a, const Sparse& b);  // A=A+B
Sparse& operator-=(Sparse& a, const Sparse& b);  // A=A-B
//Sparse& operator*=(Sparse& a, const Sparse& b);
//Sparse& operator/=(Sparse& a, const Sparse& b);

// スカラー演算
Sparse operator-(const Sparse& a);            // -A (単項)
Sparse operator+(const Sparse& a, double s);  // A+sI
Sparse operator-(const Sparse& a, double s);  // A-sI
Sparse operator*(const Sparse& a, double s);  // A*s
Sparse operator/(const Sparse& a, double s);  // A/s
Sparse operator+(double s, const Sparse& a);  // sI+A
Sparse operator-(double s, const Sparse& a);  // sI-A
Sparse operator*(double s, const Sparse& a);  // s*A
//Sparse operator/(double s, const Sparse& a);  // s*inv(A)

// 行列同士の演算
Sparse operator+(const Sparse& a, const Sparse& b);  // A+B
Sparse operator-(const Sparse& a, const Sparse& b);  // A-B
Sparse operator*(const Sparse& a, const Sparse& b);  // A*B
Matrix operator*(const Sparse& a, const Matrix& b);  // A*b (Matrix版)

// ユーティリティ
Sparse tp   (const Sparse& a);  // A^T
Matrix dense(const Sparse& a);  // 密行列変換
Matrix tpab (const Sparse& a, const Matrix& b);  // A^T * b

// 入出力
std::istream& operator>>(std::istream& s,       Sparse& a);  // テキストI/O
std::ostream& operator<<(std::ostream& s, const Sparse& a);  // 〃
std::istream& read      (std::istream& s,       Sparse& a);  // バイナリI/O
std::ostream& write     (std::ostream& s, const Sparse& a);  // 〃

// スパースソルバ: Ax=b を解く
Matrix solve_cg  (const Sparse& a, const Matrix& b, const SparseConf& cf=SparseConf());  // CG法 (A:対称正定値)
Matrix solve_bcg (const Sparse& a, const Matrix& b, const SparseConf& cf=SparseConf());  // BiCG法
Matrix solve_pbcg(const Sparse& a, const Matrix& b, const SparseConf& cf=SparseConf());  // Preconditioned BiCG

// 前処理
void   ilu0(Sparse& a);  // Aの不完全LU分解ILU(0)を上書き算出
Matrix lux   (const Sparse& lu, const Matrix& x);  // LUx
Matrix sluxb (const Sparse& lu, const Matrix& b);  // (LU)^-1 b (solve LUx=b)
Matrix slutxb(const Sparse& lu, const Matrix& b);  // (LU)^-T b (solve LU^T x=b)


}  //namespace nmlib
#endif //SPARSE_H
