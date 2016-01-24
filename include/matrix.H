// 行列クラスライブラリ


#ifndef MATRIX_H
#define MATRIX_H


#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>


/******** クラス一覧 ********/

template<class T> class matrix;  // T値行列
typedef matrix<double> Matrix;   // 実数版


/******** クラスI/F ********/

// 行列クラス
template<class T> class matrix{
public:
  explicit matrix(size_t r=0, size_t c=1);  // 指定サイズ
  matrix(T x, T y, T z);                    // 3Dベクトル

  size_t nrow(void) const;
  size_t ncol(void) const;
  size_t dim (void) const;

  T  operator()(size_t i, size_t j) const;  // 参照 t=M(i,j)
  T& operator()(size_t i, size_t j);        // 設定 M(i,j)=t
  T  operator()(size_t i) const;  // t=V(i)
  T& operator()(size_t i);        // V(i)=t

private:
  size_t row,col;      // サイズ
  std::vector<T> val;  // データ
};


/******** ユーティリティ関数I/F ********/

template<class T> matrix<T>& operator+=(matrix<T>& m, T t);  // インクリメンタル演算..
template<class T> matrix<T>& operator-=(matrix<T>& m, T t);
template<class T> matrix<T>& operator*=(matrix<T>& m, T t);
template<class T> matrix<T>& operator/=(matrix<T>& m, T t);
template<class T> matrix<T>& operator+=(matrix<T>& m1, const matrix<T>& m2);
template<class T> matrix<T>& operator-=(matrix<T>& m1, const matrix<T>& m2);
//template<class T> matrix<T>& operator*=(matrix<T>& m1, const matrix<T>& m2);
//template<class T> matrix<T>& operator/=(matrix<T>& m1, const matrix<T>& m2);

template<class T> matrix<T> operator-(const matrix<T>& m     );  // 単項マイナス
template<class T> matrix<T> operator+(const matrix<T>& m, T t);  // スカラー演算..
template<class T> matrix<T> operator-(const matrix<T>& m, T t);
template<class T> matrix<T> operator*(const matrix<T>& m, T t);
template<class T> matrix<T> operator/(const matrix<T>& m, T t);
template<class T> matrix<T> operator+(T t, const matrix<T>& m);
template<class T> matrix<T> operator-(T t, const matrix<T>& m);
template<class T> matrix<T> operator*(T t, const matrix<T>& m);
template<class T> matrix<T> operator/(T t, const matrix<T>& m);
template<class T> matrix<T> operator+(const matrix<T>& m1, const matrix<T>& m2);  // 行列同士..
template<class T> matrix<T> operator-(const matrix<T>& m1, const matrix<T>& m2);
template<class T> matrix<T> operator*(const matrix<T>& m1, const matrix<T>& m2);
//template<class T> matrix<T> operator/(const matrix<T>& m1, const matrix<T>& m2);

template<class T> T          norm  (const matrix<T>& m);                        // ノルム
template<class T> T          inner (const matrix<T>& m1, const matrix<T>& m2);  // 内積
template<class T> matrix<T>  outer (const matrix<T>& v1, const matrix<T>& v2);  // 外積(3D限定)
template<class T> matrix<T>  inv   (const matrix<T>& m);                        // 逆行列 M^-1
template<class T> matrix<T>  pow   (const matrix<T>& m, int n);                 // べき乗 M^n (n<0も可)
template<class T> matrix<T>  tp    (const matrix<T>& m);                        // 転置   M^T
template<class T> matrix<T>  orth  (const matrix<T>& m);                        // Schmidt正規直交化
template<class T> matrix<T>  eigen (const matrix<T>& m);                        // 固有ベクトル(を並べた直交行列)
template<class T> matrix<T>  getvec(const matrix<T>& m, size_t j);              // 列ベクトル取得
template<class T> void       setvec(      matrix<T>& m, size_t j, const matrix<T>& v);  // 〃 設定
template<class T> matrix<T>  hcat  (const matrix<T>& m1, const matrix<T>& m2);  // 横に連結
template<class T> matrix<T>  vcat  (const matrix<T>& m1, const matrix<T>& m2);  // 縦に連結

template<class T> std::istream& operator>>(std::istream&,       matrix<T>& m);  // 入出力
template<class T> std::ostream& operator<<(std::ostream&, const matrix<T>& m);


/******** インライン実装 ********/

// クラスメソッド
template<class T>        matrix<T>::matrix(size_t r, size_t c): row(r), col(c), val(r*c) {}
template<class T>        matrix<T>::matrix(T x, T y, T z): row(3), col(1), val(3) {  val[0]=x;  val[1]=y;  val[2]=z;  }
template<class T> size_t matrix<T>::nrow(void) const {  return row;  }
template<class T> size_t matrix<T>::ncol(void) const {  return col;  }
template<class T> size_t matrix<T>::dim (void) const {  return row*col;  }
template<class T> T      matrix<T>::operator()(size_t i, size_t j) const {  return val[i*col+j];  }
template<class T> T      matrix<T>::operator()(size_t i          ) const {  return val[i];        }
template<class T> T&     matrix<T>::operator()(size_t i, size_t j)       {  return val[i*col+j];  }
template<class T> T&     matrix<T>::operator()(size_t i          )       {  return val[i];        }

// インクリメンタル演算
template<class T> matrix<T>& operator+=(matrix<T>& m, T t){
  if(m.nrow()!=m.ncol()) throw std::domain_error("operator+=(M,t): not square");
  for(size_t i=0; i<m.nrow(); i++) m(i,i)+=t;
  return m;
}
template<class T> matrix<T>& operator-=(matrix<T>& m, T t){  m+=(-t);  return m;  }
template<class T> matrix<T>& operator*=(matrix<T>& m, T t){  for(size_t k=0; k<m.dim(); k++) m(k)*=t;  return m;  }
template<class T> matrix<T>& operator/=(matrix<T>& m, T t){  m*=(1/t);  return m;  }
template<class T> matrix<T>& operator+=(matrix<T>& m, const matrix<T>& dm){
  if(!(m.nrow()==dm.nrow() and m.ncol()==dm.ncol())) throw std::domain_error("operator+=(M,M): sizes mismatch");
  for(size_t k=0; k<m.dim(); k++) m(k)+=dm(k);
  return m;
}
template<class T> matrix<T>& operator-=(matrix<T>& m, const matrix<T>& dm){
  if(!(m.nrow()==dm.nrow() and m.ncol()==dm.ncol())) throw std::domain_error("operator-=(M,M): sizes mismatch");
  for(size_t k=0; k<m.dim(); k++) m(k)-=dm(k);
  return m;
}

// スカラー演算
template<class T> matrix<T> operator-(const matrix<T>& m0     ){  matrix<T> m=m0;  m*=-T(1); return m;  }  // 単項マイナス
template<class T> matrix<T> operator+(const matrix<T>& m0, T t){  matrix<T> m=m0;  m+=t;  return m;  }  // 行列 vs. スカラー
template<class T> matrix<T> operator-(const matrix<T>& m0, T t){  matrix<T> m=m0;  m-=t;  return m;  }
template<class T> matrix<T> operator*(const matrix<T>& m0, T t){  matrix<T> m=m0;  m*=t;  return m;  }
template<class T> matrix<T> operator/(const matrix<T>& m0, T t){  matrix<T> m=m0;  m/=t;  return m;  }
template<class T> matrix<T> operator+(T t, const matrix<T>& m0){  matrix<T> m=m0;  m+=t;  return m;  }  // スカラー vs. 行列
template<class T> matrix<T> operator-(T t, const matrix<T>& m0){  matrix<T> m=-m0; m+=t;  return m;  }
template<class T> matrix<T> operator*(T t, const matrix<T>& m0){  matrix<T> m=m0;  m*=t;  return m;  }
template<class T> matrix<T> operator/(T t, const matrix<T>& m0){  matrix<T> m=inv(m0);  m*=t;  return m;  }

// 行列同士の演算
template<class T> matrix<T> operator+(const matrix<T>& m1, const matrix<T>& m2){  matrix<T> m=m1;  m+=m2;  return m;  }
template<class T> matrix<T> operator-(const matrix<T>& m1, const matrix<T>& m2){  matrix<T> m=m1;  m-=m2;  return m;  }
template<class T> matrix<T> operator*(const matrix<T>& m1, const matrix<T>& m2){
  if(m1.ncol()!=m2.nrow()) throw std::domain_error("operator*(M,M): sizes mismatch");
  matrix<T> m(m1.nrow(),m2.ncol());
  for(size_t i=0; i<m.nrow(); i++)
    for(size_t j=0; j<m.ncol(); j++){
      T s=0;
      for(size_t k=0; k<m1.ncol(); k++) s+=m1(i,k)*m2(k,j);
      m(i,j)=s;
    }
  return m;
}


// ノルム
template<class T> T norm (const matrix<T>& m){
  return sqrt(inner(m,m));
}


// 内積
template<class T> T inner(const matrix<T>& m1, const matrix<T>& m2){
  if( !(m1.nrow()==m2.nrow() and m1.ncol()==m2.ncol()) ) throw std::domain_error("inner(M,M): sizes mismatch");
  T t=0;
  for(size_t k=0; k<m1.dim(); k++) t+=m1(k)*m2(k);
  return t;
}


// 外積(3D限定)
template<class T> matrix<T> outer(const matrix<T>& v, const matrix<T>& w){
  if( !(v.nrow()==3 and w.ncol()==1 and v.nrow()==3 and v.ncol()==1) ) throw std::domain_error("outer(V,V): not 3-vectors");
  return matrix<T>(v(1)*w(2)-w(1)*v(2), v(2)*w(0)-w(2)*v(0), v(0)*w(1)-w(0)*v(1));
}


// 逆行列 M^-1
template<class T> matrix<T> inv(const matrix<T>& m){
  if( m.nrow()!=m.ncol() ) throw std::domain_error("inv(M): not square");
  matrix<T> m1=m;
  matrix<T> m2=matrix<T>(m.nrow(),m.nrow())+T(1);  // m^-1 の入れ物

  for(size_t k=0; k<m.nrow(); k++){
    // pivote: k行 <--> i0行
    size_t i0=k;
    for(size_t i=k+1; i<m1.nrow(); i++)
      if(fabs(m1(i,k))>fabs(m1(i0,k))) i0=i;  // i0 = argmax_{i: i>k} |m1(i,k)|
    if(i0!=k){
      for(size_t j=k; j<m1.nrow(); j++) std::swap(m1(k,j),m1(i0,j));
      for(size_t j=0; j<m2.nrow(); j++) std::swap(m2(k,j),m2(i0,j));
    }

    // sweep: m1(i,k) --> δ(i,k)
    T t;
    t=m1(k,k);  // ゼロ除算注意
    for(size_t j=k; j<m1.ncol(); j++) m1(k,j)/=t;
    for(size_t j=0; j<m2.ncol(); j++) m2(k,j)/=t;
    for(size_t i=0; i<m1.nrow(); i++){
      if(i==k) continue;
      t=m1(i,k);
      for(size_t j=k; j<m1.ncol(); j++) m1(i,j)-=m1(k,j)*t;
      for(size_t j=0; j<m2.ncol(); j++) m2(i,j)-=m2(k,j)*t;
    }
  }

  return m2;
}


// べき乗 M^n (n<=0も可)
template<class T> matrix<T> pow(const matrix<T>& m, int n){
  if(m.nrow()!=m.ncol()) throw std::domain_error("pow(M,n): not square");
  if(n<0) return inv(pow(m,-n));
  matrix<T> m1=m, mn=matrix<T>(m.nrow(),m.ncol())+T(1);
  for( ; n>0; n/=2,m1=m1*m1)
    if(n%2==1) mn=mn*m1;
  return mn;
}


// 転置 M^T
template<class T> matrix<T> tp(const matrix<T>& m0){
  matrix<T> m1(m0.ncol(),m0.nrow());
  for(size_t i=0; i<m0.nrow(); i++)
    for(size_t j=0; j<m0.ncol(); j++)
      m1(j,i)=m0(i,j);
  return m1;
}


// Schmidt正規直交化
template<class T> matrix<T> orth(const matrix<T>& m0){
  if(m0.nrow()<m0.ncol()) throw std::domain_error("orth(M): nrow<ncol");
  matrix<T> m=m0, v, v1;
  for(size_t j=0; j<m.ncol(); j++){
    v=getvec(m,j);
    for(size_t j1=0; j1<j; j1++){
      v1=getvec(m,j1);
      v-=inner(v,v1)*v1;
    }
    v/=norm(v);
    setvec(m,j,v);
  }
  return m;
}


// 固有ベクトル(を並べた直交行列): U=eigen(M) => U^T M U = diagonal
template<class T> matrix<T> eigen(const matrix<T>& m0){
  if(m0.nrow()!=m0.ncol()) throw std::domain_error("eigen(M): not square");
  size_t n=m0.nrow();
  matrix<T> m=m0+tp(m0), u=matrix<T>(n,n)+T(1);

  for(int extra=m0.dim(); extra; ){
    // 絶対値最大の成分(i0,j0)
    size_t i0=0, j0=1;
    for(size_t i=0; i<n; i++)
      for(size_t j=0; j<n; j++)
	if( i!=j and fabs(m(i0,j0))<=fabs(m(i,j)) ) { i0=i; j0=j; }

    // Givens回転角θ
    T th=atan2(2*m(i0,j0),m(i0,i0)-m(j0,j0))/2;
    matrix<T> rot(n,n);
    rot+=T(1);
    rot(i0,i0)=cos(th);  rot(i0,j0)=-sin(th);
    rot(j0,i0)=sin(th);  rot(j0,j0)= cos(th);

    // 相似変換＆収束判定
    T s0=0;  for(size_t k=0; k<n; k++) s0+=m(k,k)*m(k,k);
    u=u*rot;
    m=tp(u)*m0*u;
    T s1=0;  for(size_t k=0; k<n; k++) s1+=m(k,k)*m(k,k);
    if(!(s0<s1)) extra--;  // Sは単調増加だが、誤差を見て多めにループする
  }

  return u;  // 固有ベクトルを並べた直交行列
}


// 列ベクトル取得・設定
template<class T> matrix<T> getvec(const matrix<T>& m, size_t j){
  if(!(j<m.ncol())) throw std::domain_error("getvec(M,j): bad column");
  matrix<T> v(m.nrow());
  for(size_t i=0; i<v.nrow(); i++) v(i)=m(i,j);
  return v;
}
template<class T> void setvec(matrix<T>& m, size_t j, const matrix<T>& v){
  if(!(j<m.ncol() and m.nrow()==v.nrow() and v.ncol()==1)) throw std::domain_error("setvec(M,j,V): bad column or size");
  for(size_t i=0; i<m.nrow(); i++) m(i,j)=v(i);
}


// 横・縦に連結
template<class T> matrix<T> hcat(const matrix<T>& m1, const matrix<T>& m2){
  if(m1.nrow()!=m2.nrow()) throw std::domain_error("hcat(M,M): sizes mismatch");
  matrix<T> m(m1.nrow(),m1.ncol()+m2.ncol());
  for(size_t i=0; i<m.nrow(); i++)
    for(size_t j=0; j<m.ncol(); j++)
      m(i,j) = (j<m1.ncol() ? m1(i,j) : m2(i,j-m1.ncol()));
  return m;
}
template<class T> matrix<T> vcat(const matrix<T>& m1, const matrix<T>& m2){
  if(m1.ncol()!=m2.ncol()) throw std::domain_error("vcat(M,M): sizes mismatch");
  matrix<T> m(m1.nrow()+m2.nrow(),m1.ncol());
  for(size_t i=0; i<m.nrow(); i++)
    for(size_t j=0; j<m.ncol(); j++)
      m(i,j) = (i<m1.nrow() ? m1(i,j) : m2(i-m1.nrow(),j));
  return m;
}


// 入出力
template<class T> std::istream& operator>>(std::istream& str, matrix<T>& m){
  size_t r,c;
  str >> r >> c;
  m=matrix<T>(r,c);
  for(size_t k=0; k<r*c; k++) str >> m(k);
  return str;
}
template<class T> std::ostream& operator<<(std::ostream& str, const matrix<T>& m){
  size_t r=m.nrow(), c=m.ncol();
  str << r << '\t' << c << '\t';
  for(size_t k=0; k<r*c; k++) str << m(k) << '\t';
  return str;
}


#endif //MATRIX_H
