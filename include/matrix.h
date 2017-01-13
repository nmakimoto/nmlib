// Matrix class template library


#ifndef MATRIX_H
#define MATRIX_H


#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>
namespace nmlib{


/******** Class I/F ********/

template<class T> class matrix;  // T-valued matrix
typedef matrix<double> Matrix;   // real matrix (alias)


template<class T> class matrix{
public:
  explicit matrix(void);                      // empty matrix
  explicit matrix(size_t r, size_t c=1);      // rxc zero matrix
  explicit matrix(const std::vector<T>& v);   // nx1 val=v
  explicit matrix(const T *v0, const T *vn);  // nx1 val=[v0..vn)
  explicit matrix(T x, T y, T z);             // 3-vector
  operator const std::vector  <T>& (void) const;
  operator matrix<std::complex<T> >(void) const;

  size_t nrow(void) const;  // number of rows
  size_t ncol(void) const;  // number of columns
  size_t dim (void) const;  // number of components
  matrix<T>& resize(size_t r, size_t c);
  matrix<T>& fill  (T z);

  T  operator()(size_t i, size_t j) const;  // get t=M(i,j)
  T& operator()(size_t i, size_t j);        // set M(i,j)=t
  T  operator()(size_t i) const;  // get t=V(i)
  T& operator()(size_t i);        // set V(i)=t

private:
  size_t row,col;      // dimensions
  std::vector<T> val;  // components
};


/******** Utilitiy I/F ********/

// Basic operations (incremental)
template<class T> matrix<T>& operator+=(matrix<T>& m, T t);  // M+=sI
template<class T> matrix<T>& operator-=(matrix<T>& m, T t);  // M-=sI
template<class T> matrix<T>& operator*=(matrix<T>& m, T t);  // M*=s
template<class T> matrix<T>& operator/=(matrix<T>& m, T t);  // M/=s
template<class T> matrix<T>& operator+=(matrix<T>& m1, const matrix<T>& m2);  // M1+=M2
template<class T> matrix<T>& operator-=(matrix<T>& m1, const matrix<T>& m2);  // M1-=M2
//template<class T> matrix<T>& operator*=(matrix<T>& m1, const matrix<T>& m2);  // M1*=M2 (ambiguous)
//template<class T> matrix<T>& operator/=(matrix<T>& m1, const matrix<T>& m2);  // M1/=M2 (ambiguous)

// Basic operations
template<class T> matrix<T> operator-(const matrix<T>& m     );  // -M (unary)
template<class T> matrix<T> operator+(const matrix<T>& m, T t);  // M+tI
template<class T> matrix<T> operator-(const matrix<T>& m, T t);  // M-tI
template<class T> matrix<T> operator*(const matrix<T>& m, T t);  // M*t
template<class T> matrix<T> operator/(const matrix<T>& m, T t);  // M/t
template<class T> matrix<T> operator+(T t, const matrix<T>& m);  // tI+M
template<class T> matrix<T> operator-(T t, const matrix<T>& m);  // tI-M
template<class T> matrix<T> operator*(T t, const matrix<T>& m);  // t*M
template<class T> matrix<T> operator/(T t, const matrix<T>& m);  // t*M^-1
template<class T> matrix<T> operator+(const matrix<T>& m1, const matrix<T>& m2);  // M1+M2
template<class T> matrix<T> operator-(const matrix<T>& m1, const matrix<T>& m2);  // M1-M2
template<class T> matrix<T> operator*(const matrix<T>& m1, const matrix<T>& m2);  // M1*M2
//template<class T> matrix<T> operator/(const matrix<T>& m1, const matrix<T>& m2);  // M1/M2 (ambiguous)

// Miscelanaous operations
template<class T> T          norm  (const matrix<T>& m);                        // norm |M|
template<class T> T          norm  (const matrix<std::complex<T> >& m);         // norm |M| (complex M)
template<class T> T          inner (const matrix<T>& m1, const matrix<T>& m2);  // inner product <M1,M2>
template<class T> matrix<T>  inv   (const matrix<T>& m);                        // inverse M^-1
template<class T> matrix<T>  pow   (const matrix<T>& m, int n);                 // power M^n (n<0 is ok)
template<class T> matrix<T>  tp    (const matrix<T>& m);                        // transpose M^T
template<class T> matrix<T>  orth  (const matrix<T>& m);                        // Gram-Schmidt orthonormalization of M
template<class T> matrix<T>  eigen (const matrix<T>& m, double tol=0, int iter=-1);  // eigen vectors of symmetric M
template<class T> matrix<T>  getvec(const matrix<T>& m, size_t j);              // get column vector Mj
template<class T> void       setvec(      matrix<T>& m, size_t j, const matrix<T>& v);  // set column vector Mj to V
template<class T> matrix<T>  getsub(const matrix<T>& m, size_t i0, size_t j0, size_t r, size_t c);  // get submatrix
template<class T> void       setsub(      matrix<T>& m, size_t i0, size_t j0, const matrix<T> m1);  // set submatrix
template<class T> matrix<T>  getdiag(const matrix<T>& m);
template<class T> void       setdiag(      matrix<T>& m, const matrix<T>& v);
template<class T> matrix<T>  hcat  (const matrix<T>& m1, const matrix<T>& m2);  // horizontal concatenation
template<class T> matrix<T>  vcat  (const matrix<T>& m1, const matrix<T>& m2);  // vertical concatenation
template<class T> void       sort_columns_by_value(matrix<T>& m, const matrix<T>& v);

// 3D geometry
template<class T> matrix<T> outer(const matrix<T>& v, const matrix<T>& w);  // V x W
template<class T> matrix<T> vec2rot (const matrix<T>& v);  // R^3=so(3)-->SO(3) (exp, rotation about v by angle |v|)
template<class T> matrix<T> rot2vec (const matrix<T>& r);  // SO(3)-->so(3)=R^3 (log, local inverse of the above)
template<class T> matrix<T> vec2asym(const matrix<T>& v);  // R^3-->so(3) (AX=VxX)
template<class T> matrix<T> asym2vec(const matrix<T>& a);  // so(3)-->R^3 (AX=VxX)
template<class T> matrix<T> rotabout(int k, T th);         // rotation about k-th coord axis by angle th[rad]

// Complex/real matrix conversion
template<class T> matrix<T> real(const matrix<std::complex<T> >& m);
template<class T> matrix<T> imag(const matrix<std::complex<T> >& m);
template<class T> matrix<std::complex<T> > conj(const matrix<std::complex<T> >& m);
template<class T> matrix<T> real(const matrix<T>& m);
template<class T> matrix<T> imag(const matrix<T>& m);
template<class T> matrix<T> conj(const matrix<T>& m);

// stream I/O
template<class T> std::istream& operator>>(std::istream& str,       matrix<T>& m);  // str>>M
template<class T> std::ostream& operator<<(std::ostream& str, const matrix<T>& m);  // str<<M
template<class T> class omstream;  // m << r,c,x_11,x_12,...,x_rc;


/******** Implementation ********/

// Private utils for complex/real number
template<class T> std::complex<T> nm_sign(std::complex<T> x){ return std::polar(T(1),std::arg(x)); }
template<class T>              T  nm_sign(             T  x){ return x>0 ? +T(1) : -T(1); }

// Class methods
template<class T>        matrix<T>::matrix(void): row(0), col(0), val() {}
template<class T>        matrix<T>::matrix(size_t r, size_t c): row(r), col(c), val(r*c) {}
template<class T>        matrix<T>::matrix(const std::vector<T>& v): row(v.size()), col(1), val(v) {}
template<class T>        matrix<T>::matrix(const T *v0, const T *vn): row(vn-v0), col(1), val(v0,vn) {}
template<class T>        matrix<T>::matrix(T x, T y, T z): row(3), col(1), val(3) {  val[0]=x;  val[1]=y;  val[2]=z;  }
template<class T> size_t matrix<T>::nrow(void) const {  return row;  }
template<class T> size_t matrix<T>::ncol(void) const {  return col;  }
template<class T> size_t matrix<T>::dim (void) const {  return row*col;  }
template<class T> matrix<T>& matrix<T>::resize(size_t r, size_t c) {  row=r; col=c; val.resize(r*c); return *this;  }
template<class T> matrix<T>& matrix<T>::fill  (T z)                {  std::fill(val.begin(),val.end(),z); return *this;  }
template<class T>        matrix<T>::operator const std::vector  <T>& (void) const {  return val;  }
template<class T>        matrix<T>::operator matrix<std::complex<T> >(void) const {  matrix<std::complex<T> > mc(nrow(),ncol()); for(size_t k=0; k<dim(); k++) mc(k)=val[k]; return mc;  }
template<class T> T      matrix<T>::operator()(size_t i, size_t j) const {  return val[i*col+j];  }
template<class T> T      matrix<T>::operator()(size_t i          ) const {  return val[i];        }
template<class T> T&     matrix<T>::operator()(size_t i, size_t j)       {  return val[i*col+j];  }
template<class T> T&     matrix<T>::operator()(size_t i          )       {  return val[i];        }

// Basic operations (incremental)
template<class T> matrix<T>& operator+=(matrix<T>& m, T t){
  if(m.nrow()!=m.ncol()) throw std::domain_error("operator+=(M,t): not square");
  for(size_t i=0; i<m.nrow(); i++) m(i,i)+=t;
  return m;
}
template<class T> matrix<T>& operator-=(matrix<T>& m, T t){  m+=(-t);  return m;  }
template<class T> matrix<T>& operator*=(matrix<T>& m, T t){  for(size_t k=0; k<m.dim(); k++) m(k)*=t;  return m;  }
template<class T> matrix<T>& operator/=(matrix<T>& m, T t){  m*=(T(1)/t);  return m;  }
template<class T> matrix<T>& operator+=(matrix<T>& m, const matrix<T>& dm){
  if(!(m.nrow()==dm.nrow() && m.ncol()==dm.ncol())) throw std::domain_error("operator+=(M,M): sizes mismatch");
  for(size_t k=0; k<m.dim(); k++) m(k)+=dm(k);
  return m;
}
template<class T> matrix<T>& operator-=(matrix<T>& m, const matrix<T>& dm){
  if(!(m.nrow()==dm.nrow() && m.ncol()==dm.ncol())) throw std::domain_error("operator-=(M,M): sizes mismatch");
  for(size_t k=0; k<m.dim(); k++) m(k)-=dm(k);
  return m;
}

// Basic operations (scalar vs. matrix)
template<class T> matrix<T> operator-(const matrix<T>& m0     ){  matrix<T> m=m0;  m*=-T(1); return m;  }
template<class T> matrix<T> operator+(const matrix<T>& m0, T t){  matrix<T> m=m0;  m+=t;  return m;  }
template<class T> matrix<T> operator-(const matrix<T>& m0, T t){  matrix<T> m=m0;  m-=t;  return m;  }
template<class T> matrix<T> operator*(const matrix<T>& m0, T t){  matrix<T> m=m0;  m*=t;  return m;  }
template<class T> matrix<T> operator/(const matrix<T>& m0, T t){  matrix<T> m=m0;  m/=t;  return m;  }
template<class T> matrix<T> operator+(T t, const matrix<T>& m0){  matrix<T> m=m0;  m+=t;  return m;  }
template<class T> matrix<T> operator-(T t, const matrix<T>& m0){  matrix<T> m=-m0; m+=t;  return m;  }
template<class T> matrix<T> operator*(T t, const matrix<T>& m0){  matrix<T> m=m0;  m*=t;  return m;  }
template<class T> matrix<T> operator/(T t, const matrix<T>& m0){  matrix<T> m=inv(m0);  m*=t;  return m;  }

// Basic operations (matrix vs. matrix)
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


// Norm |M|
template<class T> T norm (const matrix<T>& m){
  return sqrt(inner(m,m));
}
template<class T> T norm (const matrix<std::complex<T> >& m){
  T s=0;
  for(size_t k=0; k<m.dim(); k++) s+=std::norm(m(k));  // note: std::norm(C):=|C|^2
  return sqrt(s);
}


// Inner product <M1,M2> (aka "dot" product)
template<class T> T inner(const matrix<T>& m1, const matrix<T>& m2){
  if( !(m1.nrow()==m2.nrow() && m1.ncol()==m2.ncol()) ) throw std::domain_error("inner(M,M): sizes mismatch");
  T t=0;
  for(size_t k=0; k<m1.dim(); k++) t+=std::conj(m1(k))*m2(k);
  return t;
}


// Inverse M^-1
template<class T> matrix<T> inv(const matrix<T>& m){
  if( m.nrow()!=m.ncol() ) throw std::domain_error("inv(M): not square");
  matrix<T> m1=m;
  matrix<T> m2=matrix<T>(m.nrow(),m.nrow())+T(1);  // container of m^-1

  for(size_t k=0; k<m.nrow(); k++){
    // swapping: k-th row <--> i0-th row
    size_t i0=k;
    for(size_t i=k+1; i<m1.nrow(); i++)
      if(std::abs(m1(i,k))>std::abs(m1(i0,k))) i0=i;  // i0=argmax_{i>k} |m1(i,k)|
    if(i0!=k){
      for(size_t j=k; j<m1.nrow(); j++) std::swap(m1(k,j),m1(i0,j));
      for(size_t j=0; j<m2.nrow(); j++) std::swap(m2(k,j),m2(i0,j));
    }

    // elimination: m1(i,k) --> \delta(i,k)
    T t;
    t=m1(k,k);
    for(size_t j=k; j<m1.ncol(); j++) m1(k,j)/=t;  // (division-by-zero handling?)
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


// Power M^n (n<=0 is OK)
template<class T> matrix<T> pow(const matrix<T>& m, int n){
  if(m.nrow()!=m.ncol()) throw std::domain_error("pow(M,n): not square");
  if(n<0) return inv(pow(m,-n));
  matrix<T> m1=m, mn=matrix<T>(m.nrow(),m.ncol())+T(1);
  for( ; n>0; n/=2,m1=m1*m1)
    if(n%2==1) mn=mn*m1;
  return mn;
}


// Transpose M^T
template<class T> matrix<T> tp(const matrix<T>& m0){
  matrix<T> m1(m0.ncol(),m0.nrow());
  for(size_t i=0; i<m0.nrow(); i++)
    for(size_t j=0; j<m0.ncol(); j++)
      m1(j,i)=std::conj(m0(i,j));
  return m1;
}


// Gram-Schmidt orthonormalization of M (an orthonormal frame U s.t. U^T M = upper trianglar)
template<class T> matrix<T> orth(const matrix<T>& m0){
  if(m0.nrow()<m0.ncol()) throw std::domain_error("orth(M): nrow<ncol");
  matrix<T> m=m0, v, v1;
  for(size_t j=0; j<m.ncol(); j++){
    v=getvec(m,j);
    for(size_t j1=0; j1<j; j1++){
      v1=getvec(m,j1);
      v-=inner(v1,v)*v1;
    }
    v/=T(norm(v));  // (division-by-zero handling?)
    setvec(m,j,v);
  }
  return m;
}


// Eigen vectors of symmetric M (an orthonormal frame U s.t. U^T M U = diagonal)
template<class T> matrix<T> eigen(const matrix<T>& m0, double tol, int iter){
  if(m0.nrow()!=m0.ncol()) throw std::domain_error("eigen(M): not square");
  size_t n=m0.nrow();
  matrix<T> m=m0, u=matrix<T>(n,n)+T(1);

  while(iter--){
    // (i0,j0) = argmax_(i,j) |M(i,j)|
    size_t i0=0, j0=1;
    for(size_t i=0; i<n; i++)
      for(size_t j=0; j<n; j++)
	if( i!=j && std::abs(m(i0,j0))<=std::abs(m(i,j)) ) { i0=i; j0=j; }

    // error before rotation - will be decreased by 2|M(i0,j0)|^2 after rotation
    double err0=0;
    for(size_t k=0; k<n; k++)
      err0 += (k==i0 ? 0 : std::norm(m(k,i0))+std::norm(m(i0,k))) + (k==j0 ? 0 :std::norm(m(k,j0))+std::norm(m(j0,k)));

    // Givens rotation R - atan2(2mij,mii-mjj)/2
    T c,s,p,q;
    p=m(i0,i0)-m(j0,j0);
    q=T(2)*m(i0,j0);
    c = cos(atan2(std::abs(q),std::abs(p))/2);
    s = sin(atan2(std::abs(q),std::abs(p))/2) * nm_sign(std::conj(p)) * nm_sign(q);
    for(size_t i=0; i<n; i++){ p=u(i,i0); q=u(i,j0); u(i,i0)=c*p+std::conj(s)*q; u(i,j0)=        -s *p+c*q; }  // U=U*R
    for(size_t i=0; i<n; i++){ p=m(i,i0); q=m(i,j0); m(i,i0)=c*p+std::conj(s)*q; m(i,j0)=        -s *p+c*q; }  // M=M*R
    for(size_t j=0; j<n; j++){ p=m(i0,j); q=m(j0,j); m(i0,j)=c*p+        s *q; m(j0,j)=-std::conj(s)*p+c*q; }  // M=R^T*M

    // error after rotation - check for convergence
    double err1=0;
    for(size_t k=0; k<n; k++)
      err1 += (k==i0 ? 0 : std::norm(m(k,i0))+std::norm(m(i0,k))) + (k==j0 ? 0 :std::norm(m(k,j0))+std::norm(m(j0,k)));
    if(!(tol<err1 && err1<err0)) break;
  }

  return u;
}


// Get/set column vector
template<class T> matrix<T> getvec(const matrix<T>& m, size_t j){
  if(!(j<m.ncol())) throw std::domain_error("getvec(M,j): bad column");
  matrix<T> v(m.nrow());
  for(size_t i=0; i<v.nrow(); i++) v(i)=m(i,j);
  return v;
}
template<class T> void setvec(matrix<T>& m, size_t j, const matrix<T>& v){
  if(!(j<m.ncol() && m.nrow()==v.nrow() && v.ncol()==1)) throw std::domain_error("setvec(M,j,V): bad column or size");
  for(size_t i=0; i<m.nrow(); i++) m(i,j)=v(i);
}


// Get/set submatrix  M i0+[0,r) x j0+[0,c) <--> M1 [0,r) x [0,c)
template<class T> matrix<T>  getsub(const matrix<T>& m, size_t i0, size_t j0, size_t r, size_t c){
  if(m.nrow()<i0+r || m.ncol()<j0+c) throw std::domain_error("getsub(M,i0,j0,r,c): bad dimension");
  matrix<T> m1(r,c);
  for(size_t i=0; i<r; i++)
    for(size_t j=0; j<c; j++)
      m1(i,j)=m(i0+i,j0+j);
  return m1;
}
template<class T> void setsub(matrix<T>& m, size_t i0, size_t j0, const matrix<T> m1){
  size_t r=m1.nrow(), c=m1.ncol();
  if(m.nrow()<i0+r || m.ncol()<j0+c) throw std::domain_error("setsub(M,i0,j0,m1): bad dimension");
  for(size_t i=0; i<r; i++)
    for(size_t j=0; j<c; j++)
      m(i0+i,j0+j)=m1(i,j);
}


// Get/set diagonal part
template<class T> matrix<T> getdiag(const matrix<T>& m){
  size_t n=m.nrow();
  if(m.nrow()!=n || m.ncol()!=n) throw std::domain_error("getdiag(M): not square");
  matrix<T> v(n);
  for(size_t k=0; k<n; k++) v(k)=m(k,k);
  return v;
}
template<class T> void setdiag(matrix<T>& m, const matrix<T>& v){
  size_t n=v.dim();
  if(m.nrow()!=n || m.ncol()!=n) throw std::domain_error("setdiag(M,V): sizes mismatch");
  for(size_t k=0; k<n; k++) m(k,k)=v(k);
}


// Horizonta/vertical concatenation
template<class T> matrix<T> hcat(const matrix<T>& m1, const matrix<T>& m2){
  if(m1.nrow()!=m2.nrow()) throw std::domain_error("hcat(M,M): sizes mismatch");
  matrix<T> m(m1.nrow(),m1.ncol()+m2.ncol());  // M(r,c1+c2)
  for(size_t i=0; i<m.nrow(); i++)
    for(size_t j=0; j<m.ncol(); j++)
      m(i,j) = (j<m1.ncol() ? m1(i,j) : m2(i,j-m1.ncol()));
  return m;
}
template<class T> matrix<T> vcat(const matrix<T>& m1, const matrix<T>& m2){
  if(m1.ncol()!=m2.ncol()) throw std::domain_error("vcat(M,M): sizes mismatch");
  matrix<T> m(m1.nrow()+m2.nrow(),m1.ncol());  // M(r1+r2,c)
  for(size_t i=0; i<m.nrow(); i++)
    for(size_t j=0; j<m.ncol(); j++)
      m(i,j) = (i<m1.nrow() ? m1(i,j) : m2(i-m1.nrow(),j));
  return m;
}


template<class T> void sort_columns_by_value(matrix<T>& m, const matrix<T>& v){
  std::vector<size_t> kk(v.dim());
  std::iota(kk.begin(), kk.end(), 0);
  std::sort(kk.begin(), kk.end(), [&](int i, int j){ return std::real(v(i))<std::real(v(j)); });
  matrix<T> m1=m;
  for(size_t k=0; k<v.dim(); k++)
    if(k<m.ncol() && kk[k]<m1.ncol() && k!=kk[k]) setvec(m,k,getvec(m1,kk[k]));
}


// 3D geometry
template<class T> matrix<T> outer(const matrix<T>& v, const matrix<T>& w){  // V x W
  if( !(v.nrow()==3 && w.ncol()==1 && v.nrow()==3 && v.ncol()==1) ) throw std::domain_error("outer(V,V): not 3-vectors");
  return matrix<T>(v(1)*w(2)-w(1)*v(2), v(2)*w(0)-w(2)*v(0), v(0)*w(1)-w(0)*v(1));
}
template<class T> matrix<T> vec2rot(const matrix<T>& v){  // R^3=so(3)-->SO(3) (exp, rotation about v by angle |v|)
  if(!(v.nrow()==3 && v.ncol()==1)) throw std::domain_error("vec2rot(V): not a 3-vector");
  double th=norm(v);
  if(std::abs(th)<1.e-8){ matrix<T> r=vec2asym(v);  return r*r/T(2)+r+T(1); }
  else{ matrix<T> r=vec2asym(v/th); return ((T(1)-cos(th))*r+sin(th))*r+T(1); }  // Rodrigues
}
template<class T> matrix<T> rot2vec(const matrix<T>& r){  // SO(3)-->so(3)=R^3 (log, local inverse of exp)
  if(!(r.nrow()==3 && r.ncol()==3)) throw std::domain_error("rot2vec(R): not a 3x3 matrix");
  T c=(r(0,0)+r(1,1)+r(2,2)-1)/2;
  T th=(c<=-1 ? M_PI : 1<=c ? 0 : acos(c));  // (singularities sin(th)=0 will be treated separately)
  if(c<-0.99){
    matrix<T> sym=(r+tp(r))/2.0-1.0, axis(3);
    for(int k=0; k<3; k++) axis+=outer(getvec(sym,(k+1)%3), getvec(sym,(k+2)%3));
    return (th/norm(axis))*axis;
  }
  matrix<T> axis((r(2,1)-r(1,2))/2, (r(0,2)-r(2,0))/2, (r(1,0)-r(0,1))/2);
  if(th<1.e-6) return (1+inner(axis,axis)/6)*axis;
  else         return (th/sin(th))          *axis;
}
template<class T> matrix<T> asym2vec(const matrix<T>& a){  // so(3)-->R^3 (A-->V s.t. AX=VxX)
  if(!(a.nrow()==3 && a.ncol()==3)) throw std::domain_error("asym2vec(V): not a 3 by 3 matrix");
  return matrix<T>(a(2,1),a(0,2),a(1,0));
}
template<class T> matrix<T> vec2asym(const matrix<T>& v){  // R^3-->so(3) (V-->A s.t. AX=VxX)
  if(!(v.nrow()==3 && v.ncol()==1)) throw std::domain_error("vec2asym(V): not a 3-vector");
  matrix<T> a(3,3);
  a(2,1)=v(0);  a(1,2)=-v(0);
  a(0,2)=v(1);  a(2,0)=-v(1);
  a(1,0)=v(2);  a(0,1)=-v(2);
  return a;
}
template<class T> matrix<T> rotabout(int k, T th){  // rotation about k-th coord axis by angle th
  if(!(0<=k && k<3)) throw std::domain_error("rotabout(k,th): bad axis");
  matrix<T> r(3,3);
  int i=(k+1)%3, j=(k+2)%3;
  r(i,i)=cos(th);  r(i,j)=-sin(th);
  r(j,i)=sin(th);  r(j,j)= cos(th);
  r(k,k)=T(1);
  return r;
}


// Complex/real matrix conversion
template<class T> matrix<T> real(const matrix<std::complex<T> >& m){
  matrix<T> mr(m.nrow(),m.ncol());
  for(size_t k=0; k<m.dim(); k++) mr(k)=m(k).real();
  return mr;
}
template<class T> matrix<T> imag(const matrix<std::complex<T> >& m){
  matrix<T> mr(m.nrow(),m.ncol());
  for(size_t k=0; k<m.dim(); k++) mr(k)=m(k).imag();
  return mr;
}
template<class T> matrix<std::complex<T> > conj(const matrix<std::complex<T> >& m){
  matrix<std::complex<T> > mc(m.nrow(),m.ncol());
  for(size_t k=0; k<m.dim(); k++) mc(k)=std::conj(m(k));
  return mc;
}
template<class T> matrix<T> real(const matrix<T>& m){ return m; }
template<class T> matrix<T> imag(const matrix<T>& m){ return matrix<T>(m.nrow(),m.ncol()); }
template<class T> matrix<T> conj(const matrix<T>& m){ return m; }


// stream I/O (format: r c m00 m01... m10 m11...)
template<class T> std::istream& operator>>(std::istream& str, matrix<T>& m){
  size_t r,c;
  str >> r >> c;
  if(!str) return str;
  m=matrix<T>(r,c);
  for(size_t k=0; k<r*c; k++) str >> m(k);
  return str;
}
template<class T> std::ostream& operator<<(std::ostream& str, const matrix<T>& m){
  size_t r=m.nrow(), c=m.ncol();
  str << r << ' ' << c << '\t';
  if(!str) return str;
  for(size_t k=0; k<r*c; k++) str << m(k) << ((k+1)%c ? ' ' : '\t');
  return str;
}
template<class T> std::istream& read(std::istream& str, matrix<T>& m){
  size_t r,c;
  str.read((char*)&r,sizeof(r));
  str.read((char*)&c,sizeof(c));
  if(!str) return str;
  m=matrix<T>(r,c);
  str.read((char*)&m(0),sizeof(m(0))*r*c);  // note: assumes contiguous memory
  return str;
}
template<class T> std::ostream& write(std::ostream& str, const matrix<T>& m){
  size_t r=m.nrow(),c=m.ncol();
  str.write((char*)&r,sizeof(r));
  str.write((char*)&c,sizeof(c));
  if(!str) return str;
  const std::vector<T>& val=m;  // note: &m(0) does not work
  str.write((char*)&val[0],sizeof(val[0])*r*c);
  return str;
}

// output matrix stream (usage: m << r,c,x_11,x_12, ..., x_rc)
template<class T> omstream<T> operator<<(matrix<T>& m0, size_t r0){ return omstream<T>(m0,r0); }  // m<<r...
template<class T> class omstream{
public:
  omstream(matrix<T>& m0, size_t r0)
    : m(m0),r(r0),c(0),k(0) {}
  omstream& operator,(T x){
    if(c==0){ c=size_t(std::real<T>(x)+0.5); m.resize(r,c); return *this; }  // m<<r,c...
    if(k<r*c){ m(k++)=x; return *this; }  // m<<r,c,x_11...
    throw std::runtime_error("omstream: too many arguments");
  }
  matrix<T>& m;
  size_t r,c,k;
};


}  //namespace nmlib
#endif //MATRIX_H
