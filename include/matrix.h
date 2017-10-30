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
  matrix(void);                               // empty matrix
  explicit matrix(size_t r, size_t c=1);      // zero matrix of size (r,c)
  matrix(const std::vector<T>& v);            // column vector
  matrix(const std::initializer_list<T>& v);  // column vector M={x,y,z..}
  matrix(const std::initializer_list<std::initializer_list<T>>& v);  // column vectors M={{x1,y1,z1..},{x2,y2,z2..}..}
  template<class It> matrix         (size_t r, size_t c, const It& first, const It& last);
  template<class It> matrix<T>& init(size_t r, size_t c, const It& first, const It& last);

  bool   good(void) const;  // health check
  size_t nrow(void) const;  // number of rows
  size_t ncol(void) const;  // number of columns
  size_t dim (void) const;  // number of components
  matrix<T>& resize(size_t r, size_t c);
  matrix<T>& fill  (T z);
  matrix<T>& swap(matrix<T>& m);

  T  operator()(size_t i, size_t j) const;  // get t=M(i,j)
  T& operator()(size_t i, size_t j);        // set M(i,j)=t
  T  operator()(size_t i) const;  // get t=V(i)
  T& operator()(size_t i);        // set V(i)=t
  const std::vector<T>& as_vector(void) const;

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

// Element-wise operations
template<class T,class OP> matrix<T>& opeq(      matrix<T>& m0,const OP& op1);                       // (op1(Mij))-->M
template<class T,class OP> matrix<T>& opeq(      matrix<T>& m1,const OP& op2, const T& t );          // (op2(Mij,t))-->M
template<class T,class OP> matrix<T>& opeq(      matrix<T>& m1,const OP& op2, const matrix<T>& m2);  // (op2(Mij,Mij))-->M
template<class T,class OP> matrix<T>  op  (const matrix<T>& m0,const OP& op1);                       // (op1(Mij))
template<class T,class OP> matrix<T>  op  (const matrix<T>& m1,const OP& op2, const T& t );          // (op2(Mij,t))
template<class T,class OP> matrix<T>  op  (const matrix<T>& m1,const OP& op2, const matrix<T>& m2);  // (op2(Mij,Mij))

// Miscelanaous operations
template<class T> T          norm  (const matrix<T>& m);                        // norm |M|
template<class T> T          norm  (const matrix<std::complex<T> >& m);         // norm |M| (complex M)
template<class T> T          inner (const matrix<T>& m1, const matrix<T>& m2);  // inner product <M1,M2>
template<class T> T          det   (const matrix<T>& m);                        // determinant
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
template<class T> matrix<T>  solve (const matrix<T>& m1, const matrix<T>& m2);  // M1^-1 M2
template<class T> int        gaussian_elimination(matrix<T>& a, matrix<T>& b, bool norm_d, bool elim_u);
template<class T> void       sort_columns_by_value(matrix<T>& m, const matrix<T>& v);

// 3D geometry
template<class T> matrix<T> outer(const matrix<T>& v, const matrix<T>& w);  // V x W
template<class T> matrix<T> vec2rot (const matrix<T>& v);  // R^3=so(3)-->SO(3) (exp, rotation about v by angle |v|)
template<class T> matrix<T> rot2vec (const matrix<T>& r);  // SO(3)-->so(3)=R^3 (log, local inverse of the above)
template<class T> matrix<T> vec2asym(const matrix<T>& v);  // R^3-->so(3) (AX=VxX)
template<class T> matrix<T> asym2vec(const matrix<T>& a);  // so(3)-->R^3 (AX=VxX)
template<class T> matrix<T> rotabout(int k, T th);         // rotation about k-th coord axis by angle th[rad]

// Complex/real matrix conversion
template<class T> matrix<std::complex<T>> complex(const matrix<T>& m);
template<class T> matrix<T> real(const matrix<std::complex<T> >& m);
template<class T> matrix<T> imag(const matrix<std::complex<T> >& m);
template<class T> matrix<std::complex<T> > conj(const matrix<std::complex<T> >& m);
template<class T> matrix<T> real(const matrix<T>& m);
template<class T> matrix<T> imag(const matrix<T>& m);
template<class T> matrix<T> conj(const matrix<T>& m);

// stream I/O
template<class T> std::istream& operator>>(std::istream& str,       matrix<T>& m);  // str>>M
template<class T> std::ostream& operator<<(std::ostream& str, const matrix<T>& m);  // str<<M
template<class T> std::istream& read (std::istream& str,       matrix<T>& m);  // binary I/O
template<class T> std::ostream& write(std::ostream& str, const matrix<T>& m);


/******** Implementation ********/

// Private utils for complex/real number
template<class T> std::complex<T> nm_conj(const std::complex<T>& x){ return std::conj(x); }
template<class T>              T  nm_conj(                   T   x){ return           x;  }
template<class T> std::complex<T> nm_sign(const std::complex<T>& x){ return std::polar(T(1),std::arg(x)); }
template<class T>              T  nm_sign(                   T   x){ return x>0 ? +T(1) : -T(1); }

// Class methods
template<class T>        matrix<T>::matrix(void): row(0), col(0), val() {}
template<class T>        matrix<T>::matrix(size_t r, size_t c): row(r), col(c), val(r*c) {}
template<class T>        matrix<T>::matrix(const std::vector<T>& v): row(v.size()), col(1), val(v) {}
template<class T>        matrix<T>::matrix(const std::initializer_list<T>& v): row(v.size()), col(1), val(v) {}
template<class T>        matrix<T>::matrix(const std::initializer_list<std::initializer_list<T>>& v): row(0), col(0), val() {
  if( v.size()==0 ) return;
  col=v.size();  row=v.begin()->size();  val.resize(row*col);
  size_t j=0;
  for(auto c: v){
    if( c.size()!=row ) throw std::domain_error("matrix(): bad initializer_list (variable row size)");
    size_t i=0;
    for(auto x: c){ val[i*col+j]=x; i++; }
    j++;
  }
}
template<class T> template<class It> matrix<T>::matrix(size_t r, size_t c, const It& first, const It& last)
  : row(r),col(c),val(first,last) {
  if(row*col!=val.size()) throw std::domain_error("matrix: init: illegal size");
}
template<class T> template<class It> matrix<T>& matrix<T>::init(size_t r, size_t c, const It& first, const It& last){
  matrix<T> m(r,c,first,last);
  return this->swap(m);
}
template<class T> bool   matrix<T>::good(void) const {
  bool ok = (row*col==val.size());
  for(auto x: val)
    ok &= std::isfinite(std::real(x)) && std::isfinite(std::imag(x));
  return ok;
}
template<class T> size_t matrix<T>::nrow(void) const {  return row;  }
template<class T> size_t matrix<T>::ncol(void) const {  return col;  }
template<class T> size_t matrix<T>::dim (void) const {  return row*col;  }
template<class T> matrix<T>& matrix<T>::resize(size_t r, size_t c) {  row=r; col=c; val.resize(r*c); return *this;  }
template<class T> matrix<T>& matrix<T>::fill  (T z)                {  std::fill(val.begin(),val.end(),z); return *this;  }
template<class T> matrix<T>& matrix<T>::swap  (matrix<T>& m)       {  std::swap(row,m.row); std::swap(col,m.col); val.swap(m.val); return *this;  }
template<class T> T      matrix<T>::operator()(size_t i, size_t j) const {  return val[i*col+j];  }
template<class T> T      matrix<T>::operator()(size_t i          ) const {  return val[i];        }
template<class T> T&     matrix<T>::operator()(size_t i, size_t j)       {  return val[i*col+j];  }
template<class T> T&     matrix<T>::operator()(size_t i          )       {  return val[i];        }
template<class T> const std::vector<T>& matrix<T>::as_vector(void) const {  return val;  }

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


// Element-wise operations
template<class T,class OP> matrix<T>& opeq(matrix<T>& m0,const OP& op1){
  for(size_t k=0; k<m0.dim(); k++) m0(k)=op1(m0(k));
  return m0;
}
template<class T,class OP> matrix<T>& opeq(matrix<T>& m0,const OP& op2, const T& t){
  for(size_t k=0; k<m0.dim(); k++) m0(k)=op2(m0(k),t);
  return m0;
}
template<class T,class OP> matrix<T>& opeq(matrix<T>& m1,const OP& op2, const matrix<T>& m2){
  if( !(m1.nrow()==m2.nrow() && m1.ncol()==m2.ncol()) )
    throw std::domain_error("op(M,M): sizes mismatch");
  for(size_t k=0; k<m1.dim(); k++) m1(k)=op2(m1(k),m2(k));
  return m1;
}
template<class T,class OP> matrix<T> op(const matrix<T>& m0,const OP& op1){
  matrix<T> m(m0);
  return opeq(m,op1);
}
template<class T,class OP> matrix<T> op(const matrix<T>& m0,const OP& op2, const T& t){
  matrix<T> m(m0);
  return opeq(m,op2,t);
}
template<class T,class OP> matrix<T> op(const matrix<T>& m1,const OP& op2, const matrix<T>& m2){
  matrix<T> m(m1);
  return opeq(m,op2,m2);
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
  for(size_t k=0; k<m1.dim(); k++) t+=nm_conj(m1(k))*m2(k);
  return t;
}


// Determinant
template<class T> T det(const matrix<T>& m){
  const size_t n=m.nrow();
  if( m.nrow()!=m.ncol() ) throw std::domain_error("det(M): not square");
  if( n==2 ) return m(0,0)*m(1,1)-m(1,0)*m(0,1);
  if( n==3 ) return
	       m(0,0)*(m(1,1)*m(2,2)-m(2,1)*m(1,2)) +
	       m(1,0)*(m(2,1)*m(0,2)-m(0,1)*m(2,2)) +
	       m(2,0)*(m(0,1)*m(1,2)-m(1,1)*m(0,2));
  matrix<T> m1(m),dummy(n,0);
  T t=gaussian_elimination(m1,dummy,false,false);
  for(size_t k=0; k<n; k++) t*=m1(k,k);
  return t;
}
// Inverse M^-1
template<class T> matrix<T> inv(const matrix<T>& m){
  const size_t n=m.nrow();
  if( m.nrow()!=m.ncol() ) throw std::domain_error("inv(M): not square");
  matrix<T> m1(m), m2(n,n);
  m2+=T(1);
  gaussian_elimination(m1,m2,true,true);
  return m2;
}
// Solve M1 X = M2
template<class T> matrix<T> solve(const matrix<T>& m1, const matrix<T>& m2){
  const size_t n=m1.nrow();
  if( !(m1.nrow()==m1.ncol() && m2.nrow()==n) ) throw std::domain_error("solve(): not square");
  matrix<T> a(m1), x(m2);
  gaussian_elimination(a,x,true,true);
  return x;
}


// Gaussian elimination (solve AX=B by successive fundamentarl row operations)
template<class T> int gaussian_elimination(matrix<T>& a, matrix<T>& b, bool norm_d, bool elim_u){
  const size_t na=a.ncol(), nb=b.ncol();
  if( !(a.nrow()==a.ncol() && b.nrow()==na) ) throw std::domain_error("gaussian_elimination(): sizes mismatch");

  int sgn=+1;  // sign of permutation of rows

  for(size_t k=0; k<na; k++){
    // pivoting: k-th <--> i0-th row
    size_t i0=k;
    for(size_t i=k+1; i<na; i++)
      if(std::abs(a(i,k))>std::abs(a(i0,k))) i0=i;  // i0=argmax_{i>k} |a(i,k)|
    if(i0!=k){
      for(size_t j=k; j<na; j++) std::swap(a(k,j),a(i0,j));
      for(size_t j=0; j<nb; j++) std::swap(b(k,j),b(i0,j));
      sgn=-sgn;
    }

    // normalization: a(k,k) --> 1
    T t0=a(k,k);
    if( norm_d ){
      T t=T(1)/t0;  // (division-by-zero handling?)
      for(size_t j=k+1; j<na; j++) a(k,j)*=t;  a(k,k)=1;
      for(size_t j=0;   j<nb; j++) b(k,j)*=t;
    }

    // elimination: a(i,k) --> 0 for i!=k
    for(size_t i=0; i<na; i++){
      if( i==k ) continue;
      T t=a(i,k);
      if( (i<k && !elim_u)  ||  t==T(0) ) continue;
      if( !norm_d ) t/=t0;
      for(size_t j=k+1; j<na; j++) a(i,j)-=a(k,j)*t;  a(i,k)=0;
      for(size_t j=0;   j<nb; j++) b(i,j)-=b(k,j)*t;
    }
  }

  return sgn;
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
      m1(j,i)=nm_conj(m0(i,j));
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
    s = sin(atan2(std::abs(q),std::abs(p))/2) * nm_sign(nm_conj(p)) * nm_sign(q);
    for(size_t i=0; i<n; i++){ p=u(i,i0); q=u(i,j0); u(i,i0)=c*p+nm_conj(s)*q; u(i,j0)=        -s *p+c*q; }  // U=U*R
    for(size_t i=0; i<n; i++){ p=m(i,i0); q=m(i,j0); m(i,i0)=c*p+nm_conj(s)*q; m(i,j0)=        -s *p+c*q; }  // M=M*R
    for(size_t j=0; j<n; j++){ p=m(i0,j); q=m(j0,j); m(i0,j)=c*p+        s *q; m(j0,j)=-nm_conj(s)*p+c*q; }  // M=R^T*M

    // error after rotation - check for convergence
    double err1=0, eps=std::numeric_limits<double>::epsilon();
    for(size_t k=0; k<n; k++)
      err1 += (k==i0 ? 0 : std::norm(m(k,i0))+std::norm(m(i0,k))) + (k==j0 ? 0 :std::norm(m(k,j0))+std::norm(m(j0,k)));
    if(!(tol<err1 && err1<err0)) break;
    if( err0-err1 < (err0+err1)*eps ) break;  // workaround for infinite loop
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
  return matrix<T>({v(1)*w(2)-w(1)*v(2), v(2)*w(0)-w(2)*v(0), v(0)*w(1)-w(0)*v(1)});
}
template<class T> matrix<T> vec2rot(const matrix<T>& v){  // R^3=so(3)-->SO(3) (exp, rotation about v by angle |v|)
  if(!(v.nrow()==3 && v.ncol()==1)) throw std::domain_error("vec2rot(V): not a 3-vector");
  T th=norm(v), c=cos(th), s=sin(th);
  if( th<1.e-6 ) return c + vec2asym(v) + 0.5*v*tp(v);  // small angle
  matrix<T> u=v/th;
  return c + s*vec2asym(u) + (1-c)*u*tp(u);  // Rodrigues
}
template<class T> matrix<T> rot2vec(const matrix<T>& r){  // SO(3)-->so(3)=R^3 (log, local inverse of exp)
  if(!(r.nrow()==3 && r.ncol()==3)) throw std::domain_error("rot2vec(R): not a 3x3 matrix");
  if( (r(0,0)+r(1,1)+r(2,2)-1)/2 > 1-0.5e-12 ) return asym2vec(r-tp(r))*0.5;  // cos(th)=(trR-1)/2=1-eps

  matrix<T> m=r-T(1);  // Im(M)=rotation plane, Ker(M)=rotation axis
  matrix<T> a0=getvec(m,0),  a1=getvec(m,1),  a2=getvec(m,2);  // Im
  matrix<T> b0=outer(a1,a2), b1=outer(a2,a0), b2=outer(a0,a1); // Ker
  T         s0=norm(a0),     s1=norm(a1),     s2=norm(a2);  // for stability
  T         t0=norm(b0),     t1=norm(b1),     t2=norm(b2);

  a0 = (s0<=s2 && s1<=s2 ? a2 : s0<=s1 ? a1 : a0);
  b0 = (t0<=t2 && t1<=t2 ? b2 : t0<=t1 ? b1 : b0);
  a1 = r*a0;
  b1 = outer(a0,a1);
  T th = atan2(norm(b1),inner(a0,a1));
  if( inner(b0,b1)<0 ) th=-th;
  return (th/norm(b0))*b0;
}
template<class T> matrix<T> asym2vec(const matrix<T>& a){  // so(3)-->R^3 (A-->V s.t. AX=VxX)
  if(!(a.nrow()==3 && a.ncol()==3)) throw std::domain_error("asym2vec(V): not a 3 by 3 matrix");
  return matrix<T>({a(2,1),a(0,2),a(1,0)});
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
template<class T> matrix<std::complex<T>> complex(const matrix<T>& m){
  matrix<std::complex<T>> mc(m.nrow(),m.ncol());
  for(size_t k=0; k<m.dim(); k++) mc(k)=m(k);
  return mc;
}
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
  for(size_t k=0; k<m.dim(); k++) mc(k)=nm_conj(m(k));
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
  m.resize(r,c);
  str.read((char*)(&m(0)),sizeof(T)*m.dim());  // contiguity of v is assumed
  return str;
}
template<class T> std::ostream& write(std::ostream& str, const matrix<T>& m){
  size_t r=m.nrow(),c=m.ncol();
  str.write((char*)&r,sizeof(r));
  str.write((char*)&c,sizeof(c));
  if(!str) return str;
  const std::vector<T>& v=m.as_vector();
  str.write((char*)(v.data()),sizeof(T)*v.size());  // contiguity of v is assumed
  return str;
}


}  //namespace nmlib
#endif //MATRIX_H
